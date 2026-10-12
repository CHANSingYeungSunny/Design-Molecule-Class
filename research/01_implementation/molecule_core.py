"""An explicit molecular graph and a standard-library-only SMILES codec.

Scope: a connectivity-only SMILES subset; no specified stereochemistry.
The parser/writer preserve supplied aromatic annotations and bracket H policy.
They do NOT perform chemical sanitization, aromaticity perception, radical
assignment, canonicalization, or a general valence check. See README.md.

No external package, subprocess, network service, or stored source SMILES is used.
"""
from __future__ import annotations

import ast
import re
from collections import deque
from dataclasses import dataclass, fields, replace
from typing import Iterable


class MoleculeError(ValueError):
    """Malformed structure or malformed supported SMILES syntax."""


class UnsupportedSmilesError(MoleculeError):
    """Input requires a feature outside the declared codec scope."""


ELEMENTS = frozenset((
    'H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co '
    'Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb '
    'Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re '
    'Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es '
    'Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og *'
).split())
ORGANIC = frozenset(('B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I', '*'))
AROMATIC = frozenset(('B', 'C', 'N', 'O', 'P', 'S', 'As', 'Se'))
BARE_AROMATIC = frozenset('bcnops')
BOND_SYMBOLS = {'-': 'SINGLE', '=': 'DOUBLE', '#': 'TRIPLE', ':': 'AROMATIC'}
MAX_TEXT = 100_000
MAX_ATOMS = 10_000


@dataclass(frozen=True)
class Atom:
    id: str
    element: str
    formal_charge: int = 0
    isotope: int = 0
    aromatic: bool = False
    # None means unbracketed-style implicit-H rules apply.
    # An integer means an exact bracket H count; no additional implicit H.
    # Thus C -> None, [C] -> 0, [CH3] -> 3, and [nH] -> 1.
    explicit_hs: int | None = None
    atom_map: int | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.id, str) or not self.id:
            raise MoleculeError('Atom id must be a nonempty string.')
        if not isinstance(self.element, str) or self.element not in ELEMENTS:
            raise MoleculeError(f'Unknown element: {self.element!r}')
        if type(self.formal_charge) is not int or abs(self.formal_charge) > 100:
            raise MoleculeError('formal_charge must be an integer in [-100, 100].')
        if type(self.isotope) is not int or not 0 <= self.isotope <= 999:
            raise MoleculeError('isotope must be an integer in [0, 999].')
        if type(self.aromatic) is not bool:
            raise MoleculeError('aromatic must be a bool.')
        if self.aromatic and self.element not in AROMATIC:
            raise MoleculeError(f'Aromatic notation is unsupported for {self.element}.')
        if self.explicit_hs is not None:
            if type(self.explicit_hs) is not int or not 0 <= self.explicit_hs <= 9:
                raise MoleculeError('explicit_hs must be None or an integer in [0, 9].')
        if self.element == 'H' and self.explicit_hs not in (None, 0):
            raise MoleculeError('Hydrogen-to-hydrogen bonds require explicit H atom vertices.')
        if self.atom_map is not None:
            if type(self.atom_map) is not int or not 0 <= self.atom_map < 2**31:
                raise MoleculeError('atom_map must be None or a nonnegative 32-bit integer.')
        bare_symbol = self.element.lower() if self.aromatic else self.element
        can_be_bare = (bare_symbol in BARE_AROMATIC if self.aromatic
                       else self.element in ORGANIC)
        needs_bracket = (not can_be_bare or self.formal_charge != 0
                         or self.isotope != 0 or self.atom_map is not None)
        if needs_bracket and self.explicit_hs is None:
            raise MoleculeError(
                'This atom requires brackets. Specify explicit_hs=0, or the '
                'intended bracket H count, to avoid silently changing H semantics.')


@dataclass(frozen=True)
class Bond:
    begin: str
    end: str
    bond_type: str = 'SINGLE'

    def __post_init__(self) -> None:
        if not isinstance(self.begin, str) or not isinstance(self.end, str):
            raise MoleculeError('Bond endpoints must be atom ID strings.')
        if not isinstance(self.bond_type, str) or self.bond_type not in BOND_SYMBOLS.values():
            raise MoleculeError(f'Unsupported bond type: {self.bond_type!r}')


def validate_graph(atoms: tuple[Atom, ...], bonds: tuple[Bond, ...]) -> None:
    """Check record consistency, NOT chemical validity or realizability."""
    if not atoms or len(atoms) > MAX_ATOMS:
        raise MoleculeError(f'The graph must contain 1..{MAX_ATOMS} atoms.')
    if any(type(a) is not Atom for a in atoms) or any(type(b) is not Bond for b in bonds):
        raise MoleculeError('Expected Atom and Bond records.')
    by_id = {a.id: a for a in atoms}
    if len(by_id) != len(atoms):
        raise MoleculeError('Atom IDs must be unique.')
    pairs = set()
    for bond in bonds:
        if bond.begin not in by_id or bond.end not in by_id:
            raise MoleculeError('A bond references an unknown atom.')
        pair = frozenset((bond.begin, bond.end))
        if len(pair) != 2 or pair in pairs:
            raise MoleculeError('Self-bonds and duplicate bonds are not allowed.')
        pairs.add(pair)
        if bond.bond_type == 'AROMATIC' and not (
            by_id[bond.begin].aromatic and by_id[bond.end].aromatic
        ):
            raise MoleculeError('An aromatic bond must join two aromatic atoms.')


_BRACKET = re.compile(
    r'(?P<isotope>[0-9]*)(?P<element>\*|[A-Z][a-z]?|se|as|[bcnops])'
    r'(?P<hs>H[0-9]?)?(?P<charge>\++|-+|[+-][0-9]+)?'
    r'(?::(?P<map>[0-9]+))?'
)


def _read_atom(text: str, pos: int, atom_id: str) -> tuple[Atom, int]:
    if text[pos] == '[':
        end = text.find(']', pos + 1)
        if end < 0:
            raise MoleculeError(f'Unclosed bracket at character {pos}.')
        match = _BRACKET.fullmatch(text[pos + 1:end])
        if match is None:
            raise UnsupportedSmilesError(f'Malformed or unsupported bracket atom at {pos}.')
        symbol = match['element']
        htoken, ctoken = match['hs'], match['charge']
        if symbol == 'H' and htoken is not None:
            raise MoleculeError('A bracket hydrogen atom cannot have its own H-count token.')
        hs = 0 if htoken is None else (int(htoken[1:]) if len(htoken) > 1 else 1)
        charge = 0
        if ctoken:
            magnitude = int(ctoken[1:]) if ctoken[1:].isdigit() else len(ctoken)
            charge = magnitude if ctoken[0] == '+' else -magnitude
        return Atom(
            id=atom_id, element=symbol.capitalize() if symbol != '*' else '*',
            formal_charge=charge, isotope=int(match['isotope'] or 0),
            aromatic=symbol.islower(), explicit_hs=hs,
            atom_map=int(match['map']) if match['map'] is not None else None,
        ), end + 1
    symbol = text[pos:pos + 2] if text[pos:pos + 2] in ('Cl', 'Br') else text[pos]
    if symbol not in ORGANIC and symbol not in BARE_AROMATIC:
        raise UnsupportedSmilesError(f'Unsupported atom/token {symbol!r} at {pos}.')
    return Atom(id=atom_id, element=symbol.upper() if symbol.islower() else symbol,
                aromatic=symbol.islower()), pos + len(symbol)


def parse_smiles(text: str) -> tuple[tuple[Atom, ...], tuple[Bond, ...]]:
    """Parse the declared non-stereo SMILES subset into explicit graph records."""
    if not isinstance(text, str) or not text or len(text) > MAX_TEXT:
        raise MoleculeError(f'SMILES must be a nonempty string of at most {MAX_TEXT} characters.')
    if not text.isascii() or any(c.isspace() for c in text):
        raise UnsupportedSmilesError('Only ASCII SMILES without names or whitespace are supported.')
    if any(c in text for c in ('@', '/', '\\')):
        raise UnsupportedSmilesError('Specified stereochemistry is unsupported; it is never discarded.')
    if any(c in text for c in ('|', '>', '<', '~', '$', '?', '&', ';', '!')):
        raise UnsupportedSmilesError('CXSMILES, reaction/query syntax, and extended bonds are unsupported.')
    atoms: list[Atom] = []
    bonds: list[Bond] = []
    current: int | None = None
    pending: str | None = None
    branches: list[tuple[int, int]] = []
    rings: dict[int, tuple[int, str | None]] = {}
    pairs: set[tuple[int, int]] = set()
    last = 'start'
    atomlike = ('atom', 'ring', 'close')

    def connect(a: int, b: int, kind: str | None) -> None:
        pair = tuple(sorted((a, b)))
        if a == b or pair in pairs:
            raise MoleculeError('SMILES specifies a self-bond or duplicate bond.')
        pairs.add(pair)
        kind = kind or ('AROMATIC' if atoms[a].aromatic and atoms[b].aromatic else 'SINGLE')
        bonds.append(Bond(atoms[a].id, atoms[b].id, kind))

    i = 0
    while i < len(text):
        char = text[i]
        if char == '(':
            if current is None or last not in atomlike:
                raise MoleculeError(f'Unexpected branch opening at {i}.')
            branches.append((current, len(atoms)))
            last = 'open'
            i += 1
        elif char == ')':
            if not branches or last not in atomlike or len(atoms) == branches[-1][1]:
                raise MoleculeError(f'Unexpected or empty branch closing at {i}.')
            current, _ = branches.pop()
            last = 'close'
            i += 1
        elif char == '.':
            if current is None or last not in atomlike:
                raise MoleculeError(f'Unexpected component separator at {i}.')
            if branches:
                raise UnsupportedSmilesError('Component separators inside branches are outside this subset.')
            current, last = None, 'dot'
            i += 1
        elif char in BOND_SYMBOLS:
            if current is None or last not in (*atomlike, 'open'):
                raise MoleculeError(f'Unexpected bond symbol at {i}.')
            pending, last = BOND_SYMBOLS[char], 'bond'
            i += 1
        elif char.isdigit() or char == '%':
            if current is None or last not in (*atomlike, 'bond'):
                raise MoleculeError(f'Unexpected ring label at {i}.')
            if branches and len(atoms) == branches[-1][1]:
                raise MoleculeError('A branch must begin with an atom, not a ring label.')
            if char == '%':
                label = text[i + 1:i + 3]
                if len(label) != 2 or not label.isdigit():
                    raise UnsupportedSmilesError('Percent ring labels require exactly two digits.')
                number, i = int(label), i + 3
            else:
                number, i = int(char), i + 1
            if number in rings:
                other, opening_kind = rings.pop(number)
                if pending and opening_kind and pending != opening_kind:
                    raise MoleculeError('Conflicting bond types on a ring closure.')
                connect(other, current, pending or opening_kind)
            else:
                rings[number] = (current, pending)
            pending, last = None, 'ring'
        else:
            atom, i = _read_atom(text, i, f'a{len(atoms)}')
            atoms.append(atom)
            if len(atoms) > MAX_ATOMS:
                raise UnsupportedSmilesError('Atom count exceeds the configured limit.')
            new = len(atoms) - 1
            if current is not None:
                connect(current, new, pending)
            current, pending, last = new, None, 'atom'
    if last not in atomlike or pending is not None or branches or rings:
        raise MoleculeError('Incomplete SMILES: trailing token, unclosed branch, or unpaired ring label.')
    result = tuple(atoms), tuple(bonds)
    validate_graph(*result)
    return result


def _atom_text(atom: Atom) -> str:
    symbol = atom.element.lower() if atom.aromatic else atom.element
    if atom.explicit_hs is None:
        return symbol
    isotope = str(atom.isotope) if atom.isotope else ''
    hs = '' if atom.explicit_hs == 0 else ('H' if atom.explicit_hs == 1 else f'H{atom.explicit_hs}')
    charge = ''
    if atom.formal_charge:
        charge = '+' if atom.formal_charge > 0 else '-'
        if abs(atom.formal_charge) != 1:
            charge += str(abs(atom.formal_charge))
    mapping = '' if atom.atom_map is None else f':{atom.atom_map}'
    return f'[{isotope}{symbol}{hs}{charge}{mapping}]'


def _bond_text(bond: Bond, by_id: dict[str, Atom]) -> str:
    if bond.bond_type == 'SINGLE':
        # A missing symbol between two aromatic atoms would mean aromatic.
        return '-' if by_id[bond.begin].aromatic and by_id[bond.end].aromatic else ''
    return {'DOUBLE': '=', 'TRIPLE': '#', 'AROMATIC': ''}[bond.bond_type]


def write_smiles(atoms: tuple[Atom, ...], bonds: tuple[Bond, ...]) -> str:
    """Serialize current graph records. Output is deterministic, NOT canonical.

    Preserve input aromatic annotations; do not infer chemical validity.
    A DFS forest provides branches; non-tree edges provide ring closures.
    """
    validate_graph(atoms, bonds)
    by_id = {a.id: a for a in atoms}
    position = {a.id: i for i, a in enumerate(atoms)}
    adjacency: dict[str, list[tuple[str, int]]] = {a.id: [] for a in atoms}
    for e, bond in enumerate(bonds):
        adjacency[bond.begin].append((bond.end, e))
        adjacency[bond.end].append((bond.begin, e))
    for values in adjacency.values():
        values.sort(key=lambda item: position[item[0]])
    visited: dict[str, int] = {}
    seen_edges: set[int] = set()
    children: dict[str, list[tuple[str, int]]] = {a.id: [] for a in atoms}
    closures: list[tuple[str, str, int]] = []
    roots: list[str] = []
    for atom in atoms:
        root = atom.id
        if root in visited:
            continue
        roots.append(root)
        visited[root] = len(visited)
        stack = [(root, iter(adjacency[root]))]
        while stack:
            u, iterator = stack[-1]
            item = next(iterator, None)
            if item is None:
                stack.pop()
                continue
            v, e = item
            if e in seen_edges:
                continue
            seen_edges.add(e)
            if v in visited:
                closures.append((u, v, e))
            else:
                children[u].append((v, e))
                visited[v] = len(visited)
                stack.append((v, iter(adjacency[v])))
    if len(closures) > 99:
        raise UnsupportedSmilesError('This writer supports at most 99 non-tree edges per molecule.')
    closures.sort(key=lambda item: (min(visited[item[0]], visited[item[1]]),
                                    max(visited[item[0]], visited[item[1]])))
    labels: dict[str, list[str]] = {a.id: [] for a in atoms}
    for number, (u, v, e) in enumerate(closures, 1):
        if visited[u] > visited[v]:
            u, v = v, u
        label = str(number) if number < 10 else f'%{number:02d}'
        labels[u].append(_bond_text(bonds[e], by_id) + label)
        labels[v].append(label)
    result: list[str] = []
    for component, root in enumerate(roots):
        if component:
            result.append('.')
        actions = [('node', root)]
        while actions:
            kind, value = actions.pop()
            if kind == 'text':
                result.append(value)
                continue
            result.append(_atom_text(by_id[value]) + ''.join(labels[value]))
            following = []
            for j, (child, e) in enumerate(children[value]):
                branch = j < len(children[value]) - 1
                if branch:
                    following.append(('text', '('))
                following.extend([('text', _bond_text(bonds[e], by_id)), ('node', child)])
                if branch:
                    following.append(('text', ')'))
            actions.extend(reversed(following))
    return ''.join(result)


@dataclass(frozen=True, init=False)
class Molecule:
    atoms: tuple[Atom, ...]
    bonds: tuple[Bond, ...]

    def __init__(self, smiles: str | None = None, *,
                 atoms: Iterable[Atom] | None = None,
                 bonds: Iterable[Bond] | None = None) -> None:
        if smiles is not None:
            if atoms is not None or bonds is not None:
                raise MoleculeError('Supply either smiles or explicit records, not both.')
            atom_records, bond_records = parse_smiles(smiles)
        else:
            if atoms is None:
                raise MoleculeError('Supply smiles or atoms.')
            atom_records, bond_records = tuple(atoms), tuple(bonds) if bonds is not None else ()
        validate_graph(atom_records, bond_records)
        object.__setattr__(self, 'atoms', atom_records)
        object.__setattr__(self, 'bonds', bond_records)

    @classmethod
    def from_smiles(cls, smiles: str) -> Molecule:
        return cls(smiles)

    def decode_to_smiles(self) -> str:
        return write_smiles(self.atoms, self.bonds)

    def validate(self) -> None:
        """Validate graph records only; this is not a chemistry validator."""
        validate_graph(self.atoms, self.bonds)

    def atom(self, atom_id: str) -> Atom:
        for atom in self.atoms:
            if atom.id == atom_id:
                return atom
        raise MoleculeError(f'Unknown atom ID: {atom_id!r}')

    def neighbors(self, atom_id: str) -> tuple[str, ...]:
        self.atom(atom_id)
        return tuple(b.end if b.begin == atom_id else b.begin
                     for b in self.bonds if atom_id in (b.begin, b.end))

    def bond_between(self, begin: str, end: str) -> Bond | None:
        self.atom(begin)
        self.atom(end)
        return next((b for b in self.bonds if {b.begin, b.end} == {begin, end}), None)

    def shortest_path(self, begin: str, end: str) -> tuple[str, ...] | None:
        self.atom(begin)
        self.atom(end)
        adjacency = {a.id: [] for a in self.atoms}
        for bond in self.bonds:
            adjacency[bond.begin].append(bond.end)
            adjacency[bond.end].append(bond.begin)
        previous: dict[str, str | None] = {begin: None}
        queue = deque([begin])
        while queue:
            current = queue.popleft()
            if current == end:
                path = []
                while current is not None:
                    path.append(current)
                    current = previous[current]
                return tuple(reversed(path))
            for other in adjacency[current]:
                if other not in previous:
                    previous[other] = current
                    queue.append(other)
        return None

    @property
    def total_charge(self) -> int:
        return sum(a.formal_charge for a in self.atoms)

    def replace_atom(self, atom_id: str, **changes) -> Molecule:
        old = self.atom(atom_id)
        if 'id' in changes:
            raise MoleculeError('Atom IDs cannot be changed by replace_atom.')
        new = replace(old, **changes)
        return Molecule(atoms=(new if a.id == atom_id else a for a in self.atoms), bonds=self.bonds)

    def replace_bond(self, begin: str, end: str, *, bond_type: str) -> Molecule:
        old = self.bond_between(begin, end)
        if old is None:
            raise MoleculeError('No bond exists between these atoms.')
        return Molecule(atoms=self.atoms,
                        bonds=(replace(b, bond_type=bond_type) if b is old else b for b in self.bonds))

    def to_code(self) -> str:
        """Return a constructor expression with current, authoritative records."""
        def record(value) -> str:
            arguments = []
            for field in fields(value):
                item = getattr(value, field.name)
                if field.name in ('id', 'element', 'begin', 'end') or item != field.default:
                    arguments.append(f'{field.name}={item!r}')
            return f'{type(value).__name__}(' + ', '.join(arguments) + ')'
        atom_lines = ''.join(f'        {record(a)},\n' for a in self.atoms)
        bond_lines = ''.join(f'        {record(b)},\n' for b in self.bonds)
        return f'Molecule(\n    atoms=(\n{atom_lines}    ),\n    bonds=(\n{bond_lines}    ),\n)'

    @classmethod
    def from_code(cls, source: str) -> Molecule:
        """Read constructors and literals only; never execute Python code.

        Resource limits are modest guards, not a general untrusted-code sandbox.
        """
        if not isinstance(source, str) or not 0 < len(source) <= MAX_TEXT:
            raise MoleculeError('Code is empty, not a string, or exceeds the input limit.')
        try:
            root = ast.parse(source, mode='eval')
            if sum(1 for _ in ast.walk(root)) > 40_000:
                raise MoleculeError('Too many syntax-tree nodes.')
            constructors = {'Atom': Atom, 'Bond': Bond, 'Molecule': cls}

            def read(node):
                if isinstance(node, (ast.Tuple, ast.List)):
                    return tuple(read(value) for value in node.elts)
                if isinstance(node, ast.Call):
                    if not isinstance(node.func, ast.Name) or node.func.id not in constructors:
                        raise MoleculeError('Only Molecule, Atom, and Bond constructors are allowed.')
                    if node.args or any(k.arg is None for k in node.keywords):
                        raise MoleculeError('Only explicit keyword arguments are allowed.')
                    names = [k.arg for k in node.keywords]
                    if len(names) != len(set(names)):
                        raise MoleculeError('Duplicate constructor keyword.')
                    return constructors[node.func.id](**{k.arg: read(k.value) for k in node.keywords})
                if isinstance(node, ast.Constant) and type(node.value) in (str, int, bool, type(None)):
                    return node.value
                if (isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd))
                        and isinstance(node.operand, ast.Constant) and type(node.operand.value) is int):
                    return -node.operand.value if isinstance(node.op, ast.USub) else node.operand.value
                raise MoleculeError('Only constructor calls and simple literal values are allowed.')

            molecule = read(root.body)
            if type(molecule) is not cls:
                raise MoleculeError('The expression must produce a Molecule.')
            return molecule
        except (SyntaxError, TypeError, RecursionError) as exc:
            raise MoleculeError(f'Malformed constructor expression: {exc}') from exc
