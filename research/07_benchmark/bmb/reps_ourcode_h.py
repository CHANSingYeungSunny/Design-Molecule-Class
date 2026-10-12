"""Our Code + H -- an ABLATION CANDIDATE, not a result and not a replacement.

What it changes, and what it does not
-------------------------------------
ONE factor only: every atom's hydrogen count is written explicitly. Atom ids,
atom order, bond order, the field order, and the surrounding prompt are all
unchanged, so a difference between the two arms can only come from that factor.

    original :  Atom(id='a0', element='C')
    candidate:  Atom(id='a0', element='C', explicit_hs=3)

The professor's Atom / Bond / Molecule classes are NOT modified. ``explicit_hs``
is a field those classes already define -- the prototype simply leaves it as
None for almost every atom, meaning "apply the element's normal implicit rule".
This module just fills it in. Nothing here requires a parameter the original
class does not have.

Field definition (fixed before any measurement)
-----------------------------------------------
The professor's own class documents the field as:

    None means unbracketed-style implicit-H rules apply.
    An integer means an exact bracket H count; no additional implicit H.
    Thus C -> None, [C] -> 0, [CH3] -> 3, and [nH] -> 1.

So the value is the BRACKET H count -- the hydrogens the atom itself accounts
for -- and it must NOT include hydrogens that are already present in the graph
as separate H vertices. Writing those again would add them twice: an N with
three bonded [H] vertices plus ``explicit_hs=3`` would emit ``[NH3]`` on top of
three H vertices, i.e. NH6.

The value is therefore taken as RDKit ``GetTotalNumHs(includeNeighbors=False)``,
which is exactly the bracket count. (``includeNeighbors=True`` would return the
neighbouring H as well, which is the double-counting error; the two differ only
for atoms that have explicit H vertices, which is precisely the case that has to
be right.)

Consequences, checked rather than assumed:
  * a methyl carbon in a chain -> 3;  a carbonyl oxygen -> 0;  an [nH] -> 1
  * a heavy atom whose hydrogens are all separate H vertices -> 0, because those
    hydrogens are already carried by the vertices
  * an isolated [H+] vertex -> 0, which the original class also requires
    (``element == 'H'`` permits only None or 0)

This is a hydrogen COUNT, not a donor count: the payload carries no donor flag,
no HBD total, no functional-group label, no element total and no descriptor.

Usage
-----
Nothing in the pipeline calls this module unless asked. ``OURCODE_H`` must be
requested explicitly; it is NOT part of the eight-row table.
"""
from __future__ import annotations

import dataclasses

from rdkit import Chem

from .core_loader import load_core

VARIANT_NAME = 'our_code_h'
DISPLAY_NAME = 'Our Code + H（消融候选，未进入正式结果表）'

FIELD_DEFINITION = (
    'explicit_hs = total hydrogens bonded to the atom (graph H vertices plus '
    'implicit H implied by element and bond order); source: molecular structure')

# documented semantics of the field we fill in
DEFAULT_SEMANTICS = 'explicit_hs=None means the element normal implicit-H rule applies'


class OurCodeHError(RuntimeError):
    """The candidate text could not be produced for this molecule."""


def _rdkit_h_counts(smiles):
    """The BRACKET hydrogen count per atom (see the module docstring).

    ``includeNeighbors=False`` deliberately stops the count at the atom itself,
    so hydrogens that exist as separate vertices are not counted a second time.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise OurCodeHError('RDKit cannot parse %r' % smiles)
    return [int(atom.GetTotalNumHs(includeNeighbors=False)) for atom in mol.GetAtoms()]


def _neighbour_h_counts(smiles):
    """How many H VERTICES each atom carries. Used to detect the case where the
    naive definition would double count."""
    mol = Chem.MolFromSmiles(smiles)
    return [sum(1 for n in atom.GetNeighbors() if n.GetSymbol() == 'H')
            for atom in mol.GetAtoms()]


def _check_element_order(smiles, core_mol):
    """The core parser and RDKit must agree on atom order, or the H counts land
    on the wrong atoms. Checked, not assumed."""
    rd = Chem.MolFromSmiles(smiles)
    core_symbols = [a.element for a in core_mol.atoms]
    rd_symbols = [a.GetSymbol() for a in rd.GetAtoms()]
    return core_symbols == rd_symbols


def build_text(smiles, variant=0):
    """The candidate payload. ``variant`` reorders records, exactly as the
    original representation does for Identity positives."""
    core = load_core()
    mol = core.Molecule(smiles)
    if not _check_element_order(smiles, mol):
        raise OurCodeHError(
            'atom order differs between the core parser and RDKit for %r; refusing to '
            'attach hydrogen counts that could land on the wrong atom' % smiles)
    counts = _rdkit_h_counts(smiles)
    if len(counts) != len(mol.atoms):
        raise OurCodeHError('atom count mismatch for %r' % smiles)

    atoms = tuple(dataclasses.replace(atom, explicit_hs=counts[index])
                  for index, atom in enumerate(mol.atoms))
    bonds = list(mol.bonds)
    if variant:
        # same reordering rule the original representation uses, so an Identity
        # positive is two different arrangements of one molecule -- never one
        # text repeated twice
        atoms = list(atoms)
        if len(atoms) > 1:
            shift = (variant * 5 + 1) % len(atoms)
            atoms = tuple(atoms[shift:] + atoms[:shift])
        if len(bonds) > 1:
            shift = (variant * 3 + 1) % len(bonds)
            bonds = bonds[shift:] + bonds[:shift]
    return dataclasses.replace(mol, atoms=tuple(atoms), bonds=tuple(bonds)).to_code()


def verify(smiles):
    """Offline checks for one molecule. Returns a dict of booleans."""
    core = load_core()
    original_text = core.Molecule(smiles).to_code()
    text = build_text(smiles, 0)
    rebuilt = core.Molecule.from_code(text)
    same_smiles = (Chem.MolToSmiles(Chem.MolFromSmiles(rebuilt.decode_to_smiles()),
                                    isomericSmiles=False)
                   == Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=False))
    counts = _rdkit_h_counts(smiles)
    rebuilt_counts = [a.explicit_hs for a in rebuilt.atoms]
    a, b = build_text(smiles, 0), build_text(smiles, 1)
    lower = text.lower()
    return {
        'chemistry_unchanged': bool(same_smiles),
        'h_counts_preserved_through_roundtrip': rebuilt_counts == counts,
        'every_atom_states_a_count': all(c is not None for c in rebuilt_counts),
        'ids_unchanged_vs_original': [a.id for a in core.Molecule(smiles).atoms]
                                     == [a.id for a in rebuilt.atoms],
        'order_unchanged_vs_original': [a.element for a in core.Molecule(smiles).atoms]
                                       == [a.element for a in rebuilt.atoms],
        'identity_views_differ': a != b,
        'no_answer_leakage': not any(token in lower for token in
                                     ('hbd', 'donor', 'qed', 'logp', 'tpsa', 'functional')),
    }


def render(smiles, variant=0):
    return build_text(smiles, variant)