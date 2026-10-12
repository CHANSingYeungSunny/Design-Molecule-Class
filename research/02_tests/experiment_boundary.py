"""Stage-2 boundary experiments against the professor's baseline prototype.

Every number and message printed here was produced by running this file with
`python -S -B experiment_boundary.py`. Nothing is copied from the historical
`raw.txt` / `verification_results.json`.

Purpose: measure where the parse -> graph -> write path is consistent and where
it is not, and separate reproduced defects from declared limitations.
"""
from __future__ import annotations

import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _harness import load_prototype  # noqa: E402

proto = load_prototype()
Atom = proto.Atom
Bond = proto.Bond
Molecule = proto.Molecule

RESULTS = []


def probe(name, fn):
    try:
        value = fn()
    except Exception as exc:  # noqa: BLE001 - the failure text is the evidence
        RESULTS.append({'probe': name, 'outcome': 'raised', 'error': type(exc).__name__,
                        'message': str(exc)[:300]})
        print('[raised   ] {}: {}: {}'.format(name, type(exc).__name__, str(exc)[:200]))
    else:
        text = value if isinstance(value, str) else repr(value)
        RESULTS.append({'probe': name, 'outcome': 'returned', 'value': str(text)[:300]})
        print('[returned ] {}: {}'.format(name, str(text)[:200]))


def banner(title):
    print()
    print('=' * 78)
    print(title)
    print('=' * 78)


banner('E1  ring-closure capacity: does parse accept what the writer refuses?')
text = '.'.join(['C1CC1'] * 100)
print('input text length =', len(text))
probe('E1a parse 100 cyclopropanes', lambda: len(Molecule(text).atoms))
probe('E1b write 100 cyclopropanes', lambda: Molecule(text).decode_to_smiles())
probe('E1c parse 99 cyclopropanes', lambda: len(Molecule('.'.join(['C1CC1'] * 99)).atoms))
probe('E1d write 99 cyclopropanes', lambda: Molecule('.'.join(['C1CC1'] * 99)).decode_to_smiles()[:40])
for n in (50, 98, 99, 100):
    probe('E1e write {} cyclopropanes'.format(n),
          lambda n=n: Molecule('.'.join(['C1CC1'] * n)).decode_to_smiles()[:20])

banner('E2  ring-label numbering: reuse, sharing, renumbering')
probe('E2a two rings in one component reuse label 1', lambda: Molecule('C1CC1C1CC1').decode_to_smiles())
probe('E2b fused bicyclic', lambda: Molecule('C1CCC2CCCC2C1').decode_to_smiles())
probe('E2c writer renumbers %10 input to the lowest free label', lambda: Molecule('C%10CC%10').decode_to_smiles())
probe('E2d benzene', lambda: Molecule('c1ccccc1').decode_to_smiles())
probe('E2e decalin', lambda: Molecule('C1CCC2CCCCC2C1').decode_to_smiles())

banner('E3  constructor-code round-trip capacity (to_code -> from_code)')
for size in (100, 400, 500, 600, 1000, 2000):
    mol = Molecule('C' * size)
    code = mol.to_code()
    try:
        back = Molecule.from_code(code)
        ok = 'ok atoms={} smiles_match={}'.format(
            len(back.atoms), back.decode_to_smiles() == mol.decode_to_smiles())
    except Exception as exc:  # noqa: BLE001
        ok = '{}: {}'.format(type(exc).__name__, str(exc)[:120])
    print('  atoms={:5d}  to_code() length={:7d}  from_code() -> {}'.format(size, len(code), ok))
    RESULTS.append({'probe': 'E3 atoms={}'.format(size), 'outcome': 'info', 'value': str(ok)})
print('  MAX_TEXT (from_code text limit) =', proto.MAX_TEXT)
print('  MAX_ATOMS (parse atom limit)    =', proto.MAX_ATOMS)

banner('E4  to_code() gives no hint that its own output may be unreadable')
code2000 = Molecule('C' * 2000).to_code()
print('  to_code() on 2000 atoms succeeded, length =', len(code2000), ' MAX_TEXT =', proto.MAX_TEXT)
probe('E4a from_code() on that code', lambda: Molecule.from_code(code2000))
probe('E4b CCO still round-trips', lambda: Molecule.from_code(Molecule('CCO').to_code()).decode_to_smiles())

banner('E5  equality: record order vs graph structure')
mol = Molecule('CCO')
shuffled = Molecule(atoms=tuple(reversed(mol.atoms)), bonds=tuple(reversed(mol.bonds)))
probe('E5a original == record-shuffled copy', lambda: mol == shuffled)
probe('E5b original == re-parsed CCO', lambda: mol == Molecule('CCO'))
renamed = Molecule(atoms=(Atom('n0', 'C'), Atom('n1', 'C'), Atom('n2', 'O')),
                   bonds=(Bond('n0', 'n1'), Bond('n1', 'n2')))
probe('E5c CCO == same graph with other atom ids', lambda: mol == renamed)
probe('E5d CCO == CCS', lambda: Molecule('CCO') == Molecule('CCS'))
probe('E5e neighbors() order follows bond record order',
      lambda: (Molecule('CC(C)O').neighbors('a1'), Molecule('CCO').neighbors('a1')))
probe('E5f bond_between() is order-insensitive',
      lambda: (Molecule('CCO').bond_between('a0', 'a1').begin,
               Molecule('CCO').bond_between('a1', 'a0').begin))
print('  dataclass __eq__ compares the (atoms, bonds) tuples position by position.')

banner('E6  hydrogen policy after editing')
probe('E6a CCO -> S (implicit H path)', lambda: Molecule('CCO').replace_atom('a2', element='S').decode_to_smiles())
probe('E6b CCO -> S with explicit_hs=0', lambda: Molecule('CCO').replace_atom('a2', element='S', explicit_hs=0).decode_to_smiles())
probe('E6c CCO -> Na, explicit_hs omitted', lambda: Molecule('CCO').replace_atom('a2', element='Na').decode_to_smiles())
probe('E6d CCO -> Na, explicit_hs=0', lambda: Molecule('CCO').replace_atom('a2', element='Na', explicit_hs=0).decode_to_smiles())
probe('E6e CCO -> aromatic c (implicit H path)', lambda: Molecule('CCO').replace_atom('a2', element='C', aromatic=True).decode_to_smiles())
probe('E6f CCO -> O with charge -1', lambda: Molecule('CCO').replace_atom('a2', element='O', formal_charge=-1, explicit_hs=0).decode_to_smiles())

banner('E7  graph consistency checks vs chemical validity')
five_bonded = Molecule(
    atoms=(Atom('a0', 'C'), Atom('a1', 'C'), Atom('a2', 'C'), Atom('a3', 'C'), Atom('a4', 'C'), Atom('a5', 'C')),
    bonds=(Bond('a0', 'a1'), Bond('a0', 'a2'), Bond('a0', 'a3'), Bond('a0', 'a4'), Bond('a0', 'a5')))
probe('E7a five single bonds on one carbon are accepted', lambda: five_bonded.decode_to_smiles())
probe('E7b that record count', lambda: (len(five_bonded.atoms), len(five_bonded.bonds)))
probe('E7c aromatic flag alone is accepted', lambda: Molecule(atoms=(Atom('a0', 'C', aromatic=True),), bonds=()).decode_to_smiles())
probe('E7d aromatic bond between two non-aromatic carbons is refused', lambda: Molecule(
    atoms=(Atom('a0', 'C'), Atom('a1', 'C')),
    bonds=(Bond('a0', 'a1', 'AROMATIC'),)).decode_to_smiles())
probe('E7e aromatic bond between two aromatic carbons is accepted', lambda: Molecule(
    atoms=(Atom('a0', 'C', aromatic=True), Atom('a1', 'C', aromatic=True)),
    bonds=(Bond('a0', 'a1', 'AROMATIC'),)).decode_to_smiles())
print('  graph consistency = unique ids, resolvable endpoints, no self/duplicate bond,')
print('  aromatic bond joins two aromatic atoms.')
print('  chemical validity = valence, kekulization, radicals, aromaticity perception. NOT checked.')

banner('E8  traversal order and internal IDs')
probe('E8a write 1000-carbon chain (first 20 chars)', lambda: Molecule('C' * 1000).decode_to_smiles()[:20])
probe('E8b branched input order preserved', lambda: Molecule('CC(C)CC(C)C').decode_to_smiles())
probe('E8c to_code() text depends on the internal ids',
      lambda: Molecule('CCO').to_code() == Molecule(
          atoms=(Atom('x0', 'C'), Atom('x1', 'C'), Atom('x2', 'O')),
          bonds=(Bond('x0', 'x1'), Bond('x1', 'x2'))).to_code())
probe('E8d same graph, professor-style ids, to_code()', lambda: Molecule(
    atoms=(Atom('a0', 'C'), Atom('a1', 'C'), Atom('a2', 'O')),
    bonds=(Bond('a0', 'a1'), Bond('a1', 'a2'))).to_code())

banner('E9  validate_graph() boundaries')
probe('E9a empty atom list', lambda: Molecule(atoms=(), bonds=()))
probe('E9b exactly MAX_ATOMS (10000)', lambda: len(Molecule('C' * 10000).atoms))
probe('E9c MAX_ATOMS + 1', lambda: len(Molecule('C' * 10001).atoms))
probe('E9d Molecule() with no arguments at all', lambda: Molecule())
probe('E9e non-iterable bonds argument', lambda: Molecule(atoms=(Atom('a0', 'C'),), bonds=Bond('a0', 'a0')))

banner('E10  declared limitations re-verified this round')
probe('E10a tetrahedral @ is rejected', lambda: Molecule('C[C@H](N)O'))
probe('E10b directional / is rejected', lambda: Molecule('F/C=C/F'))
probe('E10c CXSMILES layer is rejected', lambda: Molecule('CC |$;;_AP1$|'))
probe('E10d names are rejected', lambda: Molecule('ethanol'))
probe('E10e five-bond carbon still serializes', lambda: five_bonded.decode_to_smiles())

banner('E11  atom counts vs the declared MAX_ATOMS gate in parse_smiles()')
probe('E11a 10000 atoms parse', lambda: len(Molecule('C' * 10000).atoms))
probe('E11b 10001 atoms parse', lambda: len(Molecule('C' * 10001).atoms))

out = pathlib.Path(__file__).resolve().parent.parent / '05_results' / 'boundary_probe_output.json'
out.write_text(json.dumps({'python': sys.version.split()[0], 'probes': RESULTS}, indent=2, ensure_ascii=False),
               encoding='utf-8')
print()
print('wrote', out)
