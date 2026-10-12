"""Stage-1 reproduction: what the professor's baseline prototype actually does.

Written from scratch for this round (2026-09-28). The historical log `raw.txt`
and report `verification_results.json` are NOT reused: their test script
(`test_standalone.py`) and data file (`first_5K.smi`) are absent from the
workspace, so the original experiments cannot be re-executed here.

Standard library only.
"""
from __future__ import annotations

import pathlib
import sys
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _harness import (  # noqa: E402
    EXPECTED_BASELINE_SHA256,
    PROTOTYPE_PATH,
    REPO_ROOT,
    baseline_sha256,
    load_prototype,
)

proto = load_prototype()
Atom = proto.Atom
Bond = proto.Bond
Molecule = proto.Molecule
MoleculeError = proto.MoleculeError
UnsupportedSmilesError = proto.UnsupportedSmilesError


class BaselineIntegrityTests(unittest.TestCase):
    def test_baseline_file_hash_is_unchanged(self):
        self.assertEqual(baseline_sha256(), EXPECTED_BASELINE_SHA256)

    def test_baseline_copy_matches_project_root_original(self):
        original = REPO_ROOT / 'molecule_standalone.py'
        self.assertTrue(original.is_file(), 'professor original missing from project root')
        self.assertEqual(PROTOTYPE_PATH.read_bytes(), original.read_bytes())

    def test_baseline_source_mentions_no_chemistry_package(self):
        text = PROTOTYPE_PATH.read_text(encoding='utf-8').lower()
        found = [n for n in ('rdkit', 'moljson', 'openbabel', 'pysmiles', 'networkx')
                 if n in text]
        self.assertEqual(found, [], f'baseline source mentions: {found}')

    def test_no_chemistry_module_in_sys_modules(self):
        bad = sorted(n for n in sys.modules
                     if n.split('.')[0].lower() in
                     {'rdkit', 'moljson', 'openbabel', 'pysmiles', 'networkx'})
        self.assertEqual(bad, [], f'external modules loaded: {bad}')

    def test_working_core_matches_the_baseline_byte_for_byte(self):
        """The deliverable core must be the professor original, not a variant.

        Everything added to Atom/Bond/Molecule in an earlier round was rolled
        back, so the working core and the baseline record have to stay
        identical. This test fails if anyone re-introduces a private extension
        into the core instead of putting it in the benchmark module.
        """
        core = PROTOTYPE_PATH.parent.parent / '01_implementation' / 'molecule_core.py'
        self.assertTrue(core.is_file(), 'working core missing')
        self.assertEqual(core.read_bytes(), PROTOTYPE_PATH.read_bytes())
        self.assertEqual(baseline_sha256(), EXPECTED_BASELINE_SHA256)

    def test_rolled_back_extensions_are_absent_from_the_core(self):
        text = (PROTOTYPE_PATH.parent.parent / '01_implementation' /
                'molecule_core.py').read_text(encoding='utf-8')
        for name in ('is_isomorphic', 'label_of', 'ATOM_LABEL_FIELDS',
                     '_closure_text', '_next_label'):
            self.assertNotIn(name, text, f'{name} was rolled back but is present')


class ParseAndWriteTests(unittest.TestCase):
    def test_cc_o_records(self):
        mol = Molecule('CCO')
        self.assertEqual(
            [(a.id, a.element, a.formal_charge, a.isotope, a.aromatic, a.explicit_hs, a.atom_map)
             for a in mol.atoms],
            [('a0', 'C', 0, 0, False, None, None),
             ('a1', 'C', 0, 0, False, None, None),
             ('a2', 'O', 0, 0, False, None, None)],
        )
        self.assertEqual([(b.begin, b.end, b.bond_type) for b in mol.bonds],
                         [('a0', 'a1', 'SINGLE'), ('a1', 'a2', 'SINGLE')])

    def test_cc_o_round_trips_exactly(self):
        self.assertEqual(Molecule('CCO').decode_to_smiles(), 'CCO')

    def test_branch(self):
        mol = Molecule('CC(C)O')
        self.assertEqual(mol.decode_to_smiles(), 'CC(C)O')
        self.assertEqual(mol.neighbors('a1'), ('a0', 'a2', 'a3'))

    def test_ring_closure(self):
        mol = Molecule('C1CC1')
        self.assertEqual(len(mol.bonds), 3)
        self.assertEqual(mol.bond_between('a0', 'a2').bond_type, 'SINGLE')
        self.assertEqual(mol.decode_to_smiles(), 'C1CC1')

    def test_two_digit_ring_label_is_accepted_then_renumbered_by_the_writer(self):
        # The parser accepts %NN; the writer renumbers closures from 1, so the
        # string is not required to come back byte-identical.
        mol = Molecule('C%10CC%10')
        self.assertEqual(len(mol.atoms), 3)
        self.assertEqual(mol.decode_to_smiles(), 'C1CC1')
        self.assertEqual(Molecule(mol.decode_to_smiles()).decode_to_smiles(), 'C1CC1')

    def test_disconnected_components(self):
        mol = Molecule('CCO.CC')
        self.assertEqual(len(mol.atoms), 5)
        self.assertEqual(mol.decode_to_smiles(), 'CCO.CC')

    def test_bond_types(self):
        self.assertEqual([b.bond_type for b in Molecule('C=C').bonds], ['DOUBLE'])
        self.assertEqual([b.bond_type for b in Molecule('C#N').bonds], ['TRIPLE'])
        aromatic = Molecule('c1ccccc1')
        self.assertEqual([b.bond_type for b in aromatic.bonds], ['AROMATIC'] * 6)
        self.assertTrue(all(a.aromatic for a in aromatic.atoms))

    def test_long_chain_without_recursion(self):
        # MAX_ATOMS is 10000, so 2000 works; 15000 exceeds the declared limit.
        mol = Molecule('C' * 2000)
        self.assertEqual(len(mol.atoms), 2000)
        self.assertEqual(mol.decode_to_smiles(), 'C' * 2000)
        with self.assertRaises(UnsupportedSmilesError) as ctx:
            Molecule('C' * 15000)
        self.assertIn('exceeds the configured limit', str(ctx.exception))

    def test_output_reparses_with_same_records(self):
        cases = ('CCO', 'CC(C)O', 'C1CC1', 'CCO.CC', 'C=C', 'C#N', 'c1ccccc1', 'C%10CC%10',
                 'CC(=O)[O-]', '[13CH4]', '[CH3:1][OH:2]', '[NH4+]')
        for text in cases:
            mol = Molecule(text)
            again = Molecule(mol.decode_to_smiles())
            self.assertEqual(len(again.atoms), len(mol.atoms), text)
            self.assertEqual(len(again.bonds), len(mol.bonds), text)
            for bond in again.bonds:
                original = mol.bond_between(bond.begin, bond.end)
                self.assertIsNotNone(original, text)
                self.assertEqual(bond.bond_type, original.bond_type, text)

    def test_writer_is_deterministic_for_a_fixed_ordered_graph(self):
        mol = Molecule('CC(C)C(=O)OC')
        self.assertEqual(mol.decode_to_smiles(), mol.decode_to_smiles())


class HydrogenPolicyTests(unittest.TestCase):
    def test_four_documented_policies(self):
        self.assertIsNone(Molecule('C').atoms[0].explicit_hs)
        self.assertEqual(Molecule('[C]').atoms[0].explicit_hs, 0)
        self.assertEqual(Molecule('[CH3]').atoms[0].explicit_hs, 3)
        self.assertEqual(Molecule('[nH]').atoms[0].explicit_hs, 1)

    def test_bracket_required_atom_must_declare_hydrogen_count(self):
        with self.assertRaises(MoleculeError) as ctx:
            Atom(id='x', element='C', formal_charge=1)
        self.assertIn('requires brackets', str(ctx.exception))

    def test_ammonium_from_markdown(self):
        ammonium = Molecule(atoms=(Atom(id='n', element='N', formal_charge=1, explicit_hs=4),))
        self.assertEqual(ammonium.decode_to_smiles(), '[NH4+]')
        self.assertEqual(ammonium.total_charge, 1)

    def test_bracket_hydrogen_cannot_carry_its_own_h_token(self):
        with self.assertRaises(MoleculeError):
            Molecule('[HH]')

    def test_separate_hydrogen_vertex(self):
        mol = Molecule('[H]C')
        self.assertEqual([a.element for a in mol.atoms], ['H', 'C'])
        self.assertEqual(mol.decode_to_smiles(), '[H]C')

    def test_aromatic_nh_inside_a_ring(self):
        mol = Molecule('c1cc[nH]c1')
        self.assertEqual(mol.atoms[3].explicit_hs, 1)
        self.assertTrue(mol.atoms[3].aromatic)


class ChargeIsotopeMapTests(unittest.TestCase):
    def test_charge_round_trips(self):
        for text, charge in (('[NH4+]', 1), ('[O-]', -1), ('[Fe+2]', 2), ('[Fe-2]', -2)):
            mol = Molecule(text)
            self.assertEqual(mol.atoms[0].formal_charge, charge, text)
            self.assertEqual(mol.decode_to_smiles(), text, text)

    def test_isotope_round_trips(self):
        mol = Molecule('[13CH4]')
        self.assertEqual(mol.atoms[0].isotope, 13)
        self.assertEqual(mol.decode_to_smiles(), '[13CH4]')

    def test_atom_map_round_trips(self):
        mol = Molecule('[CH3:1][OH:2]')
        self.assertEqual([a.atom_map for a in mol.atoms], [1, 2])
        self.assertEqual(mol.decode_to_smiles(), '[CH3:1][OH:2]')

    def test_total_charge(self):
        self.assertEqual(Molecule('CC(=O)[O-]').total_charge, -1)
        self.assertEqual(Molecule('CCO').total_charge, 0)


class QueryTests(unittest.TestCase):
    def test_neighbors(self):
        mol = Molecule('CCO')
        self.assertEqual(mol.neighbors('a0'), ('a1',))
        self.assertEqual(mol.neighbors('a1'), ('a0', 'a2'))
        self.assertEqual(mol.neighbors('a2'), ('a1',))

    def test_bond_between(self):
        mol = Molecule('CCO')
        self.assertIsNone(mol.bond_between('a0', 'a2'))
        self.assertEqual(mol.bond_between('a0', 'a1').bond_type, 'SINGLE')
        self.assertEqual(mol.bond_between('a1', 'a0').bond_type, 'SINGLE')

    def test_shortest_path(self):
        self.assertEqual(Molecule('CCO').shortest_path('a0', 'a2'), ('a0', 'a1', 'a2'))
        self.assertIsNone(Molecule('CCO.CC').shortest_path('a0', 'a4'))
        # In a triangle both routes are length 2; BFS returns whichever the
        # bond order yields first, so the exact tuple is not a spec guarantee.
        self.assertEqual(Molecule('C1CC1').shortest_path('a0', 'a2'), ('a0', 'a2'))


class EditSemanticsTests(unittest.TestCase):
    def test_oxygen_to_sulfur_yields_ccs_and_leaves_original(self):
        mol = Molecule('CCO')
        edited = mol.replace_atom('a2', element='S')
        self.assertEqual(edited.decode_to_smiles(), 'CCS')
        self.assertEqual(mol.decode_to_smiles(), 'CCO')

    def test_records_after_edit(self):
        edited = Molecule('CCO').replace_atom('a2', element='S')
        self.assertEqual([(a.id, a.element) for a in edited.atoms],
                         [('a0', 'C'), ('a1', 'C'), ('a2', 'S')])
        self.assertEqual([(b.begin, b.end, b.bond_type) for b in edited.bonds],
                         [('a0', 'a1', 'SINGLE'), ('a1', 'a2', 'SINGLE')])

    def test_atom_id_cannot_be_changed(self):
        with self.assertRaises(MoleculeError) as ctx:
            Molecule('C1CC1').replace_atom('a0', id='z9')
        self.assertIn('cannot be changed', str(ctx.exception))

    def test_invalid_replacement_is_rejected(self):
        mol = Molecule('CCO')
        with self.assertRaises(MoleculeError):
            mol.replace_atom('a2', element='Xx')
        with self.assertRaises(MoleculeError):
            mol.replace_atom('a2', element='Na', explicit_hs=None)
        with self.assertRaises(MoleculeError):
            mol.replace_atom('a2', element='Se', aromatic=True, explicit_hs=None)

    def test_edit_can_produce_a_bracket_form_that_changes_hydrogen_semantics(self):
        # S is in the bare organic subset, so explicit_hs=0 forces brackets and
        # the decoded string becomes CC[S] (zero explicit H) instead of CCS.
        edited = Molecule('CCO').replace_atom('a2', element='S', explicit_hs=0)
        self.assertEqual(edited.decode_to_smiles(), 'CC[S]')

    def test_edit_can_produce_a_chemically_invalid_graph(self):
        mol = Molecule('CCO')
        five_bonded_carbon = mol.replace_atom('a2', element='C', explicit_hs=0)
        self.assertEqual(five_bonded_carbon.decode_to_smiles(), 'CC[C]')

    def test_charge_edit(self):
        edited = Molecule('CCO').replace_atom('a2', element='O', formal_charge=-1, explicit_hs=0)
        self.assertEqual(edited.decode_to_smiles(), 'CC[O-]')

    def test_bond_edit(self):
        mol = Molecule('CCO')
        edited = mol.replace_bond('a0', 'a1', bond_type='DOUBLE')
        self.assertEqual(edited.decode_to_smiles(), 'C=CO')
        self.assertEqual(mol.decode_to_smiles(), 'CCO')

    def test_replace_bond_requires_an_existing_bond(self):
        with self.assertRaises(MoleculeError):
            Molecule('CCO').replace_bond('a0', 'a2', bond_type='DOUBLE')

    def test_frozen_records_and_molecule_are_hashable_immutable_values(self):
        mol = Molecule('CCO')
        with self.assertRaises(Exception):
            mol.atoms[0].element = 'N'
        with self.assertRaises(Exception):
            mol.atoms = ()


class ConstructorCodeTests(unittest.TestCase):
    def test_to_code_shape(self):
        self.assertEqual(
            Molecule('CCO').to_code(),
            "Molecule(\n"
            "    atoms=(\n"
            "        Atom(id='a0', element='C'),\n"
            "        Atom(id='a1', element='C'),\n"
            "        Atom(id='a2', element='O'),\n"
            "    ),\n"
            "    bonds=(\n"
            "        Bond(begin='a0', end='a1'),\n"
            "        Bond(begin='a1', end='a2'),\n"
            "    ),\n"
            ")",
        )

    def test_round_trip(self):
        mol = Molecule('CCO')
        rebuilt = Molecule.from_code(mol.to_code())
        self.assertEqual(rebuilt.decode_to_smiles(), 'CCO')
        for a, b in zip(mol.atoms, rebuilt.atoms):
            self.assertEqual((a.id, a.element, a.formal_charge, a.isotope, a.aromatic,
                              a.explicit_hs, a.atom_map),
                             (b.id, b.element, b.formal_charge, b.isotope, b.aromatic,
                              b.explicit_hs, b.atom_map))
        for a, b in zip(mol.bonds, rebuilt.bonds):
            self.assertEqual((a.begin, a.end, a.bond_type), (b.begin, b.end, b.bond_type))

    def test_round_trip_for_annotated_atoms(self):
        for text in ('CC(=O)[O-]', '[13CH4]', '[CH3:1][OH:2]', 'c1cc[nH]c1', '[NH4+]'):
            mol = Molecule(text)
            self.assertEqual(Molecule.from_code(mol.to_code()).decode_to_smiles(),
                             mol.decode_to_smiles(), text)

    def test_from_code_rejects_non_constructor_expressions(self):
        bad = (
            'dict()',
            'Molecule(atoms=1+1)',
            'Molecule(atoms=__import__("os").system("echo hi"))',
            "Molecule(atoms=(Atom(id='a0', element='C').id,))",
            "Molecule(atoms=(Atom(id='a0', element='C').__class__,))",
            "Molecule(atoms=(Atom(id='a0', element='C'),), bonds=(), extra=1)",
            "Molecule((Atom(id='a0', element='C'),), ())",
        )
        for source in bad:
            with self.assertRaises(MoleculeError, msg=source):
                Molecule.from_code(source)

    def test_from_code_rejects_duplicate_keywords(self):
        with self.assertRaises(MoleculeError):
            Molecule.from_code("Atom(id='a0', element='C', id='a1')")

    def test_from_code_enforces_text_length_limit(self):
        source = "Molecule(atoms=(Atom(id='a0', element='C'),), bonds=())"
        with self.assertRaises(MoleculeError) as ctx:
            Molecule.from_code(source * 2000)
        self.assertTrue(str(ctx.exception).startswith('Code is empty, not a string, or exceeds'))

    def test_from_code_enforces_ast_node_limit(self):
        # Deep nesting trips CPython's parser depth limit before the node-count
        # guard, so the surfaced error is "too many nested parentheses".
        with self.assertRaises(MoleculeError) as ctx:
            Molecule.from_code('Molecule(atoms=' + '(' * 20000 + ')' * 20000 + ')')
        self.assertIn('too many nested parentheses', str(ctx.exception))

    def test_from_code_does_not_execute_arbitrary_python(self):
        marker = pathlib.Path(__file__).resolve().parent / '__never_created__.marker'
        self.assertFalse(marker.exists(), 'stale marker from an earlier run')
        with self.assertRaises(MoleculeError):
            Molecule.from_code("Molecule(atoms=open('%s', 'w').write('x'),)" % marker)
        self.assertFalse(marker.exists(), 'from_code() executed arbitrary code')


class RejectionTests(unittest.TestCase):
    """Declared-unsupported input must be refused, never silently dropped."""

    def test_stereochemistry_is_rejected(self):
        for text in ('F/C=C/F', 'F/C=C\\F', 'C[C@H](N)O', 'C[C@@H](N)O'):
            with self.assertRaises(UnsupportedSmilesError, msg=text):
                Molecule(text)

    def test_cxsmiles_reaction_and_query_syntax_is_rejected(self):
        for text in ('CC |$;;_AP1$|', 'CC>>CC', 'C~C', 'C?C', 'C&C'):
            with self.assertRaises(UnsupportedSmilesError, msg=text):
                Molecule(text)

    def test_malformed_smiles_is_rejected(self):
        for text in ('CC(', 'C)', 'C1CC', 'CC=', '=CC', 'C..C', 'CC.', '.CC', 'C%1CC', 'C%123'):
            with self.assertRaises(MoleculeError, msg=text):
                Molecule(text)

    def test_names_and_whitespace_are_rejected(self):
        for text in ('ethanol', 'CC O', 'CC\t'):
            with self.assertRaises(UnsupportedSmilesError, msg=text):
                Molecule(text)

    def test_bad_graph_records_are_rejected(self):
        with self.assertRaises(MoleculeError):
            Molecule(atoms=(), bonds=())
        with self.assertRaises(MoleculeError):
            Molecule(atoms=(Atom(id='a0', element='C'), Atom(id='a0', element='O')))
        with self.assertRaises(MoleculeError):
            Molecule(atoms=(Atom(id='a0', element='C'),),
                     bonds=(Bond(begin='a0', end='a9'),))
        with self.assertRaises(MoleculeError):
            Molecule(atoms=(Atom(id='a0', element='C'), Atom(id='a1', element='C')),
                     bonds=(Bond(begin='a0', end='a1'), Bond(begin='a1', end='a0')))
        with self.assertRaises(MoleculeError) as ctx:
            Molecule(atoms=(Atom(id='a0', element='C'), Atom(id='a1', element='C')),
                     bonds=(Bond(begin='a0', end='a0'),))
        self.assertIn('Self-bonds and duplicate bonds', str(ctx.exception))

    def test_aromatic_bond_between_non_aromatic_atoms_is_rejected(self):
        with self.assertRaises(MoleculeError) as ctx:
            Molecule(atoms=(Atom(id='a0', element='C'), Atom(id='a1', element='C')),
                     bonds=(Bond(begin='a0', end='a1', bond_type='AROMATIC'),))
        self.assertIn('aromatic bond must join two aromatic atoms', str(ctx.exception))

    def test_aromatic_atom_of_unsupported_element_is_rejected(self):
        with self.assertRaises(MoleculeError) as ctx:
            Atom(id='x', element='Na', aromatic=True)
        self.assertIn('Aromatic notation is unsupported', str(ctx.exception))

    def test_non_iterable_bonds_argument_raises_type_error_not_molecule_error(self):
        with self.assertRaises(TypeError) as ctx:
            Molecule(atoms=(Atom(id='a0', element='C'),), bonds=Bond(begin='a0', end='a0'))
        self.assertIn('not iterable', str(ctx.exception))


if __name__ == '__main__':
    unittest.main(verbosity=2)
