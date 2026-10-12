"""Reference labels (oracle) computed by RDKit.

Split of roles, per the protocol:
  task      decides what the model is asked
  oracle    produces the reference value here (RDKit), never the model
  evaluator decides the score from a frozen rule (scoring.py)

The core class (Atom/Bond/Molecule) is deliberately absent from this file. It
may not mark its own homework: identity, structure and property labels all come
from a fixed external toolchain.

RDKit version is pinned in dataset_manifest.json. QED uses QED.qed defaults --
the weights are stored with the manifest so the reference is reproducible.
"""
from __future__ import annotations

import hashlib
import json

from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, QED, rdMolDescriptors

from . import paths

RDLogger.DisableLog('rdApp.*')

# Query elements for T3, chosen from the measured heavy-atom distribution of the
# frozen pool (dataset_manifest.json -> task_support.atom_counts), not copied from
# the paper's example set. Measured on the 2,101-molecule pool:
#   C  2005/2101 non-zero (95%), range 1..83
#   O  1851/2101 non-zero (88%), range 1..28
#   N  1166/2101 non-zero (56%), range 1..14
#   S   353/2101 non-zero (17%), range 1..9
#   Cl  284/2101 non-zero (13%)      F  103/2101 non-zero (4.9%)
# The paper's example set used F, but F is non-zero in only 4.9% of this pool,
# so a model that always answers "0 fluorine" would score well on it without
# understanding the structure. S carries real spread and is used instead.
COUNT_ELEMENTS = ('C', 'N', 'O', 'S')

# T5 descriptors. Round 1 runs HBD only; the rest are implemented and reported as
# not-yet-run rather than silently omitted.
DESCRIPTORS_FULL = ('hbd', 'hba', 'logp', 'tpsa')
DESCRIPTORS_ROUND1 = ('hbd',)


def rdkit_version() -> str:
    import rdkit
    return rdkit.__version__


def qed_weights() -> dict:
    """The exact parameters behind QED.qed's default, captured from the module.

    QED.qed(mol) with no `w` argument uses WEIGHT_MEAN over eight descriptors,
    each mapped through an ADS (adsorption/desorption) curve. Both tables are
    read out of the installed RDKit rather than transcribed, so the reference is
    reproducible from this record alone.
    """
    weight = QED.WEIGHT_MEAN
    return {
        '_note': 'RDKit QED.qed defaults; no custom weights are passed by this project.',
        'weight_name': 'WEIGHT_MEAN',
        'weights': {name: getattr(weight, name) for name in weight._fields},
        'descriptor_order': list(weight._fields),
        'ads_parameters': {
            name: {field: getattr(param, field) for field in param._fields}
            for name, param in QED.adsParameters.items()
        },
        'descriptor_version': getattr(QED.qed, 'version', None),
        'source': 'rdkit.Chem.QED (Novartis BSD-3-Clause), read from the installed package',
    }


def load_candidate_functional_groups() -> dict:
    """The raw candidate list, including labels the support rule may disable."""
    with open(paths.FUNCTIONAL_GROUPS, encoding='utf-8') as handle:
        return json.load(handle)


def load_functional_groups() -> dict:
    """The OPERATIVE label set.

    `prepare` measures positive/negative support over the frozen pool and writes
    `functional_groups.frozen.json`, applying min_positives / min_negatives.
    That file is authoritative. Falling back to the candidate file only happens
    before `prepare` has ever run.

    This indirection is the fix for a real defect: the frozen file used to be
    written and never read, so the declared support rule did not actually gate
    which labels the model was asked about.
    """
    if paths.FUNCTIONAL_GROUPS_FROZEN.is_file():
        with open(paths.FUNCTIONAL_GROUPS_FROZEN, encoding='utf-8') as handle:
            return json.load(handle)
    return load_candidate_functional_groups()


def enabled_labels() -> list:
    return [item for item in load_functional_groups()['labels'] if item.get('enabled')]


def frozen_macro_f1_labels(paths_module=None):
    """The label set the FROZEN rule admits, or None when it was never recorded.

    The support rule is applied once, before the run, over the frozen 200-item
    sample. Scoring must average over THAT set rather than re-deriving
    'exercised' labels from whatever positives happen to appear in the sample:
    re-deriving lets a label with 7 positives back in and changes the mean's
    denominator from sample to sample.
    """
    from . import paths as _paths
    module = paths_module or _paths
    manifest = module.DELIV_DIR / 'main_manifest.json'
    if not manifest.is_file():
        return None
    try:
        data = json.loads(manifest.read_text(encoding='utf-8'))
    except (OSError, ValueError):
        return None
    frozen = (data.get('scoring') or {}).get('label_rule') or {}
    names = frozen.get('label_set')
    return sorted(names) if names else None


def label_set_fingerprint() -> str:
    """Stable identity of the label set that is currently operative."""
    names = sorted(item['name'] for item in enabled_labels())
    smarts = {item['name']: item['smarts'] for item in load_functional_groups()['labels']}
    blob = json.dumps({'names': names, 'smarts': [smarts[n] for n in names]},
                      sort_keys=True)
    return hashlib.sha256(blob.encode('utf-8')).hexdigest()


def compile_smarts(name: str, smarts: str):
    pattern = Chem.MolFromSmarts(smarts)
    if pattern is None:
        raise ValueError('SMARTS for %s does not compile: %r' % (name, smarts))
    return pattern


def functional_group_labels(smiles: str, only_enabled: bool = False,
                            candidates: bool = False) -> dict:
    """Match configured SMARTS against the reference structure.

    ``only_enabled=True`` restricts to the operative label set that the task
    actually asks about. ``candidates=True`` uses the raw candidate list, which
    is what support measurement needs so that a label disabled for lack of
    support is still reported with its measured counts rather than as a zero.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError('RDKit cannot parse %r' % smiles)
    config = (load_candidate_functional_groups() if candidates
              else load_functional_groups())
    out = {}
    for item in config['labels']:
        if only_enabled and not item.get('enabled'):
            continue
        out[item['name']] = bool(mol.HasSubstructMatch(compile_smarts(item['name'], item['smarts'])))
    return out


def atom_count(smiles: str, element: str) -> int:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError('RDKit cannot parse %r' % smiles)
    return sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == element)


def atom_counts(smiles: str) -> dict:
    return {element: atom_count(smiles, element) for element in COUNT_ELEMENTS}


def qed_reference(smiles: str) -> float:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError('RDKit cannot parse %r' % smiles)
    return float(QED.qed(mol))


def descriptor_reference(smiles: str) -> dict:
    """Reference values for T5. TPSA is in square angstroms (A^2)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError('RDKit cannot parse %r' % smiles)
    return {
        'hbd': int(rdMolDescriptors.CalcNumHBD(mol)),
        'hba': int(rdMolDescriptors.CalcNumHBA(mol)),
        'logp': float(Crippen.MolLogP(mol)),
        'tpsa': float(rdMolDescriptors.CalcTPSA(mol)),
        '_units': {'hbd': 'count', 'hba': 'count', 'logp': 'dimensionless (Crippen)',
                   'tpsa': 'angstrom^2'},
    }


def is_same_molecule(smiles_a: str, smiles_b: str) -> bool:
    """Identity gold label: same marked molecular graph, judged by RDKit.

    Defined strictly on the common support of the pool: heavy-atom skeleton,
    bond orders, formal charges, aromaticity and (absent, because the pool is
    stereo-free) stereochemistry. Hydrogens are implicit and follow RDKit's
    standard valence model on both sides.
    """
    mol_a = Chem.MolFromSmiles(smiles_a)
    mol_b = Chem.MolFromSmiles(smiles_b)
    if mol_a is None or mol_b is None:
        raise ValueError('RDKit cannot parse one of %r / %r' % (smiles_a, smiles_b))
    return Chem.MolToSmiles(mol_a) == Chem.MolToSmiles(mol_b)


def heavy_atom_count(smiles: str) -> int:
    return Chem.MolFromSmiles(smiles).GetNumHeavyAtoms()


def molecular_formula(smiles: str) -> str:
    return rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(smiles))


__all__ = [
    'COUNT_ELEMENTS', 'DESCRIPTORS_FULL', 'DESCRIPTORS_ROUND1',
    'rdkit_version', 'qed_weights', 'load_functional_groups', 'enabled_labels',
    'functional_group_labels', 'atom_count', 'atom_counts', 'qed_reference',
    'descriptor_reference', 'is_same_molecule', 'heavy_atom_count',
    'molecular_formula', 'Descriptors',
]