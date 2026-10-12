"""ChEBI-20 -> frozen molecule pool -> dev/test splits -> dataset manifest.

Filtering rule (recorded, applied before any model call)
-------------------------------------------------------
A molecule enters the pool only if ALL of the following hold. Nothing is
dropped later to make a representation look better.

  F1  RDKit parses the source SMILES with sanitisation.
  F2  The source SMILES carries no stereochemistry descriptor (@, /, \\). The
      professor prototype declares stereochemistry out of scope and refuses it
      rather than discarding it, so stripping stereo would change the molecule.
      Those rows are excluded, not rewritten.
  F3  Our Code parses it.
  F4  Our Code round-trips it: write -> re-parse gives the same atom/bond count.
  F5  RDKit reads Our Code's output and agrees it is the same molecule as the
      reference (canonical SMILES equality).
  F6  RDKit's canonical SMILES for the molecule is also readable by Our Code,
      so the SMILES and Our Code arms cover exactly the same molecule list.

Splitting is molecule-level: a canonical SMILES lives wholly inside dev or test.
Identity pair endpoints and every representation variant of one molecule stay
in the same split.
"""
from __future__ import annotations

import hashlib
import json
import random
from collections import Counter

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors

from . import labels, paths, reps
from .core_loader import load_core

RDLogger.DisableLog('rdApp.*')

SOURCE_FILES = {
    'validation': 'ChEBI-20_validation.txt',
    'test': 'ChEBI-20_test.txt',
}


def _sha256(path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_source_rows():
    """Yield (source_split, cid, smiles, description) from the raw MolT5 files."""
    for split, name in SOURCE_FILES.items():
        path = paths.RAW_DIR / name
        lines = path.read_text(encoding='utf-8').splitlines()
        if not lines or not lines[0].startswith('CID'):
            raise ValueError('%s does not look like the MolT5 ChEBI-20 table' % name)
        for line in lines[1:]:
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            yield split, parts[0], parts[1], parts[2]


def build_pool():
    """Apply F1-F6 and return (pool, funnel, exclusions)."""
    core = load_core()
    funnel = Counter()
    exclusions = Counter()
    by_canonical = {}

    for source_split, cid, raw_smiles, description in read_source_rows():
        funnel['rows_read'] += 1
        mol = Chem.MolFromSmiles(raw_smiles)
        if mol is None:
            exclusions['F1_rdkit_parse_failed'] += 1
            continue
        funnel['F1_rdkit_parsed'] += 1

        if any(marker in raw_smiles for marker in ('@', '/', chr(92))):
            exclusions['F2_has_stereochemistry'] += 1
            continue
        funnel['F2_no_stereo'] += 1

        try:
            core.Molecule(raw_smiles)
        except Exception:  # noqa: BLE001
            exclusions['F3_our_code_parse_failed'] += 1
            continue
        funnel['F3_our_code_parsed'] += 1

        canonical = reps.canonical_smiles(mol)
        try:
            code = core.Molecule(raw_smiles).decode_to_smiles()
            reread = core.Molecule(code)
        except Exception:  # noqa: BLE001
            exclusions['F4_our_code_roundtrip_failed'] += 1
            continue
        original = core.Molecule(raw_smiles)
        if (len(reread.atoms), len(reread.bonds)) != (len(original.atoms), len(original.bonds)):
            exclusions['F4_our_code_roundtrip_shape'] += 1
            continue
        funnel['F4_our_code_roundtrip'] += 1

        out_mol = Chem.MolFromSmiles(code)
        if out_mol is None or reps.canonical_smiles(out_mol) != canonical:
            exclusions['F5_rdkit_disagrees_with_our_code'] += 1
            continue
        funnel['F5_rdkit_agrees'] += 1

        try:
            core.Molecule(canonical)
        except Exception:  # noqa: BLE001
            exclusions['F6_our_code_cannot_read_canonical'] += 1
            continue
        funnel['F6_canonical_readable'] += 1

        # F7 the canonical string must be a fixed point, so the reference is
        # unambiguous and identical on every arm
        reread_mol = Chem.MolFromSmiles(canonical)
        if reread_mol is None or reps.canonical_smiles(reread_mol) != canonical:
            exclusions['F7_canonical_not_idempotent'] += 1
            continue
        funnel['F7_canonical_idempotent'] += 1

        # F8 every round-1 representation must encode it, so the three arms
        # share exactly one molecule list rather than three overlapping ones
        try:
            reps.selfies_text(canonical)
        except Exception:  # noqa: BLE001
            exclusions['F8_selfies_cannot_encode'] += 1
            continue
        funnel['F8_selfies_encodable'] += 1

        # deduplicate: one entry per molecule, first citation wins
        if canonical in by_canonical:
            exclusions['duplicate_molecule'] += 1
            continue
        try:
            entry = {
                'cid': cid,
                'smiles': canonical,
                'description': description,
                'formula': labels.molecular_formula(canonical),
                'heavy_atoms': labels.heavy_atom_count(canonical),
                'source_split': source_split,
                'cluster': canonical,
            }
        except Exception:  # noqa: BLE001
            exclusions['label_computation_failed'] += 1
            continue
        by_canonical[canonical] = entry

    pool = sorted(by_canonical.values(), key=lambda item: (item['heavy_atoms'], item['smiles']))
    return pool, funnel, exclusions


def split_pool(pool, dev_fraction: float, seed: int):
    """Molecule-level split, then every variant inherits its molecule's split."""
    ordered = sorted(pool, key=lambda item: item['cluster'])
    rng = random.Random(seed)
    shuffled = ordered[:]
    rng.shuffle(shuffled)
    cut = int(round(len(shuffled) * dev_fraction))
    dev_keys = {item['cluster'] for item in shuffled[:cut]}
    dev = [item for item in ordered if item['cluster'] in dev_keys]
    test = [item for item in ordered if item['cluster'] not in dev_keys]
    return dev, test


def task_support(pool):
    """Positive/negative support per functional-group label and per element."""
    support = {'functional_groups': {}, 'atom_counts': {}, 'qed': {}, 'hbd': {}}
    for item in pool:
        group_labels = labels.functional_group_labels(item['smiles'], candidates=True)
        for name, present in group_labels.items():
            bucket = support['functional_groups'].setdefault(name, {'positive': 0, 'negative': 0})
            bucket['positive' if present else 'negative'] += 1
        counts = labels.atom_counts(item['smiles'])
        for element, value in counts.items():
            support['atom_counts'].setdefault(element, Counter())[value] += 1
        qed = labels.qed_reference(item['smiles'])
        support['qed'].setdefault('min', qed)
        support['qed']['min'] = min(support['qed']['min'], qed)
        support['qed'].setdefault('max', qed)
        support['qed']['max'] = max(support['qed']['max'], qed)
        hbd = labels.descriptor_reference(item['smiles'])['hbd']
        support['hbd'].setdefault('distribution', Counter())[hbd] += 1
    for element, counter in support['atom_counts'].items():
        support['atom_counts'][element] = dict(sorted(counter.items()))
    support['hbd']['distribution'] = dict(sorted(support['hbd']['distribution'].items()))
    support['qed'] = {k: round(v, 4) for k, v in support['qed'].items()}
    return support


def write_prepared(pool, dev, test, funnel, exclusions, manifest, quarantined=()):
    paths.PREPARED_DIR.mkdir(parents=True, exist_ok=True)
    groups = [('pool', pool), ('dev', dev), ('test', test)]
    if quarantined:
        groups.append(('quarantined', list(quarantined)))
    for name, rows in groups:
        with open(paths.prepared('%s.jsonl' % name), 'w', encoding='utf-8') as handle:
            for row in rows:
                handle.write(json.dumps(row, ensure_ascii=False) + '\n')
    with open(paths.prepared('dataset_manifest.json'), 'w', encoding='utf-8') as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
    for name, rows in groups:
        with open(paths.prepared('%s.canonical_smiles.sha256' % name), 'w', encoding='utf-8') as handle:
            digest = hashlib.sha256('\n'.join(r['smiles'] for r in rows).encode('utf-8')).hexdigest()
            handle.write(digest + '\n')


def load_prepared(name: str):
    path = paths.prepared('%s.jsonl' % name)
    with open(path, encoding='utf-8') as handle:
        return [json.loads(line) for line in handle if line.strip()]


def load_manifest():
    with open(paths.prepared('dataset_manifest.json'), encoding='utf-8') as handle:
        return json.load(handle)


# --------------------------------------------------------------------------
# complexity gate (the main-scale sample definition)
# --------------------------------------------------------------------------
#
# The reference paper states its comprehension subset keeps molecules with
# "a molecular weight >=300 Da and at least 2 rings", to exclude trivially
# simple molecules "where all representations may succeed equally". It does
# NOT say which molecular-weight measure, which ring count, or how hydrogens
# are treated, and it publishes no code for this step.
#
# The two choices below are therefore THIS PROJECT'S definitions. They are
# recorded explicitly so the difference from the paper stays visible instead
# of being implied by a shared number.
COMPLEXITY = {
    'min_molecular_weight': 300.0,
    'min_rings': 2,
    'molecular_weight_definition': ('RDKit Descriptors.MolWt -- average molecular weight, '
                                    'hydrogens included'),
    'ring_count_definition': ('RDKit rdMolDescriptors.CalcNumRings -- SSSR ring count '
                              '(rings in the smallest set of smallest rings)'),
    'confirmed_from_paper': [
        'the numeric thresholds: >=300 Da and at least 2 rings',
        'the stated purpose: exclude trivially simple molecules where all representations '
        'may succeed equally (small acyclic fragments, single amino acids)',
        'the sampling scale: 200 molecules for the comprehension tasks',
    ],
    'not_confirmed_from_paper': [
        'which molecular-weight measure (average vs exact / monoisotopic)',
        'which ring count (SSSR vs aromatic-only vs total ring count)',
        'whether hydrogens are included in the weight',
        'the exact RDKit calls, or any code for the cut-off',
    ],
    'consequence': ('a molecule near either boundary could fall on the other side of the '
                    'paper\'s own cut-off. This subset is therefore described as using '
                    'the paper\'s stated thresholds under this project\'s definitions, '
                    'not as reproducing the paper\'s subset.'),
}


def molecular_weight(smiles: str) -> float:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError('RDKit cannot parse %r' % smiles)
    return float(Descriptors.MolWt(mol))


def ring_count(smiles: str) -> int:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError('RDKit cannot parse %r' % smiles)
    return int(rdMolDescriptors.CalcNumRings(mol))


def apply_complexity_filter(rows, min_mw=None, min_rings=None):
    """Split a molecule list by the complexity gate. Returns (kept, dropped)."""
    min_mw = COMPLEXITY['min_molecular_weight'] if min_mw is None else min_mw
    min_rings = COMPLEXITY['min_rings'] if min_rings is None else min_rings
    kept, dropped = [], []
    for row in rows:
        weight = molecular_weight(row['smiles'])
        rings = ring_count(row['smiles'])
        entry = dict(row, molecular_weight=weight, rings=rings)
        if weight >= min_mw and rings >= min_rings:
            kept.append(entry)
        else:
            entry['drop_reason'] = ('below_min_molecular_weight' if weight < min_mw
                                    else 'below_min_rings')
            dropped.append(entry)
    return kept, dropped


def complexity_distribution(rows):
    """Molecular-weight and ring-count histograms for a molecule list."""
    weights = [molecular_weight(r['smiles']) for r in rows]
    rings = [ring_count(r['smiles']) for r in rows]
    bands = [(0, 150), (150, 300), (300, 400), (400, 500), (500, 600), (600, 100000)]
    histogram = {('%d-%d' % (lo, hi) if hi < 1000 else '%d+' % lo):
                 sum(1 for w in weights if lo <= w < hi) for lo, hi in bands}
    return {
        'n': len(rows),
        'molecular_weight': {
            'min': min(weights) if weights else None,
            'max': max(weights) if weights else None,
            'mean': (sum(weights) / len(weights)) if weights else None,
            'bands': histogram,
        },
        'rings': {
            'min': min(rings) if rings else None,
            'max': max(rings) if rings else None,
            'distribution': dict(sorted(Counter(rings).items())),
        },
    }