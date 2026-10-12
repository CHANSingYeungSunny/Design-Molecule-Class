"""Recompute every frozen reference label with RDKit and compare it to the stored one.

This is reproduction path B: it does not call a model and it does not reuse the
saved answer. It re-derives the reference (oracle) values from the same frozen
molecule strings with the project's own oracle module, then compares them with
the labels that were stored next to the frozen items.

A pass here means: the published reference answers are reproducible from the
frozen data plus RDKit. It says nothing about whether the model was right.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / 'research/07_benchmark'
DELIV = ROOT / 'research/08_deliverables'
sys.path.insert(0, str(BENCH))

from bmb import labels  # noqa: E402

REPS = ('smiles', 'selfies', 'our_code', 'img', 'graph2d', 'graph3d', 'cml', 'moljson')
TASKS = ('t1_identity', 't2_functional_groups', 't3_atom_counts', 't4_qed',
         't5_hydrogen_bond_donors')
QED_TOLERANCE = 1e-12


def rows(path):
    with path.open(encoding='utf-8') as stream:
        for line in stream:
            if line.strip():
                yield json.loads(line)


def pair_key(example_id, representation):
    return example_id.replace('-' + representation + '-', '-', 1)


def load_frozen():
    main_manifest = json.loads((DELIV / 'main_manifest.json').read_text(encoding='utf-8'))
    extension_manifest = json.loads((DELIV / 'extension_manifest.json').read_text(encoding='utf-8'))
    frozen = {}
    for rep in REPS:
        for task in TASKS:
            if rep in REPS[:3]:
                path = BENCH / 'data/prepared/main_items' / ('%s-%s.jsonl' % (rep, task))
            else:
                entry = extension_manifest['item_lists'][task][rep]
                path = BENCH / 'data/prepared' / entry['pool']
            items = list(rows(path))
            # example_id differs per representation; the pair key is the shared
            # question identity across the eight arms.
            by_id = {}
            for item in items:
                by_id[pair_key(item['example_id'], rep)] = item
            frozen[(rep, task)] = by_id
    return frozen, main_manifest, extension_manifest


def compare_functional_groups(item):
    smiles = item['cluster']
    recomputed = labels.functional_group_labels(smiles, only_enabled=True)
    stored = item['label']
    missing = sorted(set(recomputed) - set(stored))
    differing = sorted(name for name in recomputed
                       if name in stored and bool(stored[name]) != recomputed[name])
    return {'stored_keys': len(stored), 'recomputed_keys': len(recomputed),
            'missing_in_stored': missing, 'differing': differing}


def compare_atom_counts(item):
    smiles = item['cluster']
    recomputed = labels.atom_counts(smiles)
    stored = item['label']
    differing = {element: [stored.get(element), recomputed[element]]
                 for element in recomputed if stored.get(element) != recomputed[element]}
    return {'recomputed': recomputed, 'differing': differing}


def compare_qed(item):
    recomputed = labels.qed_reference(item['cluster'])
    stored = float(item['label']['qed'])
    delta = abs(recomputed - stored)
    return {'stored': stored, 'recomputed': recomputed, 'delta': delta,
            'within_tolerance': delta <= QED_TOLERANCE}


def compare_hbd(item):
    recomputed = labels.descriptor_reference(item['cluster'])['hbd']
    stored = int(item['label']['hbd'])
    return {'stored': stored, 'recomputed': recomputed, 'matching': stored == recomputed}


def compare_identity(item):
    recomputed = labels.is_same_molecule(*item['clusters'])
    stored = bool(item['label']['same_molecule'])
    return {'stored': stored, 'recomputed': recomputed, 'matching': stored == recomputed}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT / 'reproduced')
    args = parser.parse_args()
    frozen, _main, _ext = load_frozen()

    per_task = {}
    per_item = []
    failures = []
    for task in TASKS:
        common = sorted(set.intersection(*(set(frozen[(rep, task)]) for rep in REPS)))
        checked = mismatched = 0
        for key in common:
            group = [frozen[(rep, task)][key] for rep in REPS]
            reference = group[0]
            for rep, item in zip(REPS, group):
                if item['cluster'] != reference['cluster']:
                    failures.append('%s: molecule drift across arms (%s)' % (key, rep))
                if json.dumps(item['label'], sort_keys=True) != json.dumps(
                        reference['label'], sort_keys=True):
                    failures.append('%s: label drift across arms (%s)' % (key, rep))
            if task == TASKS[0]:
                result = compare_identity(reference)
                ok = result['matching']
            elif task == TASKS[1]:
                result = compare_functional_groups(reference)
                ok = not result['differing'] and not result['missing_in_stored']
            elif task == TASKS[2]:
                result = compare_atom_counts(reference)
                ok = not result['differing']
            elif task == TASKS[3]:
                result = compare_qed(reference)
                ok = result['within_tolerance']
            else:
                result = compare_hbd(reference)
                ok = result['matching']
            checked += 1
            mismatched += 0 if ok else 1
            if not ok:
                failures.append('%s: recomputed label differs: %s' % (key, result))
            per_item.append({'task': task, 'pair_key': key, 'molecule': reference['cluster'],
                             'stored_label': reference['label'], 'recomputed': result,
                             'matching': ok})
        per_task[task] = {'checked': checked, 'mismatched': mismatched}

    args.out.mkdir(parents=True, exist_ok=True)
    report = {
        'mode': 'reference-label rebuild with RDKit; no model requests',
        'rdkit_version': labels.rdkit_version(),
        'python_version': sys.version.split()[0],
        'qed_tolerance': QED_TOLERANCE,
        'functional_groups_asked': sorted(entry['name'] for entry in labels.enabled_labels()),
        'functional_groups_scored': labels.frozen_macro_f1_labels(),
        'per_task': per_task,
        'total_checked': sum(v['checked'] for v in per_task.values()),
        'total_mismatched': sum(v['mismatched'] for v in per_task.values()),
        'all_labels_reproduced': not failures,
        'failures': failures[:20],
    }
    (args.out / 'label_rebuild.json').write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    with (args.out / 'label_rebuild_items.jsonl').open('w', encoding='utf-8') as stream:
        for row in per_item:
            stream.write(json.dumps(row, ensure_ascii=False) + '\n')
    for task in TASKS:
        print('%-24s checked %4d  mismatched %d'
              % (task, per_task[task]['checked'], per_task[task]['mismatched']))
    print('RDKit', labels.rdkit_version(), '| total', report['total_checked'],
          'mismatched', report['total_mismatched'])
    print('Outputs:', args.out)
    return 1 if failures else 0


if __name__ == '__main__':
    raise SystemExit(main())
