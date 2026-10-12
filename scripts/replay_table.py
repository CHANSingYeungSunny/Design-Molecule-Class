"""Recompute the published 8 x 5 table from saved responses; standard library only.

This is historical replay, not new model inference or regeneration of RDKit labels.
The original scoring.py supplies the metrics. Frozen lists, raw JSON, prompts,
labels, image bytes, item pairing and the expected table are independently checked.
"""
from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import importlib.util
import json
from pathlib import Path, PureWindowsPath

ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / 'research/07_benchmark'
DELIV = ROOT / 'research/08_deliverables'
REPS = ('smiles', 'selfies', 'our_code', 'img', 'graph2d', 'graph3d', 'cml', 'moljson')
TASKS = ('t1_identity', 't2_functional_groups', 't3_atom_counts', 't4_qed',
         't5_hydrogen_bond_donors')
CORE_SHA = '55d314aa312707b498536d339ff59dc89bbbe4956cab54a74a4746cd385af57b'
RAW_SHA = {
    'main-02': 'b7811e016f7e7e1e44f17182a54794aeb08efe5895d836c74e823de1ef3df282',
    'extension-01': 'e15eb4fef2b21bb1a9fe808396fa7965beabaa6180ed8869ee976e597edd1062',
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8'))


def rows(path):
    with path.open(encoding='utf-8') as stream:
        for line in stream:
            if line.strip():
                yield json.loads(line)


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_list_hash(items):
    text = '\n'.join(json.dumps(row, ensure_ascii=False, sort_keys=True) for row in items)
    return hashlib.sha256(text.encode('utf-8')).hexdigest()


def pair_key(record):
    return record['example_id'].replace('-' + record['representation'] + '-', '-', 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT / 'reproduced')
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location('historical_scoring', BENCH / 'bmb/scoring.py')
    scoring = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scoring)
    for path in [ROOT / 'molecule_standalone.py',
                 ROOT / 'research/00_baseline/molecule_standalone.py',
                 ROOT / 'research/01_implementation/molecule_core.py']:
        require(sha(path) == CORE_SHA, f'Baseline changed: {path}')

    main_manifest = read_json(DELIV / 'main_manifest.json')
    extension_manifest = read_json(DELIV / 'extension_manifest.json')
    label_set = main_manifest['scoring']['label_rule']['label_set']
    require(len(label_set) == 10, 'Expected the historical fixed 10-label set')
    system = (ROOT / 'prompts/system.txt').read_text(encoding='utf-8').rstrip('\n')
    frozen, image_paths = {}, {}
    image_checks = 0
    for rep in REPS:
        for task in TASKS:
            if rep in REPS[:3]:
                path = BENCH / f'data/prepared/main_items/{rep}-{task}.jsonl'
                expected = main_manifest['item_list_hashes'][task][rep]
                allowed_prompts = main_manifest['prompt_hashes'][rep][task]
            else:
                entry = extension_manifest['item_lists'][task][rep]
                path = BENCH / 'data/prepared' / entry['pool']
                expected = entry['sha256']
                allowed_prompts = extension_manifest['prompt_hashes'][rep][task]
            items = list(rows(path))
            require(canonical_list_hash(items) == expected,
                    f'Frozen content hash mismatch: {rep}/{task}')
            by_id = {}
            for item in items:
                require(item['example_id'] not in by_id, f'Duplicate frozen id: {item["example_id"]}')
                require((item['rep'], item['task'], item['split']) == (rep, task, 'test'),
                        f'Wrong frozen cell: {item["example_id"]}')
                by_id[item['example_id']] = item
                for slot, value in item['payload'].items():
                    if isinstance(value, dict) and 'image_path' in value:
                        image_path = BENCH.joinpath(*PureWindowsPath(value['image_path']).parts)
                        require(sha(image_path) == value['sha256'], f'Image changed: {image_path}')
                        image_paths[(item['example_id'], slot)] = image_path
                        image_checks += 1
            frozen[(rep, task)] = (by_id, set(allowed_prompts))

    latest, counts = {}, {}
    for run in RAW_SHA:
        path = BENCH / 'runs' / run / 'raw_responses.jsonl'
        require(sha(path) == RAW_SHA[run], f'Historical response file changed: {run}')
        count = ignored = 0
        for record in rows(path):
            count += 1
            cell = (record['representation'], record['task'])
            ids, allowed = frozen[cell]
            if record['example_id'] not in ids:
                ignored += 1
                continue
            item = ids[record['example_id']]
            require(record['label'] == item['label'], f'Label drift: {record["example_id"]}')
            require(record['cluster'] == item['cluster'], f'Molecule drift: {record["example_id"]}')
            if record['status'] == 'ok':
                parsed = json.loads(record['raw_response'])
                require(all(k in parsed and parsed[k] == v for k, v in record['parsed'].items()),
                        f'Saved parsed answer disagrees with raw JSON: {record["example_id"]}')
                # The replay uses freshly parsed raw JSON, not cached per-item scores.
                record['parsed'] = {k: parsed[k] for k in record['parsed']}
                messages = [{'role': 'system', 'content': system},
                            {'role': 'user', 'content': record['prompt_text']}]
                prompt_hash = hashlib.sha256(json.dumps(messages, ensure_ascii=False,
                                                        sort_keys=True).encode('utf-8')).hexdigest()
                require(prompt_hash == record['prompt_hash'] and prompt_hash in allowed,
                        f'Prompt drift: {record["example_id"]}')
                if cell[0] == 'img':
                    blocks = record['prompt_text']
                    image_blocks = [x for x in blocks if x['type'] == 'image_url']
                    require(len(image_blocks) == len(item['payload']), 'Image slot count mismatch')
                    for slot, block in zip(item['payload'], image_blocks):
                        blob = base64.b64decode(block['image_url']['url'].split(',', 1)[1], validate=True)
                        require(hashlib.sha256(blob).hexdigest() == item['payload'][slot]['sha256'],
                                f'Request image differs from frozen image: {record["example_id"]}')
                else:
                    for value in item['payload'].values():
                        require(value in record['prompt_text'], f'Frozen payload missing from prompt: {record["example_id"]}')
            latest[(cell, record['example_id'])] = record
        counts[run] = {'raw_rows': count, 'superseded_rows_outside_current_lists': ignored}

    cells, coverage = {}, {}
    for (rep, task), (ids, _) in frozen.items():
        rr = [latest[((rep, task), identifier)] for identifier in ids
              if ((rep, task), identifier) in latest]
        require(len(rr) == len(ids), f'Missing scheduled answers: {rep}/{task}')
        require(all(r['status'] == 'ok' for r in rr), f'Incomplete final cell: {rep}/{task}')
        keyed = {pair_key(r): r for r in rr}
        require(len(keyed) == len(rr), f'Pair-key collision: {rep}/{task}')
        cells[(rep, task)] = keyed
        coverage[f'{rep}|{task}'] = {'planned': len(ids), 'valid': len(rr)}

    expected = read_json(DELIV / 'analysis_common_items.json')
    results, denominators, per_item = {}, {}, []
    for task in TASKS:
        common = sorted(set.intersection(*(set(cells[(rep, task)]) for rep in REPS)))
        require(len(common) == expected['denominators'][task], f'Wrong denominator: {task}')
        denominators[task] = len(common)
        for key in common:
            group = [cells[(rep, task)][key] for rep in REPS]
            require(all(r['label'] == group[0]['label'] for r in group), f'Cross-arm labels differ: {key}')
            require(all(r['cluster'] == group[0]['cluster'] for r in group), f'Cross-arm molecule differs: {key}')
            if task == TASKS[0]:
                require(all(r['clusters'] == group[0]['clusters'] for r in group), f'Identity endpoints differ: {key}')
        for rep in REPS:
            rr = [cells[(rep, task)][key] for key in common]
            value = scoring._primary_value(task, rr, label_set, 0.05)
            target = expected['per_task'][task]['values'][rep]['value']
            require(value is not None and abs(value - target) < 1e-12,
                    f'Replayed value differs from reference: {rep}/{task}: {value} vs {target}')
            results.setdefault(rep, {})[task] = value
            for record in rr:
                per_item.append({'representation': rep, 'task': task, 'example_id': record['example_id'],
                                 'score_or_error': scoring.score_item(record),
                                 'parsed_from_raw': record['parsed'], 'reference': record['label']})

    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out / 'summary.csv').open('w', encoding='utf-8', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(['representation'] + list(TASKS))
        writer.writerows([rep] + [results[rep][task] for task in TASKS] for rep in REPS)
    report = {'mode': 'historical offline replay; no model requests; stored RDKit labels',
              'raw_counts': counts, 'n_frozen_lists': len(frozen), 'n_image_slots_verified': image_checks,
              'coverage': coverage, 'denominators': denominators, 'scoring_label_set': label_set,
              'all_40_values_match': True, 'results': results}
    (args.out / 'verification.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    with (args.out / 'per_item.jsonl').open('w', encoding='utf-8') as stream:
        for row in per_item: stream.write(json.dumps(row, ensure_ascii=False) + '\n')
    print('Historical replay passed: 40 frozen lists, all 40 table values match.')
    print('Denominators:', denominators)
    for rep in REPS: print(rep, ' '.join(f'{results[rep][task]:.3f}' for task in TASKS))
    print('Outputs:', args.out)


if __name__ == '__main__':
    main()
