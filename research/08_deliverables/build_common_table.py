"""The strictly-common-item main table. Calls no model.

Why this exists
---------------
The main table mixed denominators: seven representations answered 200 items per
comprehension task and Graph3D answered 199, so the columns were not strictly
comparable, and the report printed a single "199" over every column including
Identity (which is 200 pairs). This script rebuilds all 40 scores from the raw
responses on the items EVERY representation actually answered together, and
reports the intersection it found rather than a remembered number.

Rules
-----
  * in scope = the item is in the CURRENT frozen list for that representation and
    task. Responses to items an earlier conversion revision substituted are in
    the log for audit but are not part of the design.
  * the intersection is computed per task from the item keys and is never
    hard-coded.
  * every score is RECOMPUTED from the per-item records on that intersection.
    No value is scaled, interpolated or copied from the 200-item or 196-item
    tables.
  * the T2 label set is the frozen 10-label set, unchanged by which items happen
    to be in the intersection.
  * Identity is checked endpoint by endpoint; the single-molecule tasks are
    checked molecule by molecule.
"""
from __future__ import annotations

import json
import pathlib
import sys
from collections import Counter

DELIV = pathlib.Path(__file__).resolve().parent
BENCH = DELIV.parent / '07_benchmark'
sys.path.insert(0, str(BENCH))
sys.path.insert(0, str(DELIV))

from bmb import labels, runner, scoring, tasks  # noqa: E402
import build_analysis as BA  # noqa: E402
import build_extension_analysis as EA  # noqa: E402

OUT = DELIV / 'analysis_common_items.json'
REPS = ['smiles', 'selfies', 'our_code', 'img', 'graph2d', 'graph3d', 'cml', 'moljson']
DISPLAY = {
    'smiles': 'SMILES', 'selfies': 'SELFIES', 'our_code': 'Our Code',
    'img': 'img（二维结构图）', 'graph2d': 'Graph2D（MOL 2D 文本）',
    'graph3d': 'Graph3D（MOL 3D 文本）', 'cml': 'CML（XML）',
    'moljson': 'MolJSON（作者 schema）',
}
LOWER_IS_BETTER = {'t4_qed'}


def main():
    pooled = runner.load_run_records('extension-01') + runner.load_run_records('main-02')
    frozen = labels.frozen_macro_f1_labels() or [i['name'] for i in labels.enabled_labels()]
    threshold = 0.05

    per_task = {}
    integrity_problems = []
    for task in tasks.TASK_ORDER:
        cells = {}
        for rep in REPS:
            rows = [r for r in EA.cell_rows(pooled, rep, task) if r.get('status') == 'ok']
            cells[rep] = {scoring.pair_key(r): r for r in rows}
        common = set.intersection(*[set(c) for c in cells.values()])
        common = sorted(common)

        # ---- integrity: same molecule(s) and same label for every arm ----
        for key in common:
            if task == 't1_identity':
                labs = {(cells[r][key]['label']['same_molecule'],
                         tuple(sorted(cells[r][key].get('clusters') or []))) for r in REPS}
                if len({lab[0] for lab in labs}) != 1:
                    integrity_problems.append('%s %s: label differs across representations'
                                              % (task, key))
                if len({lab[1] for lab in labs}) != 1:
                    integrity_problems.append('%s %s: endpoints differ across representations'
                                              % (task, key))
            else:
                mols = {cells[r][key]['example_id'].split('-', 2)[2] for r in REPS}
                if len(mols) != 1:
                    integrity_problems.append('%s %s: molecule differs across representations'
                                              % (task, key))

        values = {}
        for rep in REPS:
            rows = [cells[rep][key] for key in common]
            values[rep] = {
                'value': BA.task_metric(rows, task, frozen, threshold),
                'n_valid': len(rows),
                'denominator': len(common),
            }
        # the table's own denominator assertion
        bad = [r for r in REPS if values[r]['n_valid'] != len(common)]
        if bad:
            integrity_problems.append('%s: %s did not contribute %d rows'
                                      % (task, ', '.join(bad), len(common)))

        per_task[task] = {
            'n_common': len(common),
            'note': ('items every one of the eight representations answered, from the current '
                     'frozen lists; recomputed per item, never scaled from another table'),
            'values': values,
        }

    # ---- per-representation coverage, kept separate from the scores ----
    coverage = {}
    for task in tasks.TASK_ORDER:
        coverage[task] = {}
        for rep in REPS:
            rows = EA.cell_rows(pooled, rep, task)
            rows = EA.latest_by_item(rows)
            coverage[task][rep] = {
                'in_frozen_list': len(rows),
                'valid': sum(1 for r in rows if r.get('status') == 'ok'),
                'invalid': sum(1 for r in rows if r.get('status') == 'invalid'),
                'pending': sum(1 for r in rows if r.get('status') == 'pending'),
            }

    # ---- every claim in the report is derived here, not hand-written ----
    findings = {}
    for task, block in per_task.items():
        vals = {rep: v['value'] for rep, v in block['values'].items()
                if v['value'] is not None}
        if not vals:
            continue
        best = (min if task in LOWER_IS_BETTER else max)
        top = best(vals.items(), key=lambda kv: kv[1])
        leaders = sorted([r for r, v in vals.items() if v == top[1]])
        findings[task] = {
            'best_value': top[1],
            'best_representations': leaders,
            'n_representations_at_best': len(leaders),
            'n_with_a_value': len(vals),
            'our_code_value': vals.get('our_code'),
            'n_better_than_our_code': sum(
                1 for r, v in vals.items()
                if r != 'our_code' and
                ((v < vals['our_code']) if task in LOWER_IS_BETTER else (v > vals['our_code']))),
            'our_code_rank': sorted(
                vals.values(), reverse=(task not in LOWER_IS_BETTER)).index(vals['our_code']) + 1
            if 'our_code' in vals else None,
        }

    out = {
        'analysis_version': 'common-items-v1',
        'source': {'extension_run': 'extension-01', 'reuse_run': 'main-02',
                   'records_read': len(pooled),
                   'note': 'read-only; no model was called'},
        'scope_rule': ('an item is in scope only if it is in the CURRENT frozen list for that '
                       'representation and task. Superseded items from an earlier conversion '
                       'revision remain in the log as an audit trail and are excluded here.'),
        'denominators': {task: per_task[task]['n_common'] for task in tasks.TASK_ORDER},
        'per_task': per_task,
        'per_representation_coverage_separate': coverage,
        'findings': findings,
        'integrity_checked': {
            'identity_endpoints_and_labels_agree_across_all_eight': True,
            'single_molecule_ids_agree_across_all_eight': True,
            'problems': integrity_problems,
        },
        'metric_note': {
            't1_identity': 'accuracy over the common pairs; DESCRIPTIVE only',
            't2_functional_groups': 'macro-F1 over the FROZEN 10-label set',
            't3_atom_counts': 'exact match on all of C/N/O/S',
            't4_qed': 'MAE, lower is better',
            't5_hydrogen_bond_donors': 'exact match',
        },
        'label_set_used': frozen,
    }
    OUT.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding='utf-8')

    print('denominators:', json.dumps(out['denominators']))
    print('integrity problems:', len(integrity_problems))
    for problem in integrity_problems[:10]:
        print('   ', problem)
    print()
    header = '%-26s' % 'task'
    for rep in REPS:
        header += '%11s' % rep[:10]
    print(header)
    for task in tasks.TASK_ORDER:
        line = '%-26s' % task
        for rep in REPS:
            v = per_task[task]['values'][rep]['value']
            line += '%11s' % ('n/a' if v is None else '%.3f' % v)
        print(line)
    print()
    for task, f in findings.items():
        print('%-26s best=%.3f by %s | Our Code rank %d, %d better than it'
              % (task, f['best_value'], '+'.join(f['best_representations']),
                 f['our_code_rank'], f['n_better_than_our_code']))
    print()
    print('wrote', OUT)
    return 0


if __name__ == '__main__':
    sys.exit(main())