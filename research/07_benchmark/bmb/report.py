"""Tables rebuilt from saved runs. No model calls, no fabricated numbers.

A cell shows one of:
  a number        measured in the named run
  待测 (not run)  implemented but never executed against a model
  未实现 (N/A)    out of this round's scope
  模型不支持       the model has no image interface
  N/A             the label had too little support to score
A zero is written only when a real zero was measured.
"""
from __future__ import annotations

import csv
import json

from . import labels, paths, reps, tasks

NOT_RUN = '待测'
NOT_IMPLEMENTED = '未实现'
MODEL_UNSUPPORTED = '模型不支持'
NA = 'N/A'

MAIN_COLUMNS = [
    ('t1_identity', 'Identity Acc ↑'),
    ('t2_functional_groups', '子结构 Macro-F1 ↑'),
    ('t3_atom_counts', '元素计数 Acc ↑'),
    ('t4_qed', 'QED MAE ↓'),
    ('t5_hydrogen_bond_donors', 'HBD Acc ↑'),
]

PRIMARY_KEY = {
    't1_identity': 'accuracy',
    't2_functional_groups': 'macro_f1_valid_subset',
    't3_atom_counts': 'exact_match_accuracy',
    't4_qed': 'mae',
    't5_hydrogen_bond_donors': 'hbd_exact_match_accuracy',
}


def _fmt(value):
    if value is None:
        return NOT_RUN
    return '%.3f' % value


def _cell(metrics, representation, task_name, split_name):
    block = metrics.get('%s|%s|%s' % (representation, task_name, split_name))
    if block is None:
        return NOT_RUN
    key = PRIMARY_KEY[task_name]
    if block.get(key) is None:
        return NA
    return _fmt(block[key])


# The eight rows of the final comparison table, in display order. The round-1
# registry listed only six names and omitted CML and MolJSON entirely.
ALL_REPRESENTATION_NAMES_EIGHT = ('smiles', 'selfies', 'our_code', 'img', 'graph2d',
                                  'graph3d', 'cml', 'moljson')


def _lookup(name):
    """Registry entry for a representation, from either registry.

    The five added representations live in ``reps_ext``; the three round-1 ones
    live in ``reps``. Looking only in ``reps`` raised KeyError('cml') the first
    time an extension run tried to write its summary, AFTER all the responses
    had already been obtained and persisted.
    """
    from . import reps_ext
    # reps_ext FIRST. The round-1 registry still carries img / graph2d / graph3d /
    # moljson as `implemented: False` placeholders from before they existed, so
    # preferring it made the summary print "模型不支持" and "未实现" for four
    # representations that had just produced hundreds of real responses.
    if name in reps_ext.RENDERERS:
        return type('ExtRep', (), {
            'name': name,
            'kind': reps_ext.REPRESENTATION_KIND.get(name, 'text'),
            'implemented': True,
            'note': ''})()
    if name in reps.REPRESENTATIONS:
        return reps.REPRESENTATIONS[name]
    return None


def main_table(metrics, representations, split_name):
    rows = []
    for name in representations:
        rep = _lookup(name)
        row = {'表示': name}
        for task_name, column in MAIN_COLUMNS:
            if rep is None:
                row[column] = NOT_IMPLEMENTED
            elif not rep.implemented:
                row[column] = MODEL_UNSUPPORTED if name == 'img' and 'image' in rep.note else NOT_IMPLEMENTED
            else:
                row[column] = _cell(metrics, name, task_name, split_name)
        rows.append(row)
    return rows


def companion_table(metrics, representations, split_name):
    rows = []
    for name in representations:
        for task_name, column in MAIN_COLUMNS:
            block = metrics.get('%s|%s|%s' % (name, task_name, split_name))
            if block is None:
                continue
            ci = block.get('ci95_primary') or {}
            rows.append({
                'representation': name,
                'task': task_name,
                'metric': column,
                'value': block.get(PRIMARY_KEY[task_name]),
                'n_total': block.get('n_total'),
                'n_valid': block.get('n_valid'),
                'n_invalid': block.get('n_invalid'),
                'n_pending': block.get('n_pending'),
                'valid_answer_rate': block.get('valid_answer_rate'),
                'ci95_low': ci.get('low'),
                'ci95_high': ci.get('high'),
                'complete': block.get('complete'),
                'qed_success_rate_all_samples': block.get('success_rate_all_samples'),
                'qed_success_threshold': block.get('success_threshold'),
                'qed_spearman': block.get('spearman'),
                'macro_f1_all_samples_penalised': block.get('macro_f1_all_samples_penalised'),
            })
    return rows


def write_summary_csv(run_id, metrics, representations, split_name):
    path = paths.run_dir(run_id) / 'summary.csv'
    rows = main_table(metrics, representations, split_name)
    with open(path, 'w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['表示'] + [c for _, c in MAIN_COLUMNS])
        writer.writeheader()
        writer.writerows(rows)
    comp = paths.run_dir(run_id) / 'companion.csv'
    companion = companion_table(metrics, representations, split_name)
    with open(comp, 'w', encoding='utf-8-sig', newline='') as handle:
        if companion:
            writer = csv.DictWriter(handle, fieldnames=list(companion[0].keys()))
            writer.writeheader()
            writer.writerows(companion)
    return path, comp


def progress_matrix(metrics=None, runs=None, scale=None):
    """representation x task: implemented? run? measured?

    ``metrics`` maps '<rep>|<task>|<split>' to a metric block for one run;
    ``runs`` maps (rep, task) to the run id that produced it. Without them every
    implemented cell is 'not_run', which is the honest state before any call.
    """
    metrics = metrics or {}
    runs = runs or {}
    rows = []
    for name in ALL_REPRESENTATION_NAMES_EIGHT:
        rep = _lookup(name)
        for task_name in tasks.TASK_ORDER:
            run_id = runs.get((name, task_name))
            block = metrics.get('%s|%s|test' % (name, task_name)) or \
                metrics.get('%s|%s|dev' % (name, task_name))
            if rep is None or not rep.implemented:
                state, detail = 'not_implemented', (rep.note if rep else 'unknown representation')
            elif block is None:
                state, detail = 'not_run', 'no model call recorded'
            elif block.get('n_pending', 0) > block.get('n_total', 0) * 0.5:
                state, detail = 'incomplete', 'most items did not return'
            else:
                state, detail = 'measured', 'value present in summary.csv'
            rows.append({
                'representation': name,
                'task': task_name,
                'task_cn': tasks.TASK_TITLES_CN[task_name],
                'implemented': 'yes' if rep.implemented else 'no',
                'run_state': state,
                # 'measured' alone hides how many items it was measured on, and
                # a smoke cell must never be read as a main-scale result
                'scale': scale_for(run_id) if block is not None else '',
                'n_items': block.get('n_total') if block else '',
                'run_id': run_id or '',
                'detail': detail,
            })
    return rows


def scale_for(run_id):
    if not run_id:
        return ''
    prefix = run_id.split('-')[0]
    if prefix == 'smoke':
        return 'smoke (pipeline check only)'
    if prefix == 'main':
        return 'main (formal)'
    if prefix == 'probe':
        return 'probe (single cell)'
    return 'unknown'


def write_progress_matrix(metrics=None, runs=None):
    paths.DELIV_DIR.mkdir(parents=True, exist_ok=True)
    path = paths.DELIV_DIR / 'progress_matrix.csv'
    rows = progress_matrix(metrics, runs)
    with open(path, 'w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return path


def run_state_map(run_id, representations, task_names):
    """(representation, task) -> run_id, for every cell a run was asked to cover."""
    return {(rep, task): run_id for rep in representations for task in task_names}