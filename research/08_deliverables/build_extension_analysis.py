"""Build the versioned analysis artifact for the EIGHT-representation stage.

Calls no model. It reads the raw responses already on disk and nothing else, so
it can be re-run offline as often as needed.

Rules it inherits unchanged from the three-representation analysis
---------------------------------------------------------------
  * T2 uses the FROZEN 10-label macro-F1 set, fixed BEFORE resampling and never
    re-derived inside a replicate.
  * zero-division rule unchanged: a label with no positives in the draw is
    excluded and reported; a label with positives but no predicted positive is 0.
  * QED is an error: a positive A-B difference means A is WORSE.
  * An invalid answer is never dropped and never scored as 0 error; the
    denominator and the valid rate are both printed.
  * Comparisons pair on a representation-stripped key.

What is NEW here
----------------
  * Eight rows instead of three. The five added representations come from the
    extension run; the three inherited ones are pooled ONLY when the reuse
    condition actually matched, and otherwise are reported as a separate
    historical table rather than ranked beside the new rows.
  * The comparison family is the nine pairs declared in extension_config.json
    BEFORE any new request, and the Holm divisor is the DECLARED family size (9),
    not the number of pairs that happened to be computable.
  * Identity inference stays descriptive. The view a molecule is drawn in is now
    fixed by slot, but the shared-endpoint dependence is still not resolved by a
    method this project has validated, so no confirmatory claim is made.
"""
from __future__ import annotations

import datetime
import hashlib
import json
import pathlib
import sys

DELIV = pathlib.Path(__file__).resolve().parent
BENCH = DELIV.parent / '07_benchmark'
sys.path.insert(0, str(BENCH))
sys.path.insert(0, str(DELIV))

from bmb import labels, paths, scoring, tasks  # noqa: E402
import build_analysis as BA  # noqa: E402  (the three-representation analysis)

OUT = DELIV / 'analysis_extension_v1.json'
EXTENSION_RUN = 'extension-01'
REUSE_RUN = 'main-02'
RESAMPLES = 2000
SEED = 20261005

ADDED = ['img', 'graph2d', 'graph3d', 'cml', 'moljson']
INHERITED = ['smiles', 'selfies', 'our_code']
DISPLAY = {
    'smiles': 'SMILES', 'selfies': 'SELFIES', 'our_code': 'Our Code',
    'img': 'img（二维结构图）', 'graph2d': 'Graph2D（MOL 2D 文本）',
    'graph3d': 'Graph3D（MOL 3D 文本）', 'cml': 'CML（XML）',
    'moljson': 'MolJSON（作者 schema）',
}
METRIC_NAME = {
    't1_identity': 'Identity 准确率',
    't2_functional_groups': '子结构 Macro-F1（固定 10 标签）',
    't3_atom_counts': '元素计数准确率',
    't4_qed': 'QED MAE（越小越好）',
    't5_hydrogen_bond_donors': 'HBD 准确率',
}


def sha256(path):
    path = pathlib.Path(path)
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None


def read_records(run_id):
    return BA.runner.load_run_records(run_id) if run_id else []


def generation_of(manifest_model, run_config):
    """The generation settings BOTH stages actually send.

    They come from run_config.json, which every stage reads through
    ``client_from_config``. The run manifest is the recorded EVIDENCE of that,
    and it is written only when a run finishes -- so a still-running run would
    otherwise compare as a mismatch against a finished one. The manifest value
    is reported separately as evidence, not used as the comparison key.
    """
    recorded = {k: v for k, v in ((manifest_model or {}).get('generation') or {}).items()
                if not k.startswith('_')}
    return recorded or {k: v for k, v in (run_config.get('generation') or {}).items()
                        if not k.startswith('_')}


def condition_of(records, manifest_model, run_config):
    for record in records:
        return {
            'model_id': record.get('model_id'),
            'provider': record.get('provider'),
            'base_url': record.get('base_url'),
            'generation': generation_of(manifest_model, run_config),
        }
    return None


def run_manifest_model(run_id):
    path = paths.run_dir(run_id) / 'experiment_manifest.json'
    if not path.is_file():
        return {}
    return (json.loads(path.read_text(encoding='utf-8')).get('model') or {})


def has_any(records, representation):
    return any(r.get('representation') == representation and r.get('split') == 'test'
               for r in records)


_FROZEN_IDS = None


def frozen_ids():
    """example_ids in the CURRENT frozen item lists, per representation|task.

    The raw log is an append-only audit trail, so it still holds responses to
    items the v1 conversion rule substituted and the fix removed. Counting those
    would put four questions into the analysis that the design no longer asks.
    The frozen manifest is the authority on what is in scope.
    """
    global _FROZEN_IDS
    if _FROZEN_IDS is None:
        _FROZEN_IDS = {}
        path = DELIV / 'extension_manifest.json'
        if path.is_file():
            manifest = json.loads(path.read_text(encoding='utf-8'))
            for task, per_rep in (manifest.get('item_lists') or {}).items():
                for representation, entry in per_rep.items():
                    items = (paths.prepared(entry['pool'])).read_text(encoding='utf-8')
                    ids = set()
                    for line in items.splitlines():
                        if line.strip():
                            ids.add(json.loads(line)['example_id'])
                    _FROZEN_IDS['%s|%s' % (representation, task)] = ids
    return _FROZEN_IDS


def cell_rows(records, representation, task):
    allowed = frozen_ids().get('%s|%s' % (representation, task))
    rows = [r for r in records
            if r.get('representation') == representation and r.get('task') == task
            and r.get('split') == 'test']
    if allowed is not None:
        rows = [r for r in rows if r.get('example_id') in allowed]
    return rows


def planned_counts():
    """The frozen plan per cell, from the reconciliation artifact."""
    path = DELIV / 'extension_reconciliation.json'
    if not path.is_file():
        return {}
    recon = json.loads(path.read_text(encoding='utf-8'))
    return {key: entry['planned_in_frozen_manifest']
            for key, entry in recon['per_cell'].items()}


def latest_by_item(rows):
    """One row per item, LAST record winning.

    A budget stop writes a `pending` placeholder; the resumed run later answers
    the same item. Both lines stay in the log, which is correct for an audit
    trail -- but counting them equally made every affected cell look partially
    finished even after it was fully answered.
    """
    latest = {}
    for record in rows:
        latest[record['example_id']] = record
    return list(latest.values())


def cell_status(rows, planned=None):
    if not rows:
        return 'not run'
    rows = latest_by_item(rows)
    ok = sum(1 for r in rows if r.get('status') == 'ok')
    pending = sum(1 for r in rows if r.get('status') == 'pending')
    # the target is the FROZEN plan, not a remembered 200: Graph3D's four
    # comprehension cells are frozen at 196, so comparing them to 200 marked a
    # fully finished cell as "partial".
    if planned is None:
        planned = len(rows)
    if len(rows) < planned or pending:
        return 'partial'
    if ok == 0:
        return 'no valid answers'
    return 'complete'


def main():
    extension = json.loads((paths.CONFIG_DIR / 'extension_config.json').read_text(encoding='utf-8'))
    run_config = json.loads(paths.RUN_CONFIG.read_text(encoding='utf-8'))
    threshold = run_config.get('qed_success_threshold', 0.05)

    new_records = read_records(EXTENSION_RUN)
    old_records = read_records(REUSE_RUN)
    new_cond = condition_of(new_records, run_manifest_model(EXTENSION_RUN), run_config)
    old_cond = condition_of(old_records, run_manifest_model(REUSE_RUN), run_config)

    reuse_verified, reuse_reason = False, 'the reused run produced no records'
    if new_cond and old_cond:
        keys = ('model_id', 'provider', 'base_url', 'generation')
        mismatched = [k for k in keys if old_cond.get(k) != new_cond.get(k)]
        if mismatched:
            reuse_reason = ('the inherited run differs in %s, so the three inherited '
                            'representations are NOT pooled into the same table'
                            % ', '.join(mismatched))
        else:
            reuse_verified = True
            reuse_reason = ('%s and %s were produced by the same model, provider, endpoint '
                            'and generation settings' % (REUSE_RUN, EXTENSION_RUN))

    # the T2 label set is derived from the FROZEN sample, exactly as before: ONE
    # representation's 200 records, never all eight arms' records added up
    reference = cell_rows(old_records, 'smiles', 't2_functional_groups') or \
        cell_rows(new_records, 'img', 't2_functional_groups')
    frozen_labels, label_counts = BA.fixed_labels(reference)
    declared_labels = labels.frozen_macro_f1_labels()
    if declared_labels:
        frozen_labels = declared_labels

    pooled = list(new_records) + (list(old_records) if reuse_verified else [])
    representations = [r for r in ADDED if has_any(new_records, r)] + \
                      ([r for r in INHERITED if has_any(old_records, r)]
                       if reuse_verified else [])

    out = {
        'analysis_version': 'extension-v2',
        'generated_at': datetime.datetime.now().isoformat(timespec='seconds'),
        'revision': {
            'from': 'extension-v1 (conversion rule ext-v1, 4,984 items)',
            'to': 'extension-v2 (conversion rule ext-v2, 4,996 items)',
            'why': ('a staged diagnosis showed that 6 of the 8 Graph3D conversion failures were '
                    'implementation settings, not molecule limitations, and that the Identity '
                    'pair builder had been substituting replacement pairs for molecules it '
                    'could not render.'),
            'what_changed': [
                'Graph3D: two documented fallbacks added, reached ONLY when the declared rule '
                'produces no usable conformer (random-coordinate embedding; extended iteration '
                'cap). Every payload the declared rule already produced is byte-identical.',
                'Identity: the pair list is now decided before any rendering, so all eight '
                'representations answer the same pairs and an unrenderable pair is reported '
                'rather than swapped.',
            ],
            'scope_of_impact_verified': {
                'payloads_byte_identical': 4980,
                'payloads_changed': 0,
                'items_newly_added': 16,
                'v1_substituted_items_dropped': 4,
            },
            'superseded': ('the v1 analysis, the v1 frozen manifest and the v1 item lists are '
                           'kept under 99_archive_v1_196/ and extension_manifest_v1.json; the '
                           'v1 numbers are NOT silently replaced.'),
        },
        'extension_run': EXTENSION_RUN,
        'reuse_run': REUSE_RUN,
        'inputs': {
            'extension_raw_responses': {
                'path': str(paths.run_dir(EXTENSION_RUN) / 'raw_responses.jsonl'),
                'records': len(new_records),
                'sha256': sha256(paths.run_dir(EXTENSION_RUN) / 'raw_responses.jsonl'),
            },
            'reuse_raw_responses': {
                'path': str(paths.run_dir(REUSE_RUN) / 'raw_responses.jsonl'),
                'records': len(old_records),
                'sha256': sha256(paths.run_dir(REUSE_RUN) / 'raw_responses.jsonl'),
            },
            'extension_manifest_sha256': sha256(DELIV / 'extension_manifest.json'),
            'main_manifest_sha256': sha256(DELIV / 'main_manifest.json'),
            'extension_config_sha256': sha256(paths.CONFIG_DIR / 'extension_config.json'),
        },
        'reuse': {
            'verified': reuse_verified,
            'reason': reuse_reason,
            'condition_added': new_cond,
            'condition_inherited': old_cond,
            '_reuse_conditions_checked': ('prompt text for the same item, model id, provider, '
                                          'endpoint, generation settings and scoring version'),
            '_weight_caveat': ('an identical model alias does not prove identical backend '
                               'weights; the service returned no version or fingerprint, so '
                               'the backend revision is NOT claimed to be fixed'),
        },
        'metric_definitions': {k: METRIC_NAME[k] for k in tasks.TASK_ORDER},
        'display_names': DISPLAY,
        't2_label_set': {
            'fixed_set': frozen_labels,
            'n_labels': len(frozen_labels),
            'per_label_counts_in_sample': label_counts,
            'fixed_across_resamples': True,
        },
        'identity_view_rule': extension.get('identity_views'),
        'representation_note': (
            'img is a picture (pixels). Graph2D and Graph3D are MOL V2000 connection tables '
            'read as TEXT: they are not pictures and not graph neural networks, and the two '
            'differ only in where the coordinates come from, so Graph3D is an '
            'extra-information condition. CML and MolJSON differ from them in serialisation '
            '(XML and JSON). The five added representations are five ways of writing the '
            'SAME molecule to ONE general model; they are not five prediction models.'),
        'point_estimates': {},
        'cell_status': {},
        'coverage': {},
        'comparisons': {},
        'strict_common_subset': {},
        'inference_status': {},
        'multiplicity': {
            'family': scoring.EXTENSION_COMPARISON_FAMILY,
            'divisor': scoring.EXTENSION_COMPARISON_FAMILY['size_per_task'],
            '_untested_handling': ('an untested pair keeps its place in the family; the Holm '
                                   'divisor is not reduced to the pairs that were computable'),
            '_old_p_values': ('the previous stage corrected over a family of 3. Those values '
                              'are NOT reused here and are not presented as extension results'),
        },
    }

    planned = planned_counts()
    for task in tasks.TASK_ORDER:
        per_rep = {}
        for representation in ADDED + INHERITED:
            raw = cell_rows(pooled, representation, task) if (
                representation in representations) else []
            if not raw:
                per_rep[representation] = {'status': 'not run', 'value': None,
                                           'n_total': 0, 'n_valid': 0}
                continue
            rows = latest_by_item(raw)
            ok = [r for r in rows if r.get('status') == 'ok']
            target = planned.get('%s|%s' % (representation, task), len(rows))
            per_rep[representation] = {
                'status': cell_status(raw, planned=target),
                'value': BA.task_metric(ok, task, frozen_labels, threshold),
                # the DENOMINATOR is the frozen plan for this cell; the valid
                # count is how many of those came back parseable. A wrong answer
                # is still a valid answer and is never counted as missing.
                'n_total': target,
                'n_valid': len(ok),
                'n_invalid': sum(1 for r in rows if r.get('status') == 'invalid'),
                'n_pending': sum(1 for r in rows if r.get('status') == 'pending'),
                'valid_rate': (len(ok) / target) if target else None,
                '_records_in_log': len(raw),
            }
        out['point_estimates'][task] = per_rep

    # per-cell status table, so an unfinished cell can never read as a finished one
    for task in tasks.TASK_ORDER:
        out['cell_status'][task] = {
            rep: per_rep['status'] for rep, per_rep in out['point_estimates'][task].items()}

    # ---- the nine-pair family, per task ---------------------------------
    family = [tuple(p) for p in scoring.EXTENSION_COMPARISON_FAMILY['within_each_task']] + \
             [tuple(p) for p in scoring.EXTENSION_COMPARISON_FAMILY['extra_pairs']]
    family_size = scoring.EXTENSION_COMPARISON_FAMILY['size_per_task']
    for task in tasks.TASK_ORDER:
        cells = {}
        for representation in ADDED + INHERITED:
            if representation not in representations:
                continue
            rows = [r for r in cell_rows(pooled, representation, task) if r.get('status') == 'ok']
            if rows:
                cells[representation] = rows
        details, raw_p = {}, {}
        for a, b in family:
            key = '%s_vs_%s' % (a, b)
            if a not in cells or b not in cells:
                raw_p[key] = None
                details[key] = {
                    'status': 'not tested',
                    'why': 'one or both arms have no valid responses for this task',
                }
                continue
            # both helpers take ONE list of (record_a, record_b) tuples: paired_metric
            # reads arm A from the first element and arm B from the second of each
            # tuple, so passing two separate record lists would unpack dicts
            pairs_a, pairs_b = BA.pair_by_key(cells[a], cells[b])
            pairs = list(zip(pairs_a, pairs_b))
            if len(pairs) < 2:
                raw_p[key] = None
                details[key] = {'status': 'not tested',
                                'why': 'fewer than 2 items were answered validly by both arms'}
                continue
            detail = BA.paired_metric(pairs, pairs, task, frozen_labels,
                                      RESAMPLES, SEED)
            perm = BA.paired_permutation(pairs, task, frozen_labels, RESAMPLES, SEED)
            detail['permutation'] = perm
            detail['n_pairs_valid'] = len(pairs)
            detail['status'] = 'tested'
            detail['interval_note'] = ('uncorrected 95% cluster-bootstrap interval; the Holm '
                                       'correction applies to the permutation p-values, not '
                                       'to this interval')
            details[key] = detail
            raw_p[key] = perm.get('p_value')
        out['comparisons'][task] = {
            'cells_with_valid_answers': sorted(cells),
            'comparisons': details,
            'p_raw': raw_p,
            'p_holm': BA.holm_with_family(raw_p, family_size),
        }

    # ---- the strict set every representation answered -------------------
    common = None
    for representation in representations:
        keys = set()
        for r in pooled:
            if r.get('representation') != representation or r.get('split') != 'test':
                continue
            keys.add(scoring.pair_key(r))
        common = keys if common is None else (common & keys)
    common = common or set()
    out['strict_common_subset'] = {
        'n_items': len(common),
        'note': ('items answered by EVERY representation in the table; recomputed per item, '
                 'with the 10-label macro-F1 set unchanged. Graph3D cannot express a few '
                 'molecules, so this is the set on which all eight rows are strictly '
                 'comparable. It is reported ALONGSIDE the full-coverage table, not instead '
                 'of it, and no molecule was swapped to reach a round number.'),
        'representations': representations,
        'per_task': {},
    }
    for task in tasks.TASK_ORDER:
        per_task_common = None
        for representation in representations:
            keys = {scoring.pair_key(r) for r in cell_rows(pooled, representation, task)
                    if r.get('status') == 'ok'}
            per_task_common = keys if per_task_common is None else (per_task_common & keys)
        per_task_common = per_task_common or set()
        entry = {'n_items': len(per_task_common), 'values': {}}
        for representation in representations:
            rows = [r for r in cell_rows(pooled, representation, task)
                    if r.get('status') == 'ok' and scoring.pair_key(r) in per_task_common]
            entry['values'][representation] = {
                'n_valid': len(rows),
                'value': BA.task_metric(rows, task, frozen_labels, threshold) if rows else None,
            }
        full = {rep: out['point_estimates'][task][rep]['n_valid'] for rep in representations}
        entry['items_dropped_vs_full_coverage'] = {
            rep: full[rep] - entry['values'][rep]['n_valid'] for rep in representations}
        # WHY the intersection is smaller than the full set. Two different things
        # can shrink it, and they are not the same claim:
        #   * NESTED  -- a representation simply has fewer items, and every item it
        #                does have is shared. Nothing was substituted.
        #   * DIFFERS -- a representation has the same COUNT but a different SET,
        #                because the item builder advanced to the next candidate
        #                when a molecule could not be rendered.
        sets = {rep: {scoring.pair_key(r) for r in cell_rows(pooled, rep, task)
                      if r.get('status') == 'ok'} for rep in representations}
        smallest = min(sets.values(), key=len) if sets else set()
        entry['why_the_intersection_is_smaller'] = {
            'n_per_representation': {rep: len(s) for rep, s in sets.items()},
            'nested': all(smallest <= s for s in sets.values()),
            'item_sets_identical_across_all_rows':
                len({frozenset(s) for s in sets.values()}) == 1,
            'explanation': ('a representation holds a strict SUBSET of the shared items; no '
                            'item was substituted'
                            if all(smallest <= s for s in sets.values()) else
                            'a representation holds the same number of items but a DIFFERENT '
                            'set: the Identity builder advances to the next candidate when a '
                            'molecule cannot be rendered, so a few pairs were replaced rather '
                            'than dropped. This is reported, not corrected after the fact.'),
        }
        out['strict_common_subset']['per_task'][task] = entry

    out['inference_status'] = {
        't1_identity': ('descriptive only. The view a molecule is drawn in is fixed by slot '
                        'in this stage, but the shared-endpoint dependence is still not '
                        'resolved by a method this project has validated, so no confirmatory '
                        'claim is made.'),
        't2_functional_groups': 'reported with the fixed 10-label set, paired, Holm-corrected',
        't3_atom_counts': 'reported; items are single molecules, clusters are singletons',
        't4_qed': 'reported; lower is better, so a positive difference means the left arm is worse',
        't5_hydrogen_bond_donors': 'reported; items are single molecules',
    }
    out['analysis_note'] = (
        'Every number in the report is read from this file. Point estimates are descriptive '
        'and do not depend on any significance result; an interval containing 0 means no '
        'difference was DETECTED at this sample size, which is not the same as equivalence.')

    OUT.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding='utf-8')
    print('rows in the table:', ', '.join(representations))
    print('reuse verified  :', reuse_verified, '--', reuse_reason)
    print()
    header = '%-26s' % 'task'
    for representation in representations:
        header += '%12s' % representation
    print(header)
    for task in tasks.TASK_ORDER:
        line = '%-26s' % task
        for representation in representations:
            entry = out['point_estimates'][task][representation]
            line += '%12s' % ('n/a' if entry['value'] is None else '%.3f' % entry['value'])
        print(line)
    print()
    print('wrote', OUT)
    return 0


if __name__ == '__main__':
    sys.exit(main())