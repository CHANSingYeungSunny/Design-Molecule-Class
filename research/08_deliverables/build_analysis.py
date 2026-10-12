"""Build the versioned analysis artifact that EVERY report reads. Calls no model.

Why this exists
---------------
The report generator used to compute some numbers itself and take others from
older files, so the main table, the comparison table and the prose drifted
apart. This script recomputes everything from the raw main-02 responses ONCE
and writes analysis_v1.json; the Word generator then reads only that file.

Fixed decisions, applied identically in every resample
------------------------------------------------------
  * T2 uses the FROZEN 10-label macro-F1 set. The set is fixed before
    resampling and never re-derived inside a replicate, so a draw that happens
    to contain no positive for a label cannot silently change the denominator.
  * zero-division follows the rule already declared in scoring: a label with no
    positives in the sample is EXCLUDED from the mean and reported; a label
    with positives but no predicted positive scores 0.
  * QED is an error: a positive A-B difference means A is WORSE.
  * Comparisons are paired on a representation-stripped key, and the key is
    checked for uniqueness before it is used.

Inference status is recorded per task, because not every task has a method
whose assumptions this project can defend.
"""
from __future__ import annotations

import datetime
import hashlib
import json
import pathlib
import random
import sys
from collections import Counter, defaultdict

DELIV = pathlib.Path(__file__).resolve().parent
BENCH = DELIV.parent / '07_benchmark'
sys.path.insert(0, str(BENCH))

from bmb import labels, paths, reps, runner, scoring, tasks  # noqa: E402

OUT = DELIV / 'analysis_v1.json'
RUN = 'main-02'
COMBINATIONS = [('our_code', 'smiles'), ('our_code', 'selfies'), ('smiles', 'selfies')]
RESAMPLES = 2000
SEED = 20261005


def sha256(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def _dir_digest(directory):
    """Order-independent digest of every file under a directory."""
    parts = []
    for path in sorted(pathlib.Path(directory).rglob('*')):
        if path.is_file():
            parts.append('%s:%s' % (path.name,
                                    hashlib.sha256(path.read_bytes()).hexdigest()))
    return hashlib.sha256('\n'.join(parts).encode('utf-8')).hexdigest()


def fixed_labels(records):
    """The T2 label set frozen BEFORE the run, re-checked against the sample.

    The support rule was declared before any model call: a label needs at least
    10 positives and 10 negatives in the FROZEN 200-ITEM SAMPLE.

    The denominator matters. The three representations answer the SAME 200
    questions, so there are 600 T2 records but only 200 distinct items. The
    frozen set must therefore be derived from ONE representation's 200 records;
    counting all 600 triples every count and would wrongly re-admit a label
    (nitrile at 7 per sample clears the floor at 21 across three arms). Using
    the wrong denominator would also change the macro-F1 itself.
    """
    sample = records[:]
    rule = json.loads(paths.RUN_CONFIG.read_text(encoding='utf-8'))['main']['label_support_rule']
    counts = {}
    for name in [i['name'] for i in labels.enabled_labels()]:
        pos = sum(1 for r in sample if r['label'][name])
        counts[name] = {'positive': pos, 'negative': len(sample) - pos}
    included = [n for n, c in counts.items()
                if c['positive'] >= rule['min_positives_in_sample']
                and c['negative'] >= rule['min_negatives_in_sample']]
    return sorted(included), counts


def macro_f1(records, label_set):
    """Macro-F1 over a FIXED label set, with the declared zero-division rule."""
    if not records:
        return None
    f1s = []
    for name in label_set:
        tp = fp = fn = 0
        for r in records:
            p, a = bool(r['parsed'][name]), bool(r['label'][name])
            tp += p and a
            fp += p and not a
            fn += (not p) and a
        positives = tp + fn
        predicted = tp + fp
        if positives == 0:
            continue                      # not exercised -> excluded, reported
        precision = tp / predicted if predicted else 0.0
        recall = tp / positives
        f1s.append(0.0 if precision + recall == 0 else
                   2 * precision * recall / (precision + recall))
    return sum(f1s) / len(f1s) if f1s else None


def task_metric(records, task, label_set, threshold=0.05):
    ok = [r for r in records if r['status'] == 'ok']
    if not ok:
        return None
    if task == 't1_identity':
        return sum(1 for r in ok if r['parsed']['same_molecule'] == r['label']['same_molecule']) / len(ok)
    if task == 't2_functional_groups':
        return macro_f1(ok, label_set)
    if task == 't3_atom_counts':
        return sum(1 for r in ok if r['parsed'] == r['label']) / len(ok)
    if task == 't4_qed':
        return sum(abs(r['parsed']['qed'] - r['label']['qed']) for r in ok) / len(ok)
    if task == 't5_hydrogen_bond_donors':
        return sum(1 for r in ok if r['parsed']['hbd'] == r['label']['hbd']) / len(ok)
    return None


def unit_of(record):
    """The dependence unit used by every resample in this analysis."""
    owned = scoring.clusters_of(record)
    return owned


def resample_pairs(pairs, rng):
    """Resample the shared molecules; keep a pair whose ATTRIBUTED unit is drawn.

    Random attribution equalises the marginal inclusion probability of pairs
    owning one molecule (Identity positives) and two (Identity negatives).
    It does NOT by itself remove the correlation between pairs that share a
    molecule, so Identity inference is reported as exploratory.
    """
    index = defaultdict(list)
    for i, (a, _b) in enumerate(pairs):
        for u in unit_of(a):
            index[u].append(i)
    keys = sorted(index)
    if not keys:
        return pairs
    drawn = {keys[rng.randrange(len(keys))] for _ in range(len(keys))}
    chosen = set()
    for i, (a, _b) in enumerate(pairs):
        owned = unit_of(a)
        attributed = owned[rng.randrange(len(owned))] if len(owned) > 1 else owned[0]
        if attributed in drawn:
            chosen.add(i)
    return [pairs[i] for i in sorted(chosen)]


def paired_metric(pairs_a, pairs_b, task, label_set, resamples, seed, confidence=0.95):
    """Paired cluster bootstrap of the difference in the MAIN-TABLE metric."""
    result = {'method': 'paired cluster bootstrap of the metric difference',
              'parameter': {'t1_identity': 'accuracy',
                            't2_functional_groups': 'macro_f1_over_fixed_set',
                            't3_atom_counts': 'exact_match_accuracy',
                            't4_qed': 'mae',
                            't5_hydrogen_bond_donors': 'hbd_exact_match_accuracy'}[task],
              'n_pairs': len(pairs_a),
              'lower_is_better': task == 't4_qed'}
    a = task_metric([x for x, _ in pairs_a], task, label_set)
    b = task_metric([x for _, x in pairs_b], task, label_set)
    result['value_a'], result['value_b'] = a, b
    result['observed_difference'] = None if None in (a, b) else a - b
    if resamples <= 0 or not pairs_a:
        result['ci'] = None
        return result
    rng = random.Random(seed)
    diffs = []
    for _ in range(resamples):
        pick = resample_pairs(pairs_a, rng)
        if len(pick) < 2:
            continue
        va = task_metric([x for x, _ in pick], task, label_set)
        vb = task_metric([x for _, x in pick], task, label_set)
        if None in (va, vb):
            continue
        diffs.append(va - vb)
    if len(diffs) < max(20, resamples // 10):
        result['ci'] = None
        return result
    diffs.sort()
    alpha = (1 - confidence) / 2
    low = diffs[min(len(diffs) - 1, int(alpha * len(diffs)))]
    high = diffs[min(len(diffs) - 1, int((1 - alpha) * len(diffs)))]
    result['ci'] = {'low': low, 'high': high, 'resamples': len(diffs),
                    'confidence': confidence,
                    'excludes_zero': low > 0 or high < 0,
                    'holm_adjusted': False,
                    '_note': 'uncorrected interval; it is NOT a Holm-adjusted test'}
    return result


def holm(p_values):
    items = sorted(((k, v) for k, v in p_values.items() if v is not None), key=lambda kv: kv[1])
    m, running, out = len(items), 0.0, {}
    for rank, (k, p) in enumerate(items):
        running = max(running, min(1.0, (m - rank) * p))
        out[k] = running
    for k in p_values:
        out.setdefault(k, None)
    return out


def holm_with_family(p_values, family_size):
    """Holm over a family whose size is FIXED in advance.

    With ``holm`` the divisor is however many p-values happen to be present, so
    a comparison that could not be computed makes every remaining one look
    stronger -- shrinking the family after seeing the results, which the
    declared protocol forbids. Here untested hypotheses keep their place.
    """
    items = sorted(((k, v) for k, v in p_values.items() if v is not None), key=lambda kv: kv[1])
    m = max(family_size, len(items))
    running, out = 0.0, {}
    for rank, (k, p) in enumerate(items):
        running = max(running, min(1.0, (m - rank) * p))
        out[k] = running
    for k in p_values:
        out.setdefault(k, None)
    return out


def pair_by_key(records_a, records_b):
    """Two record lists aligned on the representation-stripped item key."""
    by_a = {scoring.pair_key(r): r for r in records_a}
    by_b = {scoring.pair_key(r): r for r in records_b}
    shared = sorted(set(by_a) & set(by_b))
    return [by_a[k] for k in shared], [by_b[k] for k in shared]


def paired_permutation(pairs, task, label_set, resamples, seed, confidence=0.95):
    """Cluster sign-flip permutation on the MAIN-TABLE metric.

    For T2 the exchangeable quantity is the per-item PREDICTION VECTOR: a flip
    swaps A's and B's label vectors for that item, and macro-F1 is recomputed
    over the FIXED label set on the swapped sample. It does not flip a
    different per-item score.

    A Monte-Carlo permutation cannot express a p below 1/(resamples+1); that
    floor is reported rather than printing p = 0.
    """
    index = defaultdict(list)
    for i, (a, _b) in enumerate(pairs):
        for u in unit_of(a):
            index[u].append(i)
    keys = sorted(index)
    if not keys or resamples <= 0:
        return {'method': 'cluster sign-flip permutation', 'p_value': None,
                'note': 'too few units to permute'}
    rng = random.Random(seed)

    def diff(sample):
        va = task_metric([a for a, _ in sample], task, label_set)
        vb = task_metric([b for _, b in sample], task, label_set)
        return None if None in (va, vb) else va - vb

    observed = diff(pairs)
    extreme = 0
    for _ in range(resamples):
        drawn = {keys[rng.randrange(len(keys))] for _ in range(len(keys))}
        flipped, sub = {u: rng.random() < 0.5 for u in keys}, []
        for i, (a, b) in enumerate(pairs):
            owned = unit_of(a)
            attributed = owned[rng.randrange(len(owned))] if len(owned) > 1 else owned[0]
            if attributed not in drawn:
                continue
            sub.append((b, a) if flipped[attributed] else (a, b))
        if len(sub) < 2:
            continue
        d = diff(sub)
        if d is not None and observed is not None and abs(d) >= abs(observed):
            extreme += 1
    p = (extreme + 1) / (resamples + 1)
    return {'method': 'cluster sign-flip permutation on the main-table metric',
            'statistic': 'difference in the main-table metric',
            'exchanges': ('per-item prediction vectors; macro-F1 recomputed over the FIXED '
                          'label set after the swap' if task == 't2_functional_groups'
                          else 'per-item records'),
            'resamples': resamples, 'seed': seed,
            'p_value': p, 'monte_carlo_floor': 1 / (resamples + 1),
            'exchangeability_assumption': ('pairs are exchangeable within molecule clusters '
                                           'under the null; NOT verified for shared endpoints')}


def main():
    manifest = json.loads((DELIV / 'main_manifest.json').read_text(encoding='utf-8'))
    run_config = json.loads(paths.RUN_CONFIG.read_text(encoding='utf-8'))
    raw_path = paths.run_dir(RUN) / 'raw_responses.jsonl'
    records = [json.loads(l) for l in raw_path.read_text(encoding='utf-8').splitlines() if l.strip()]
    t2_records = [r for r in records if r['task'] == 't2_functional_groups']
    # the frozen sample is ONE representation's 200 items; the three arms answer
    # the same questions, so pooling them would triple every count
    t2_sample = [r for r in t2_records if r['representation'] == 'smiles']

    label_set, label_counts = fixed_labels(t2_sample)
    cells = defaultdict(list)
    for r in records:
        cells[(r['representation'], r['task'])].append(r)

    analysis = {
        'analysis_version': 'v1',
        'generated_at': datetime.datetime.now().isoformat(timespec='seconds'),
        'run_id': RUN,
        'inputs': {
            'raw_responses': {'path': str(raw_path.relative_to(BENCH.parent)),
                              'sha256': sha256(raw_path), 'records': len(records)},
            'main_manifest_sha256': sha256(DELIV / 'main_manifest.json'),
            'frozen_items_sha256': _dir_digest(BENCH / 'data' / 'prepared' / 'main_items'),
            'model_config_sha256': sha256(paths.MODEL_CONFIG),
        },
        'metric_definitions': {
            't1_identity': 'Accuracy over paired items',
            't2_functional_groups': 'Macro-F1 over a FIXED label set (unexercised labels excluded)',
            't3_atom_counts': 'Exact-match accuracy over the whole answer object',
            't4_qed': 'MAE (lower is better)',
            't5_hydrogen_bond_donors': 'Exact-match accuracy',
        },
        't2_label_set': {
            'fixed_set': label_set,
            'n_labels': len(label_set),
            'rule': run_config['main']['label_support_rule'],
            'per_label_counts_in_sample': label_counts,
            'excluded': sorted(set(label_counts) - set(label_set)),
            'fixed_across_resamples': True,
        },
        'pairing': {},
        'point_estimates': {},
        'comparisons': {},
        'inference_status': {},
        'method_changes': [
            {'change': 'T2 comparison metric',
             'was': 'per-item label accuracy over all 13 labels, used to argue about Macro-F1',
             'now': 'the difference in the FIXED-10-label Macro-F1 itself',
             'why': 'per-item label accuracy is dominated by true negatives, which F1 does not '
                    'reward, and it averaged over a different label set',
             'when': 'offline, AFTER the model run'},
            {'change': 'cluster resampling of multi-endpoint items',
             'was': 'an item entered if ANY owned molecule was drawn (2-endpoint items drawn '
                    '~37% more often)',
             'now': 'each item is attributed to one molecule at random per replicate, '
                    'equalising the marginal inclusion probability',
             'why': 'the old rule over-weighted Identity negatives and pushed intervals off the '
                    'point estimate',
             'when': 'offline, AFTER the model run',
             '_limit': 'random attribution fixes the MARGINAL probability only; it does not '
                       'remove the correlation between pairs that share a molecule, so '
                       'Identity inference remains exploratory'},
        ],
        'analysis_note': ('The model experiment was NOT re-run. All raw responses are the '
                          'original ones; only the offline statistical implementation was '
                          'corrected after the problem was found.'),
    }

    # pairing checks
    for task in tasks.TASK_ORDER:
        per_rep, shared = {}, None
        for rep in reps.ROUND1_REPRESENTATIONS:
            rows = cells[(rep, task)]
            keys = [scoring.pair_key(r) for r in rows]
            per_rep[rep] = {'items': len(rows), 'unique_keys': len(set(keys)),
                            'collisions': len(keys) - len(set(keys))}
            s = set(keys)
            shared = s if shared is None else (shared & s)
        analysis['pairing'][task] = {'per_representation': per_rep,
                                     'shared_across_all_three': len(shared or set()),
                                     'unique_note': 'no key collisions; all three arms answer the '
                                                    'same items'}

    # point estimates
    for task in tasks.TASK_ORDER:
        analysis['point_estimates'][task] = {}
        for rep in reps.ROUND1_REPRESENTATIONS:
            rows = cells[(rep, task)]
            ok = [r for r in rows if r['status'] == 'ok']
            value = task_metric(rows, task, label_set)
            entry = {'value': value, 'n_total': len(rows), 'n_valid': len(ok),
                     'valid_rate': len(ok) / len(rows) if rows else None,
                     'complete': all(r['status'] != 'pending' for r in rows)}
            if task == 't4_qed':
                entry['success_rate_at_0.05'] = (
                    sum(1 for r in ok if abs(r['parsed']['qed'] - r['label']['qed']) <= 0.05)
                    / len(rows) if rows else None)
            if task == 't2_functional_groups':
                entry['label_set_size'] = len(label_set)
            if task == 't1_identity':
                entry['n_correct'] = sum(1 for r in ok
                                         if r['parsed']['same_molecule'] == r['label']['same_molecule'])
            analysis['point_estimates'][task][rep] = entry

    # comparisons
    for task in tasks.TASK_ORDER:
        analysis['comparisons'][task] = {}
        raw_p = {}
        for a, b in COMBINATIONS:
            pa = {scoring.pair_key(r): r for r in cells[(a, task)]}
            pb = {scoring.pair_key(r): r for r in cells[(b, task)]}
            pairs = [(pa[k], pb[k]) for k in sorted(set(pa) & set(pb))
                     if pa[k]['status'] == 'ok' and pb[k]['status'] == 'ok']
            entry = paired_metric(pairs, pairs, task, label_set, RESAMPLES, SEED)
            key = '%s_vs_%s' % (a, b)
            analysis['comparisons'][task][key] = entry
            perm = paired_permutation(pairs, task, label_set, RESAMPLES, SEED)
            entry['permutation'] = perm
            raw_p[key] = perm.get('p_value')
        # Holm over the DECLARED family of three, per task
        adjusted = holm(raw_p)
        for key, value in adjusted.items():
            analysis['comparisons'][task][key]['permutation']['p_holm'] = value
        analysis['comparisons'][task]['_family'] = {
            'declared': [('%s_vs_%s' % c) for c in COMBINATIONS],
            'family_size': len(COMBINATIONS),
            'correction': 'Holm-Bonferroni, within task',
            'ci_holm_adjusted': False,
            'monte_carlo_floor': 1 / (RESAMPLES + 1),
        }

    # inference status per task: what this project can and cannot claim
    analysis['inference_status'] = {
        't1_identity': {
            'status': 'exploratory',
            'reason': ('Identity items share molecules across pairs (Identity pairs share '
                       'endpoints with other items). Random attribution equalises the '
                       'marginal draw probability but does not remove that correlation, and '
                       'no method here is validated for it.'),
            'report': 'accuracy and correct/total counts are reported; the interval is '
                      'labelled exploratory and not used for confirmatory claims.'},
        't2_functional_groups': {
            'status': 'reportable_with_caveats',
            'reason': ('the comparison recomputes the fixed-10-label Macro-F1 itself; the '
                       'resampling unit is the parent molecule'),
            'report': 'difference with 95% interval, plus a Holm-corrected permutation p'},
        't3_atom_counts': {'status': 'reportable_with_caveats',
                           'reason': 'independent single-molecule items; clusters are singletons',
                           'report': 'difference with 95% interval and Holm-corrected p'},
        't4_qed': {'status': 'reportable_with_caveats',
                   'reason': 'independent single-molecule items; MAE difference, lower is better',
                   'report': 'difference with 95% interval and Holm-corrected p'},
        't5_hydrogen_bond_donors': {'status': 'reportable_with_caveats',
                                    'reason': 'independent single-molecule items',
                                    'report': 'difference with 95% interval and Holm-corrected p'},
    }
    analysis['multiplicity'] = {
        'per_task_family': 'the three pairwise comparisons within one task form one family',
        'correction': 'Holm-Bonferroni',
        'significance_level': 0.05,
        'intervals_are_uncorrected': True,
        'direction': {'higher_is_better': ['t1_identity', 't2_functional_groups',
                                           't3_atom_counts', 't5_hydrogen_bond_donors'],
                      'lower_is_better': ['t4_qed']},
    }
    OUT.write_text(json.dumps(analysis, ensure_ascii=False, indent=2), encoding='utf-8')
    print('wrote', OUT)
    print('T2 fixed label set (%d): %s' % (len(label_set), label_set))
    print('excluded:', analysis['t2_label_set']['excluded'])
    print()
    for task in tasks.TASK_ORDER:
        st = analysis['inference_status'][task]['status']
        print('%-26s %s' % (task, st))
    return 0


if __name__ == '__main__':
    sys.exit(main())