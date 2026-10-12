"""Evaluator: frozen scoring rules, applied after the responses exist.

Rules fixed before any test item is run (protocol section 6):

  status       every scheduled item ends as 'ok' or 'invalid'. A transport
               failure that survives the declared retries is 'pending' and the
               cell is reported as incomplete, never silently dropped.
  invalid      an unparseable or out-of-range answer is 'invalid'. It counts as
               a failure in the all-sample metrics. It is never coerced into a
               correct value, never clipped into range, never dropped.
  missing      a missing answer is NOT treated as zero error. It is a failure in
               the all-sample success rate, and excluded (with n stated) from
               the conditional error metrics.
  t1           Accuracy, Balanced Accuracy, confusion matrix, valid rate.
  t2           Macro-F1 on the valid subset AND a penalised Macro-F1 over all
               samples where an invalid answer is scored as all-false (a fixed
               penalty, declared here, not a guess at the model's intent).
  t3           Exact-match accuracy over the whole answer object, plus per-element.
  t4           MAE and Spearman over valid numbers with n_valid/n_total stated;
               all-sample success rate at the pre-declared 0.05 threshold.
  t5           Exact-match accuracy.

Confidence intervals are 95% bootstrap, resampled by parent molecule (cluster),
because several items derive from one molecule and are not independent. F1 and
the correlation are recomputed inside each resample rather than averaged from
per-item scores.
"""
from __future__ import annotations

import math
import random


def _mean(values):
    return sum(values) / len(values) if values else None


def spearman(pairs):
    """Spearman rank correlation; returns None when it is undefined."""
    n = len(pairs)
    if n < 3:
        return None
    xs = [p[0] for p in pairs]
    ys = [p[1] for p in pairs]

    def ranks(values):
        order = sorted(range(len(values)), key=lambda i: values[i])
        out = [0.0] * len(values)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
                j += 1
            average = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                out[order[k]] = average
            i = j + 1
        return out

    rx, ry = ranks(xs), ranks(ys)
    mx, my = _mean(rx), _mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    if dx == 0 or dy == 0:
        return None
    return num / (dx * dy)


def f1(tp, fp, fn):
    if tp == 0 and (fp or fn):
        return 0.0
    if tp == 0 and not fp and not fn:
        return None
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    if precision + recall == 0:
        return 0.0
    return 2 * precision * recall / (precision + recall)


# --------------------------------------------------------------------------
# per-item scoring
# --------------------------------------------------------------------------

def score_item(record):
    """Return (correct_or_error, per-item detail) for one response record."""
    task = record['task']
    label = record['label']
    parsed = record.get('parsed')
    if record.get('status') != 'ok' or parsed is None:
        return None
    if task == 't1_identity':
        return bool(parsed['same_molecule'] == label['same_molecule'])
    if task == 't3_atom_counts':
        return bool(parsed == label)
    if task == 't4_qed':
        return abs(parsed['qed'] - label['qed'])
    if task == 't5_hydrogen_bond_donors':
        return bool(parsed['hbd'] == label['hbd'])
    return None


# --------------------------------------------------------------------------
# metric blocks
# --------------------------------------------------------------------------

def _identity_metrics(records):
    tp = fp = tn = fn = 0
    for record in records:
        if record['status'] != 'ok':
            continue
        predicted = record['parsed']['same_molecule']
        actual = record['label']['same_molecule']
        if predicted and actual:
            tp += 1
        elif predicted and not actual:
            fp += 1
        elif not predicted and actual:
            fn += 1
        else:
            tn += 1
    total = tp + fp + tn + fn
    accuracy = (tp + tn) / total if total else None
    recall_pos = tp / (tp + fn) if (tp + fn) else None
    recall_neg = tn / (tn + fp) if (tn + fp) else None
    balanced = (_mean([r for r in (recall_pos, recall_neg) if r is not None])
                if (recall_pos is not None or recall_neg is not None) else None)
    return {
        'accuracy': accuracy,
        'balanced_accuracy': balanced,
        'confusion': {'tp': tp, 'fp': fp, 'tn': tn, 'fn': fn},
        'n_scored': total,
    }


def _functional_group_metrics(records, label_names):
    stats = {name: {'tp': 0, 'fp': 0, 'fn': 0, 'tn': 0} for name in label_names}
    for record in records:
        if record['status'] != 'ok':
            continue
        predicted = record['parsed']
        actual = record['label']
        for name in label_names:
            p, a = bool(predicted[name]), bool(actual[name])
            key = {(True, True): 'tp', (True, False): 'fp',
                   (False, True): 'fn', (False, False): 'tn'}[(p, a)]
            stats[name][key] += 1

    per_label = {}
    for name, s in stats.items():
        tp, fp, fn, tn = s['tp'], s['fp'], s['fn'], s['tn']
        positives = tp + fn
        predicted_positive = tp + fp
        # Declared convention, stated here so nothing depends on an implicit
        # library default:
        #   recall  undefined when the sample holds no positive for that label
        #   precision undefined when the model predicted no positive
        #   a label with no positives at all is NOT EXERCISED -> F1 None and it
        #     is excluded from macro-F1 rather than scored as a perfect 1.0
        #   recall defined but precision undefined (positives exist, none
        #     predicted) is a real miss -> precision taken as 0 -> F1 0.0
        recall = (tp / positives) if positives else None
        precision = (tp / predicted_positive) if predicted_positive else None
        if recall is None:
            f1_value = None
        elif precision is None:
            f1_value = 0.0
        elif precision + recall == 0:
            f1_value = 0.0
        else:
            f1_value = 2 * precision * recall / (precision + recall)
        per_label[name] = {
            'f1': f1_value,
            'precision': precision,
            'recall': recall,
            'tp': tp, 'fp': fp, 'fn': fn, 'tn': tn,
            'positive': positives,
            'negative': tn + fp,
            'predicted_positive': predicted_positive,
            'exercised': recall is not None,
        }
    valid_f1 = [v['f1'] for v in per_label.values() if v['f1'] is not None]
    macro_valid = _mean(valid_f1)
    macro_rule = {
        'macro_f1_valid_subset': ('unweighted mean of per-label F1 over EXERCISED labels only '
                                  '(those with at least one positive in the scored sample); '
                                  'unexercised labels are excluded and listed, never counted as 1.0'),
        'macro_f1_all_samples_penalised': ('same labels; an invalid or missing answer is scored '
                                           'as predicting False for every label'),
        'zero_division': 'declared explicitly in this function, not inherited from a library',
    }

    # penalised variant over ALL scheduled samples: an invalid answer is scored
    # as predicting False for every group (fixed rule, declared in the header)
    penal = {name: {'tp': 0, 'fp': 0, 'fn': 0} for name in label_names}
    for record in records:
        if record['status'] == 'ok':
            source = {name: bool(record['parsed'][name]) for name in label_names}
        else:
            source = {name: False for name in label_names}
        for name in label_names:
            actual = bool(record['label'][name])
            if source[name] and actual:
                penal[name]['tp'] += 1
            elif source[name] and not actual:
                penal[name]['fp'] += 1
            elif not source[name] and actual:
                penal[name]['fn'] += 1
    penal_f1 = [f1(v['tp'], v['fp'], v['fn']) for v in penal.values()]
    macro_penalised = _mean([v for v in penal_f1 if v is not None])

    return {
        'macro_f1_valid_subset': macro_valid,
        'macro_f1_all_samples_penalised': macro_penalised,
        'per_label': per_label,
        'macro_f1_rule': macro_rule,
        'n_labels_exercised': len(valid_f1),
        'labels_not_exercised': sorted(n for n, v in per_label.items() if not v['exercised']),
    }


def _atom_count_metrics(records):
    exact = 0
    scored = 0
    per_element = {}
    for record in records:
        if record['status'] != 'ok':
            continue
        scored += 1
        predicted, actual = record['parsed'], record['label']
        exact += int(predicted == actual)
        for element in actual:
            bucket = per_element.setdefault(element, {'correct': 0, 'n': 0})
            bucket['n'] += 1
            bucket['correct'] += int(predicted[element] == actual[element])
    for element, bucket in per_element.items():
        bucket['accuracy'] = bucket['correct'] / bucket['n'] if bucket['n'] else None
    return {
        'exact_match_accuracy': exact / scored if scored else None,
        'per_element_accuracy': per_element,
        'n_scored': scored,
    }


def _qed_metrics(records, threshold):
    pairs = []
    errors = []
    successes = 0
    for record in records:
        if record['status'] != 'ok':
            continue
        error = abs(record['parsed']['qed'] - record['label']['qed'])
        pairs.append((record['label']['qed'], record['parsed']['qed']))
        errors.append(error)
        successes += int(error <= threshold)
    total = len(records)
    return {
        'mae': _mean(errors),
        'spearman': spearman(pairs),
        'n_valid': len(errors),
        'n_total': total,
        'success_rate_all_samples': (successes / total) if total else None,
        'success_threshold': threshold,
    }


def _descriptor_metrics(records):
    exact = 0
    scored = 0
    for record in records:
        if record['status'] != 'ok':
            continue
        scored += 1
        exact += int(record['parsed']['hbd'] == record['label']['hbd'])
    return {
        'hbd_exact_match_accuracy': exact / scored if scored else None,
        'n_scored': scored,
    }


# --------------------------------------------------------------------------
# metric selection + bootstrap
# --------------------------------------------------------------------------

PRIMARY_METRIC = {
    't1_identity': ('accuracy', True),
    't2_functional_groups': ('macro_f1_valid_subset', True),
    't3_atom_counts': ('exact_match_accuracy', True),
    't4_qed': ('mae', False),
    't5_hydrogen_bond_donors': ('hbd_exact_match_accuracy', True),
}


def _primary_value(task, records, label_names, threshold):
    if task == 't1_identity':
        return _identity_metrics(records)['accuracy']
    if task == 't2_functional_groups':
        return _functional_group_metrics(records, label_names)['macro_f1_valid_subset']
    if task == 't3_atom_counts':
        return _atom_count_metrics(records)['exact_match_accuracy']
    if task == 't4_qed':
        return _qed_metrics(records, threshold)['mae']
    if task == 't5_hydrogen_bond_donors':
        return _descriptor_metrics(records)['hbd_exact_match_accuracy']
    return None


# --------------------------------------------------------------------------
# interval helpers
# --------------------------------------------------------------------------
#
# Two families of interval are produced, and they answer different questions.
#
#   * Wilson / Clopper-Pearson treat the cell as k successes out of n
#     independent Bernoulli trials. They are reported ONLY for the tasks whose
#     primary metric is a proportion (t1, t3, t5), and ONLY as a reference for
#     what an independent-sampling interval would look like. They are NOT valid
#     when items share a parent molecule or a pair endpoint.
#
#   * The cluster bootstrap resamples the actual dependence unit. It recomputes
#     the statistic inside every resample, so Macro-F1 and the correlation are
#     never averaged from per-item scores.
#
# Neither family can rescue a degenerate cell. When every scored item is
# correct (or every one wrong) at small n, the percentile bootstrap collapses
# to a point interval. That is a property of the data, not a bug, and the
# result records it as such instead of dressing the number up.

PROPORTION_METRIC = {
    't1_identity': 'accuracy',
    't3_atom_counts': 'exact_match_accuracy',
    't5_hydrogen_bond_donors': 'hbd_exact_match_accuracy',
}


def wilson_interval(k, n, confidence=0.95):
    """Wilson score interval for a binomial proportion. Independent samples."""
    if n <= 0:
        return None
    z = _z_for(confidence)
    phat = k / n
    denominator = 1.0 + z * z / n
    centre = (phat + z * z / (2 * n)) / denominator
    half = (z * math.sqrt(phat * (1 - phat) / n + z * z / (4 * n * n))) / denominator
    return {'low': max(0.0, centre - half), 'high': min(1.0, centre + half),
            'method': 'wilson', 'assumes': 'independent Bernoulli items'}


def _z_for(confidence):
    table = {0.90: 1.6448536269514722, 0.95: 1.959963984540054, 0.99: 2.5758293035489004}
    if confidence in table:
        return table[confidence]
    # bisection on the standard normal CDF via erf
    target = 1.0 - (1.0 - confidence) / 2.0
    low, high = 0.0, 10.0
    for _ in range(80):
        mid = (low + high) / 2.0
        if 0.5 * (1.0 + math.erf(mid / math.sqrt(2.0))) < target:
            low = mid
        else:
            high = mid
    return (low + high) / 2.0


def _betacf(a, b, x, max_iter=300, eps=3e-16, tiny=1e-300):
    """Continued fraction for the incomplete beta function (Lentz's method)."""
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < tiny:
        d = tiny
    d = 1.0 / d
    h = d
    for m in range(1, max_iter + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def _betainc(a, b, x):
    """Regularised incomplete beta I_x(a, b)."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lbeta = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
    front = math.exp(lbeta + a * math.log(x) + b * math.log(1.0 - x))
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    return 1.0 - math.exp(lbeta + b * math.log(1.0 - x) + a * math.log(x)) * _betacf(b, a, 1.0 - x) / b


def _beta_quantile(p, a, b, iterations=120):
    """Inverse of _betainc by bisection. p in (0, 1)."""
    low, high = 0.0, 1.0
    for _ in range(iterations):
        mid = (low + high) / 2.0
        if _betainc(a, b, mid) < p:
            low = mid
        else:
            high = mid
    return (low + high) / 2.0


def clopper_pearson_interval(k, n, confidence=0.95):
    """Exact (Clopper-Pearson) interval for a binomial proportion."""
    if n <= 0:
        return None
    alpha = 1.0 - confidence
    low = 0.0 if k == 0 else _beta_quantile(alpha / 2.0, k, n - k + 1)
    high = 1.0 if k == n else _beta_quantile(1.0 - alpha / 2.0, k + 1, n - k)
    return {'low': low, 'high': high, 'method': 'clopper-pearson',
            'assumes': 'independent Bernoulli items'}


def _proportion_counts(task, records):
    """(k, n) for a proportion-type task, or (None, None)."""
    if task not in PROPORTION_METRIC:
        return None, None
    scored = [r for r in records if r['status'] == 'ok']
    if task == 't1_identity':
        correct = sum(1 for r in scored
                      if r['parsed']['same_molecule'] == r['label']['same_molecule'])
    elif task == 't3_atom_counts':
        correct = sum(1 for r in scored if r['parsed'] == r['label'])
    elif task == 't5_hydrogen_bond_donors':
        correct = sum(1 for r in scored if r['parsed']['hbd'] == r['label']['hbd'])
    else:
        return None, None
    return correct, len(scored)


# --------------------------------------------------------------------------
# cluster definition and cluster bootstrap
# --------------------------------------------------------------------------

def clusters_of(record):
    """The dependence units this item belongs to.

    An Identity item is a PAIR. Both endpoint molecules can be shared with other
    items, so a pair is attached to BOTH endpoints; clustering on one endpoint
    alone would ignore the other. Single-molecule tasks carry one cluster.
    """
    clusters = record.get('clusters')
    if clusters:
        return tuple(sorted(set(clusters)))   # a positive pair repeats one molecule
    single = record.get('cluster')
    if single:
        return (single,)
    return (record.get('example_id', 'unknown'),)


def cluster_summary(records):
    """What the resampling unit actually is, for the record."""
    units = set()
    multi = 0
    for record in records:
        owned = clusters_of(record)
        units.update(owned)
        if len(owned) > 1:
            multi += 1
    return {
        'definition': ('each distinct endpoint molecule; a pair item belongs to both of its '
                       'endpoints and enters a resample when at least one endpoint is drawn'),
        'n_clusters': len(units),
        'n_items_sharing_a_cluster': multi,
        'items_are_independent_singletons': len(units) == len(records) and multi == 0,
    }


def _resample(records, clusters, keys, rng):
    """Draw clusters with replacement; keep an item whose ATTRIBUTED cluster was drawn.

    Each item is attributed to ONE of its clusters, chosen uniformly at random
    per resample. An item owning k clusters then has inclusion probability
    1-(1-1/N)^N, the SAME for every item.

    Why this matters: the obvious rule "include the item if ANY owned cluster
    was drawn" gives an item owning 1 cluster p = 1-(1-1/N)^N (~0.633 at N=235)
    but an item owning 2 clusters p = 1-(1-2/N)^N (~0.866). Identity positives
    own one endpoint and negatives own two, so that rule over-weights negatives
    by ~37% and pushes the bootstrap distribution off the point estimate -- it
    is why a 0.540 accuracy could carry a [0.57, 0.65] interval.

    Random attribution still uses BOTH endpoints: across resamples each is
    chosen about half the time, so neither is ignored. It is a declared choice
    among several weighting schemes; it affects only the interval, never the
    point estimate.
    """
    drawn = set()
    for _ in range(len(keys)):
        drawn.add(keys[rng.randrange(len(keys))])
    sample = []
    for record in records:
        owned = clusters_of(record)
        attributed = owned[rng.randrange(len(owned))] if len(owned) > 1 else owned[0]
        if attributed in drawn:
            sample.append(record)
    return sample


def _bootstrap(records, task, label_names, threshold, resamples, seed, confidence=0.95):
    """Percentile cluster bootstrap. The statistic is recomputed every resample."""
    if resamples <= 0 or not records:
        return None
    clusters = {}
    for record in records:
        for unit in clusters_of(record):
            clusters.setdefault(unit, []).append(record)
    keys = sorted(clusters)
    if len(keys) < 2:
        return None
    rng = random.Random(seed)
    values = []
    for _ in range(resamples):
        sample = _resample(records, clusters, keys, rng)
        if not sample:
            continue
        value = _primary_value(task, sample, label_names, threshold)
        if value is not None:
            values.append(value)
    if len(values) < max(20, resamples // 10):
        return None
    values.sort()
    alpha = (1.0 - confidence) / 2.0
    low = values[min(len(values) - 1, int(alpha * len(values)))]
    high = values[min(len(values) - 1, int((1.0 - alpha) * len(values)))]
    return {'low': low, 'high': high, 'resamples': len(values), 'unit': 'cluster',
            'method': 'percentile cluster bootstrap',
            'degenerate': bool(low == high)}


def _interpretation(task, counts, bootstrap, summary):
    """Plain-language note that travels with the interval into the report."""
    n = counts.get('n_scored') if counts else None
    if not summary['items_are_independent_singletons']:
        return ('Items share parent molecules, so the Wilson and Clopper-Pearson '
                'intervals below are shown only as an independence reference; the '
                'cluster bootstrap is the interval that respects the design.')
    if bootstrap and bootstrap.get('degenerate'):
        return ('Every scored item fell on the same side, so the bootstrap collapses '
                'to a point interval. With n=%s this does NOT establish that the '
                'population accuracy is 100%% (or 0%%); it means this sample carries '
                'no information about the spread. Widen the sample before reading '
                'anything into it.' % n)
    return 'Items are independent singletons in this cell; all three intervals are comparable.'


def evaluate(records, task, label_names, threshold, resamples=2000, seed=0):
    """Full metric block for one (representation, task) cell of one run."""
    if task == 't1_identity':
        block = _identity_metrics(records)
    elif task == 't2_functional_groups':
        block = _functional_group_metrics(records, label_names)
    elif task == 't3_atom_counts':
        block = _atom_count_metrics(records)
    elif task == 't4_qed':
        block = _qed_metrics(records, threshold)
    elif task == 't5_hydrogen_bond_donors':
        block = _descriptor_metrics(records)
    else:
        raise KeyError(task)

    total = len(records)
    ok = sum(1 for r in records if r['status'] == 'ok')
    pending = sum(1 for r in records if r['status'] == 'pending')
    block.update({
        'task': task,
        'n_total': total,
        'n_valid': ok,
        'n_invalid': total - ok - pending,
        'n_pending': pending,
        'valid_answer_rate': (ok / total) if total else None,
        'complete': pending == 0,
    })

    summary = cluster_summary(records)
    block['cluster_summary'] = summary

    k, n = _proportion_counts(task, records)
    if k is not None:
        block['proportion'] = {'correct': k, 'scored': n}
        block['ci95_wilson'] = wilson_interval(k, n)
        block['ci95_clopper_pearson'] = clopper_pearson_interval(k, n)
    block['ci95_cluster_bootstrap'] = _bootstrap(records, task, label_names, threshold,
                                                 resamples, seed)
    block['ci95_primary'] = block['ci95_cluster_bootstrap']
    block['interval_interpretation'] = _interpretation(task, block, block['ci95_primary'],
                                                       summary)
    return block


# --------------------------------------------------------------------------
# paired comparison between two representations
# --------------------------------------------------------------------------
#
# Two intervals that overlap (or fail to) say nothing about whether the
# representations differ, because the same item is answered by both arms. The
# test below keeps the pairing: it looks only at items present in BOTH cells and
# resamples by sign-flipping whole clusters, which preserves the dependence
# structure. Multiplicity is handled by a declared family, not by picking the
# one comparison that happens to look good.

def per_item_score(task, record):
    """A higher-is-better per-item score, matched by example_id across arms."""
    if record.get('status') != 'ok':
        return None
    parsed, label = record.get('parsed'), record.get('label')
    if parsed is None:
        return None
    if task == 't1_identity':
        return 1.0 if parsed['same_molecule'] == label['same_molecule'] else 0.0
    if task == 't3_atom_counts':
        return 1.0 if parsed == label else 0.0
    if task == 't5_hydrogen_bond_donors':
        return 1.0 if parsed['hbd'] == label['hbd'] else 0.0
    if task == 't4_qed':
        return -abs(parsed['qed'] - label['qed'])
    if task == 't2_functional_groups':
        names = [k for k in label if not k.startswith('_')]
        if not names:
            return None
        hits = sum(1 for n in names if parsed.get(n) == label[n])
        return hits / len(names)
    return None


def pair_key(record):
    """A representation-INDEPENDENT identity for the item behind a record.

    example_id embeds the representation ("t3-smiles-12621609"), so two arms
    answering the same question have DIFFERENT ids and a naive join finds zero
    pairs. Stripping the representation token recovers the shared identity, so
    the same molecule asked the same question really is one pair.
    """
    rid = record.get('example_id') or ''
    representation = record.get('representation') or ''
    parts = rid.split('-')
    if len(parts) < 2:
        return rid
    prefix = parts[0] + '-'
    rest = rid[len(prefix):]
    if representation and rest.startswith(representation + '-'):
        rest = rest[len(representation) + 1:]
    return prefix + rest


def paired_comparison(records_a, records_b, task, resamples=2000, seed=0, confidence=0.95):
    """Cluster sign-flip permutation test over items answered by both arms."""
    by_id_a = {pair_key(r): r for r in records_a}
    by_id_b = {pair_key(r): r for r in records_b}
    shared = sorted(set(by_id_a) & set(by_id_b))
    pairs = []
    for example_id in shared:
        score_a = per_item_score(task, by_id_a[example_id])
        score_b = per_item_score(task, by_id_b[example_id])
        if score_a is None or score_b is None:
            continue                       # an invalid answer cannot be paired
        pairs.append((example_id, score_a - score_b,
                      clusters_of(by_id_a[example_id])))
    if not pairs:
        return {'n_pairs': 0, 'note': 'no item was answered validly by both arms'}

    observed = _mean([d for _, d, _ in pairs])
    units = sorted({u for _, _, owned in pairs for u in owned})
    result = {
        'n_pairs': len(pairs),
        'n_items_only_in_a': len(set(by_id_a) - set(by_id_b)),
        'n_items_only_in_b': len(set(by_id_b) - set(by_id_a)),
        'n_unpaired_because_invalid': len(shared) - len(pairs),
        'observed_mean_difference': observed,
        'score_definition': ('higher-is-better per-item score; exact-match 0/1 for t1/t3/t5, '
                             'negative absolute error for t4, per-item label accuracy for t2'),
        'cluster_summary': {'n_clusters': len(units)},
    }
    if resamples <= 0 or len(units) < 2:
        result['p_value'] = None
        result['note'] = 'too few clusters to permute'
        return result

    rng = random.Random(seed)
    diffs = [d for _, d, _ in pairs]
    owners = [owned for _, _, owned in pairs]
    extreme = 0
    counter = 0
    for _ in range(resamples):
        flips = {u: (1.0 if rng.random() < 0.5 else -1.0) for u in units}
        total = 0.0
        for diff, owned in zip(diffs, owners):
            sign = 1.0
            for u in owned:
                sign *= flips[u]      # a pair flips only if BOTH endpoints flip
            total += diff * sign
        statistic = total / len(diffs)
        counter += 1
        if abs(statistic) >= abs(observed):
            extreme += 1
    result['p_value'] = (extreme + 1) / (counter + 1) if counter else None
    result['method'] = 'cluster sign-flip permutation, two-sided'
    result['resamples'] = counter
    return result


def paired_metric_bootstrap(records_a, records_b, task, label_names, threshold=0.05,
                            resamples=2000, seed=0, confidence=0.95):
    """Paired bootstrap of the TASK'S OWN metric difference.

    The permutation test above works on a per-item score, which for a
    corpus-level metric like Macro-F1 is a DIFFERENT quantity: per-item label
    accuracy is dominated by the many easy absent labels, so it sits near 0.9
    even when Macro-F1 is near 0.5. Reporting its p-value as evidence about
    Macro-F1 would answer the wrong question.

    This function instead resamples clusters ONCE per replicate, computes the
    metric for arm A and arm B on that same resample, and takes the difference.
    The interval is therefore for exactly the number printed in the main table,
    and the pairing is preserved because both arms see the same draw.

    A difference interval excluding 0 is the evidence for a real difference;
    an interval containing 0 is not evidence of no difference, only of no
    detectable one at this sample size.
    """
    by_a = {pair_key(r): r for r in records_a}
    by_b = {pair_key(r): r for r in records_b}
    pairs = [(by_a[k], by_b[k]) for k in sorted(set(by_a) & set(by_b))
             if by_a[k].get('status') == 'ok' and by_b[k].get('status') == 'ok']
    if not pairs:
        return {'n_pairs': 0, 'note': 'no item was answered validly by both arms'}

    def metric(records):
        block = evaluate(records, task, label_names, threshold, resamples=0)
        return {'t1_identity': block.get('accuracy'),
                't2_functional_groups': block.get('macro_f1_valid_subset'),
                't3_atom_counts': block.get('exact_match_accuracy'),
                't4_qed': block.get('mae'),
                't5_hydrogen_bond_donors': block.get('hbd_exact_match_accuracy')}[task]

    observed_a = metric([a for a, _ in pairs])
    observed_b = metric([b for _, b in pairs])
    observed = (observed_a - observed_b) if None not in (observed_a, observed_b) else None

    clusters = {}
    for a, b in pairs:
        for unit in clusters_of(a):
            clusters.setdefault(unit, []).append((a, b))
    keys = sorted(clusters)
    result = {
        'n_pairs': len(pairs),
        'metric': {'t1_identity': 'accuracy', 't2_functional_groups': 'macro_f1_valid_subset',
                   't3_atom_counts': 'exact_match_accuracy', 't4_qed': 'mae',
                   't5_hydrogen_bond_donors': 'hbd_exact_match_accuracy'}[task],
        'value_a': observed_a, 'value_b': observed_b, 'observed_difference': observed,
        'lower_is_better': task == 't4_qed',
        'cluster_summary': {'n_clusters': len(keys)},
    }
    if resamples <= 0 or len(keys) < 2:
        result['ci'] = None
        return result

    rng = random.Random(seed)
    diffs = []
    for _ in range(resamples):
        sample = _resample([a for a, _ in pairs], clusters, keys, rng)
        picked = {id(a) for a in sample}
        sub = [(a, b) for a, b in pairs if id(a) in picked]
        if len(sub) < 2:
            continue
        va = metric([a for a, _ in sub])
        vb = metric([b for _, b in sub])
        if va is None or vb is None:
            continue
        diffs.append(va - vb)
    if len(diffs) < max(20, resamples // 10):
        result['ci'] = None
        return result
    diffs.sort()
    alpha = (1.0 - confidence) / 2.0
    low = diffs[min(len(diffs) - 1, int(alpha * len(diffs)))]
    high = diffs[min(len(diffs) - 1, int((1.0 - alpha) * len(diffs)))]
    result['ci'] = {'low': low, 'high': high, 'resamples': len(diffs),
                    'excludes_zero': (low > 0) or (high < 0)}
    return result


def holm_adjust(p_values):
    """Holm-Bonferroni adjustment for a pre-declared family of comparisons."""
    items = [(name, p) for name, p in p_values.items() if p is not None]
    items.sort(key=lambda kv: kv[1])
    m = len(items)
    adjusted = {}
    running = 0.0
    for rank, (name, p) in enumerate(items):
        value = min(1.0, (m - rank) * p)
        running = max(running, value)
        adjusted[name] = running
    for name, p in p_values.items():
        adjusted.setdefault(name, None)
    return adjusted


def paired_metric_permutation(records_a, records_b, task, label_names, threshold=0.05,
                              resamples=2000, seed=0):
    """Cluster sign-flip permutation whose statistic IS the main-table metric.

    ``paired_comparison`` permutes a per-item score, which for a corpus-level
    metric is a different quantity: for T2 a per-item label accuracy sits near
    0.9 while Macro-F1 is near 0.5, so its p-value would answer the wrong
    question. Here each replicate swaps the WHOLE prediction vector of every
    item in the flipped clusters and then recomputes the task's own metric on
    both arms, so the null distribution is for exactly the number the main table
    prints. For T2 that means Macro-F1 is recomputed over the same fixed label
    set after the swap, never replaced by label accuracy.

    Clusters are the dependence units: for an Identity pair, BOTH endpoint
    molecules, and a pair flips only when its ownership parity is odd, which is
    what keeps shared endpoints from being treated as independent.
    """
    by_a = {pair_key(r): r for r in records_a}
    by_b = {pair_key(r): r for r in records_b}
    pairs = [(by_a[k], by_b[k]) for k in sorted(set(by_a) & set(by_b))
             if by_a[k].get('status') == 'ok' and by_b[k].get('status') == 'ok']
    if not pairs:
        return {'n_pairs': 0, 'status': 'not tested',
                'note': 'no item was answered validly by both arms'}

    def metric(records):
        block = evaluate(records, task, label_names, threshold, resamples=0)
        return {'t1_identity': block.get('accuracy'),
                't2_functional_groups': block.get('macro_f1_valid_subset'),
                't3_atom_counts': block.get('exact_match_accuracy'),
                't4_qed': block.get('mae'),
                't5_hydrogen_bond_donors': block.get('hbd_exact_match_accuracy')}[task]

    observed = metric([a for a, _ in pairs]) - metric([b for _, b in pairs])
    owners = [sorted(set(clusters_of(a)) | set(clusters_of(b))) for a, b in pairs]
    units = sorted({u for owned in owners for u in owned})
    result = {
        'n_pairs': len(pairs),
        'metric': {'t1_identity': 'accuracy', 't2_functional_groups': 'macro_f1_valid_subset',
                   't3_atom_counts': 'exact_match_accuracy', 't4_qed': 'mae',
                   't5_hydrogen_bond_donors': 'hbd_exact_match_accuracy'}[task],
        'observed_difference': observed,
        'lower_is_better': task == 't4_qed',
        'cluster_summary': {'n_clusters': len(units)},
        'method': 'cluster sign-flip permutation, whole-item prediction vectors swapped',
    }
    if resamples <= 0 or len(units) < 2:
        result['p_value'] = None
        result['status'] = 'not tested'
        result['note'] = 'too few clusters to permute'
        return result

    unit_index = {unit: i for i, unit in enumerate(units)}
    owner_idx = [[unit_index[u] for u in owned] for owned in owners]
    rng = random.Random(seed)
    extreme, counter = 0, 0
    for _ in range(resamples):
        flips = [rng.random() < 0.5 for _ in units]
        left, right = [], []
        for (a, b), idxs in zip(pairs, owner_idx):
            swap = False
            for i in idxs:
                if flips[i]:
                    swap = not swap
            left.append(b if swap else a)
            right.append(a if swap else b)
        va, vb = metric(left), metric(right)
        if va is None or vb is None:
            continue
        counter += 1
        if abs(va - vb) >= abs(observed) - 1e-12:
            extreme += 1
    result['status'] = 'tested'
    result['resamples'] = counter
    # a Monte-Carlo permutation cannot express a p below 1/(resamples+1); the
    # floor is reported rather than printing p = 0
    result['p_value'] = (extreme + 1) / (counter + 1) if counter else None
    result['p_floor'] = 1.0 / (counter + 1) if counter else None
    return result


def holm_adjust_family(p_values, family_size):
    """Holm-Bonferroni where the family size is FIXED in advance.

    ``holm_adjust`` divides by however many p-values happen to be present, so a
    comparison that could not be computed quietly makes every remaining
    comparison look stronger. That is exactly the "shrink the family after
    seeing the results" move the declared protocol forbids. Here the divisor
    comes from the declared family, so an untested hypothesis stays untested
    and the others are corrected as if it were still there.
    """
    items = [(name, p) for name, p in p_values.items() if p is not None]
    items.sort(key=lambda kv: kv[1])
    m = max(family_size, len(items))
    adjusted = {}
    running = 0.0
    for rank, (name, p) in enumerate(items):
        value = min(1.0, (m - rank) * p)
        running = max(running, value)
        adjusted[name] = running
    for name in p_values:
        adjusted.setdefault(name, None)
    return adjusted


DECLARED_COMPARISON_FAMILY = {
    'comparisons': [
        ('our_code', 'smiles'), ('our_code', 'selfies'), ('smiles', 'selfies'),
    ],
    '_rule': ('Declared BEFORE the formal run. Within one task the three pairwise '
              'comparisons form one family of three, adjusted with Holm-Bonferroni. '
              'Comparisons are run only at main scale; smoke is a pipeline check and '
              'is not tested for significance.'),
}

# The EXTENSION stage family. Declared in config/extension_config.json before any
# new request; mirrored here so a comparison can never be run against a family
# that was invented after the fact.
EXTENSION_COMPARISON_FAMILY = {
    'within_each_task': [
        ('our_code', 'smiles'), ('our_code', 'selfies'), ('our_code', 'img'),
        ('our_code', 'graph2d'), ('our_code', 'graph3d'), ('our_code', 'cml'),
        ('our_code', 'moljson'),
    ],
    'extra_pairs': [('graph2d', 'graph3d'), ('smiles', 'img')],
    'size_per_task': 9,
    'correction': 'Holm-Bonferroni within the task family, family size fixed at 9',
    'significance_level': 0.05,
}