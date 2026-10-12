"""Command line entry point for the molecular representation benchmark.

    python run_benchmark.py prepare      # build the pool, splits and manifest
    python run_benchmark.py selftest     # scoring checks with known/wrong answers
    python run_benchmark.py dry-run      # build every prompt, call no model
    python run_benchmark.py smoke        # real model calls, smoke size
    python run_benchmark.py main         # real model calls, main size
    python run_benchmark.py score --run <run_id>
    python run_benchmark.py tables --run <run_id>

Run it with the project virtual environment, which is where rdkit and selfies
live:

    "C:/Users/Asus/Desktop/Design Molecule Class/.venv/Scripts/python.exe" run_benchmark.py prepare

`smoke` and `main` refuse to start until the endpoint is configured. See
README.md for the three settings a run needs and how to set them locally.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import pathlib
import random
import sys
from collections import Counter

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

# The console here is cp950; the provider notice quotes Chinese text.
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    sys.stderr.reconfigure(encoding='utf-8', errors='replace')
except Exception:  # noqa: BLE001
    pass

from bmb import client as client_mod          # noqa: E402
from bmb import dataset, labels, paths, report, reps, runner, scoring, tasks  # noqa: E402
from bmb.core_loader import core_sha256, core_matches_baseline  # noqa: E402


# Runs whose consumed molecules are quarantined from the formal evaluation set.
DEBUG_RUNS = ('probe-01', 'probe-02', 'smoke-02')

# Identity pair construction for the EXTENSION stage: slot A is always drawn in
# view 0 and slot B in view 1, for a positive pair and a negative pair alike, so
# the two views cannot themselves signal whether the answer is "same".
# The three inherited representations are not rebuilt under this rule; their
# frozen main-02 responses are reused and the difference is declared.
EXTENSION_PAIR_VIEWS = (0, 1)


def load_local_env():
    """Load config/local.env into the process environment if it exists.

    This is a convenience for repeated local runs. The file is gitignored and
    is never read into any report, log or run record -- only the variable names
    documented in README.md appear anywhere outside the file itself.
    """
    path = paths.CONFIG_DIR / 'local.env'
    if not path.is_file():
        return
    for line in path.read_text(encoding='utf-8').splitlines():
        line = line.strip()
        if not line or line.startswith('#') or '=' not in line:
            continue
        name, _, value = line.partition('=')
        name, value = name.strip(), value.strip().strip('"').strip("'")
        if name and value and name not in os.environ:
            os.environ[name] = value


def load_json(path):
    with open(path, encoding='utf-8') as handle:
        return json.load(handle)


def load_configs():
    return (load_json(paths.MODEL_CONFIG), load_json(paths.RUN_CONFIG))


def cmd_prepare(args):
    print('reading raw ChEBI-20 splits from', paths.RAW_DIR)
    pool, funnel, exclusions = dataset.build_pool()
    run_config = load_json(paths.RUN_CONFIG)
    dev, test = dataset.split_pool(pool, run_config['split']['dev_fraction'],
                                   run_config['seed'])

    # A molecule that any debug run already answered is no longer usable for the
    # formal evaluation, because the prompt budget and the thinking mode were
    # adjusted in response to those very answers. Quarantine it rather than
    # leave it in the evaluation set.
    raw_quarantine = getattr(args, 'quarantine_runs', None)
    if raw_quarantine is None:
        debug_runs = DEBUG_RUNS
    else:
        debug_runs = tuple(name.strip() for name in raw_quarantine.split(',') if name.strip())
    quarantined_keys = runner.debug_used_clusters(debug_runs, run_config)
    quarantined = [row for row in pool if row['cluster'] in quarantined_keys]
    dev = [row for row in dev if row['cluster'] not in quarantined_keys]
    test = [row for row in test if row['cluster'] not in quarantined_keys]
    if quarantined:
        print('quarantined %d molecule(s) already used by %s' % (len(quarantined), ', '.join(debug_runs)))

    support = dataset.task_support(pool)
    print('pool=%d dev=%d test=%d' % (len(pool), len(dev), len(test)))

    # freeze the functional-group label set on measured support
    fg_config = labels.load_candidate_functional_groups()
    frozen, disabled = [], []
    for item in fg_config['labels']:
        stats = support['functional_groups'].get(item['name'], {'positive': 0, 'negative': 0})
        enough = (stats['positive'] >= fg_config['min_positives']
                  and stats['negative'] >= fg_config['min_negatives'])
        entry = dict(item, support=stats, enabled=bool(item['enabled'] and enough))
        (frozen if entry['enabled'] else disabled).append(entry)
    with open(paths.CONFIG_DIR / 'functional_groups.frozen.json', 'w', encoding='utf-8') as handle:
        json.dump({'labels': frozen, 'excluded': disabled}, handle, ensure_ascii=False, indent=2)

    manifest = {
        'dataset_version': run_config['dataset_version'] + ('-q1' if quarantined else ''),
        'base_dataset_version': run_config['dataset_version'],
        'built_at': datetime.datetime.now().isoformat(timespec='seconds'),
        'source': {
            'name': 'ChEBI-20 (MolT5 release)',
            'url': 'https://github.com/blender-nlp/MolT5/tree/main/ChEBI-20_data',
            'license': 'BSD-3-Clause (MolT5 repository)',
            'files': {name: {'sha256': dataset._sha256(paths.RAW_DIR / name)}
                      for name in sorted(os.listdir(paths.RAW_DIR)) if name.endswith('.txt')},
            'note': 'Chemical structures originate from ChEBI; this is the MolT5 redistribution of the ChEBI-20 pairing.',
        },
        'filter_rule': {
            'F1': 'RDKit parses the source SMILES',
            'F2': 'no stereochemistry descriptor in the source SMILES (excluded, never stripped)',
            'F3': 'Our Code parses it',
            'F4': 'Our Code write -> re-parse keeps atom and bond counts',
            'F5': "RDKit reads Our Code's output as the same molecule",
            'F6': 'RDKit canonical SMILES is also readable by Our Code, so all arms share one molecule list',
            'F7': 'the canonical SMILES is a fixed point under re-parse (reference is unambiguous)',
            'F8': 'SELFIES encodes it, so all three round-1 arms share exactly one molecule list',
        },
        'funnel': dict(funnel),
        'exclusions': dict(exclusions),
        'counts': {'pool': len(pool), 'dev': len(dev), 'test': len(test)},
        'splits': {
            'rule': 'molecule-level; a canonical SMILES lies wholly inside dev or test',
            'dev_fraction': run_config['split']['dev_fraction'],
            'seed': run_config['seed'],
            'provenance_statement': ('THIS PROJECT\'S OWN SPLIT. dev/test are a fresh molecule-level '
                                     'partition of a pool this project assembled. It does not '
                                     'reproduce, and must not be described as reproducing, the '
                                     'reference paper\'s or ChEBI-20\'s standard evaluation protocol.'),
            'quarantined_from_debug_runs': {
                'runs': list(debug_runs),
                'n_molecules': len(quarantined),
                'reason': ('these molecules were answered by a debug run whose observations drove '
                           'a change to max_tokens and to the thinking mode, so keeping them in '
                           'the evaluation set would make it non-independent of the tuning'),
                'clusters': sorted(quarantined_keys) if len(quarantined_keys) <= 50 else
                            sorted(quarantined_keys)[:50] + ['...'],
            },
        },
        'task_support': support,
        'toolchain': {
            'rdkit': labels.rdkit_version(),
            'python': sys.version.split()[0],
            'qed_weights': labels.qed_weights(),
        },
        'scoring_rules': {
            'qed_success_threshold': run_config['qed_success_threshold'],
            'macro_f1_penalty': 'invalid answer scored as all-false (fixed, declared)',
            'missing_answer': 'failure in all-sample rates; excluded from conditional errors with n stated',
            'bootstrap': run_config['bootstrap'],
        },
    }
    import hashlib
    for name, rows in (('pool', pool), ('dev', dev), ('test', test)):
        manifest['splits']['%s_smiles_sha256' % name] = hashlib.sha256(
            '\n'.join(r['smiles'] for r in rows).encode('utf-8')).hexdigest()
    dataset.write_prepared(pool, dev, test, funnel, exclusions, manifest,
                           quarantined=quarantined)
    print('wrote', paths.prepared('dataset_manifest.json'))
    print('enabled functional groups:', [i['name'] for i in frozen])
    return 0


def cmd_selftest(args):
    """Scoring checks with known answers and deliberately wrong outputs.

    These outputs are test fixtures, not model results, and are written to the
    selftest folder. They never enter a benchmark table.
    """
    run_config = load_json(paths.RUN_CONFIG)
    failures = []
    total = []

    def _refused(fn):
        """True when the call raised a refusal, False when it went through."""
        try:
            fn()
            return False
        except (client_mod.EndpointNotPermitted, client_mod.ModelConfigMissing):
            return True

    def check(name, got, want):
        total.append(name)
        ok = got == want
        print('  %-56s %s' % (name, 'ok' if ok else 'FAIL got=%r want=%r' % (got, want)))
        if not ok:
            failures.append(name)

    label_names = [item['name'] for item in labels.enabled_labels()]

    def rec(task, parsed, label, status='ok', cluster='c', clusters=None,
            example_id=None):
        record = {'task': task, 'parsed': parsed, 'label': label, 'status': status,
                  'cluster': cluster,
                  'example_id': example_id or ('%s-%s' % (task, cluster))}
        if clusters is not None:
            record['clusters'] = list(clusters)
        return record

    # T1: two right, two wildly wrong (both off-diagonal cells)
    rows = [rec('t1_identity', {'same_molecule': True}, {'same_molecule': True}),
            rec('t1_identity', {'same_molecule': False}, {'same_molecule': False}),
            rec('t1_identity', {'same_molecule': True}, {'same_molecule': False}),
            rec('t1_identity', {'same_molecule': False}, {'same_molecule': True})]
    block = scoring.evaluate(rows, 't1_identity', label_names, 0.05, resamples=0)
    check('t1 accuracy on 2/4', block['accuracy'], 0.5)
    check('t1 balanced accuracy on 2/4', block['balanced_accuracy'], 0.5)
    check('t1 confusion', block['confusion'], {'tp': 1, 'fp': 1, 'tn': 1, 'fn': 1})

    # T1: an invalid answer must not be counted as correct or as a class default
    rows = [rec('t1_identity', None, {'same_molecule': True}, status='invalid')]
    block = scoring.evaluate(rows, 't1_identity', label_names, 0.05, resamples=0)
    check('t1 invalid answer leaves n_scored at 0', block['n_scored'], 0)
    check('t1 invalid answer counted, valid rate 0', block['valid_answer_rate'], 0.0)

    # T2: perfect prediction scores 1.0; a wrong-but-valid answer scores 0.0
    perfect = [rec('t2_functional_groups', {n: True for n in label_names},
                   {n: True for n in label_names})]
    block = scoring.evaluate(perfect, 't2_functional_groups', label_names, 0.05, resamples=0)
    check('t2 perfect macro-F1 (valid subset)', block['macro_f1_valid_subset'], 1.0)
    inverted = [rec('t2_functional_groups', {n: False for n in label_names},
                    {n: True for n in label_names})]
    block = scoring.evaluate(inverted, 't2_functional_groups', label_names, 0.05, resamples=0)
    check('t2 valid but all-wrong macro-F1', block['macro_f1_valid_subset'], 0.0)

    # an invalid answer is penalised in the all-sample variant, and the
    # valid-subset variant is left undefined rather than invented
    blank = [rec('t2_functional_groups', None, {n: True for n in label_names}, status='invalid')]
    block = scoring.evaluate(blank, 't2_functional_groups', label_names, 0.05, resamples=0)
    check('t2 no valid answers -> valid-subset metric undefined',
          block['macro_f1_valid_subset'], None)
    check('t2 no valid answers -> penalised metric 0', block['macro_f1_all_samples_penalised'], 0.0)

    # T3: partial correctness must not count as exact match
    rows = [rec('t3_atom_counts', {'C': 1, 'N': 0, 'O': 9, 'F': 0}, {'C': 1, 'N': 0, 'O': 2, 'F': 0})]
    block = scoring.evaluate(rows, 't3_atom_counts', label_names, 0.05, resamples=0)
    check('t3 one wrong element -> not exact', block['exact_match_accuracy'], 0.0)
    check('t3 per-element still counted', block['per_element_accuracy']['O']['correct'], 0)

    # T4: an out-of-range answer must not be clipped into a small error
    rows = [rec('t4_qed', {'qed': 1.0}, {'qed': 0.0}),
            rec('t4_qed', None, {'qed': 0.0}, status='invalid')]
    block = scoring.evaluate(rows, 't4_qed', label_names, 0.05, resamples=0)
    check('t4 mae over valid only', round(block['mae'], 6), 1.0)
    check('t4 n_valid/n_total stated', (block['n_valid'], block['n_total']), (1, 2))
    check('t4 success rate uses all samples', block['success_rate_all_samples'], 0.0)

    # Spearman: strictly increasing -> 1.0, strictly decreasing -> -1.0
    increasing = [(1, 10), (2, 20), (3, 30), (4, 40)]
    decreasing = [(1, 40), (2, 30), (3, 20), (4, 10)]
    check('spearman increasing', round(scoring.spearman(increasing), 6), 1.0)
    check('spearman decreasing', round(scoring.spearman(decreasing), 6), -1.0)
    # ties must not fake a perfect correlation
    check('spearman with a tie is not 1.0',
          scoring.spearman([(1, 5), (1, 6), (2, 7), (3, 9)]) < 1.0, True)

    # ---- interval behaviour on small samples -----------------------------
    def identity_block(correct_flags):
        rows = []
        for index, correct in enumerate(correct_flags):
            rows.append(rec('t1_identity', {'same_molecule': correct},
                            {'same_molecule': True}, cluster='m%d' % index,
                            clusters=['m%d' % index]))
        return scoring.evaluate(rows, 't1_identity', label_names, 0.05, resamples=400)

    block = identity_block([True] * 8)
    check('8/8 does not claim certainty (bootstrap flagged degenerate)',
          block['ci95_cluster_bootstrap']['degenerate'], True)
    check('8/8 Wilson lower bound is well below 1',
          block['ci95_wilson']['low'] < 0.75, True)

    block = identity_block([False] * 8)
    check('0/8 Wilson upper bound is well below 1',
          block['ci95_wilson']['high'] < 0.45, True)

    block = identity_block([True] * 6 + [False] * 2)
    check('mixed result gives a non-degenerate bootstrap',
          block['ci95_cluster_bootstrap']['degenerate'], False)
    check('mixed result intervals bracket the point estimate',
          block['ci95_wilson']['low'] <= 0.75 <= block['ci95_wilson']['high']
          and block['ci95_clopper_pearson']['low'] <= 0.75, True)

    # ---- dependence detection -------------------------------------------
    shared = [rec('t1_identity', {'same_molecule': True}, {'same_molecule': True},
                  clusters=['mA', 'mB']),
              rec('t1_identity', {'same_molecule': True}, {'same_molecule': True},
                  clusters=['mB', 'mC'])]
    summary = scoring.cluster_summary(shared)
    check('a pair is attached to both endpoints', scoring.clusters_of(shared[0]),
          ('mA', 'mB'))
    check('shared endpoints mean the items are not independent singletons',
          summary['items_are_independent_singletons'], False)
    check('shared endpoints are counted as items that own more than one cluster',
          summary['n_items_sharing_a_cluster'], len(shared))
    check('a positive pair repeating one molecule counts as one cluster',
          scoring.clusters_of(rec('t1_identity', {}, {}, clusters=['mX', 'mX'])),
          ('mX',))

    # a pair enters a resample when EITHER endpoint is drawn, a single-molecule
    # item only when its own molecule is drawn
    class _FixedRng:
        def randrange(self, _n):
            return 0

    pair = rec('t1_identity', {}, {}, clusters=['m0', 'm1'])
    single = rec('t3_atom_counts', {}, {}, clusters=['m1'])
    drawn = scoring._resample([pair, single], {'m0': [pair], 'm1': [pair, single]},
                              ['m0', 'm1'], _FixedRng())
    check('pair included when only its first endpoint is drawn',
          pair in drawn, True)
    check('single-molecule item excluded when its cluster is not drawn',
          single in drawn, False)

    # every item must have the SAME inclusion probability regardless of how
    # many endpoints it owns, or the bootstrap over-weights the items with more
    _many = []
    for i in range(60):
        _many.append(rec('t1_identity', {'same_molecule': True}, {'same_molecule': True},
                         clusters=['p%d' % i], example_id='E1-%d' % i))
    for i in range(60):
        _many.append(rec('t1_identity', {'same_molecule': True}, {'same_molecule': True},
                         clusters=['n%da' % i, 'n%db' % i], example_id='E2-%d' % i))
    _idx = {}
    for r in _many:
        for u in scoring.clusters_of(r):
            _idx.setdefault(u, []).append(r)
    _keys = sorted(_idx)
    _rng = random.Random(11)
    _counts = {'one': 0, 'two': 0}
    for _ in range(300):
        for r in scoring._resample(_many, _idx, _keys, _rng):
            _counts['one' if len(scoring.clusters_of(r)) == 1 else 'two'] += 1
    _p1 = _counts['one'] / 300 / 60
    _p2 = _counts['two'] / 300 / 60
    check('one-endpoint and two-endpoint items are drawn at the same rate',
          abs(_p1 - _p2) < 0.05, True)
    check('that rate is the single-cluster inclusion probability, not the double one',
          abs(_p2 - (1 - (1 - 1 / len(_keys)) ** len(_keys))) < 0.05, True)

    # ---- undefined statistics -------------------------------------------
    check('spearman on a constant array is undefined, not 0',
          scoring.spearman([(1, 5), (2, 5), (3, 5), (4, 5)]), None)
    check('spearman on too few points is undefined',
          scoring.spearman([(1, 2)]), None)
    empty = scoring.evaluate([rec('t4_qed', None, {'qed': 0.5}, status='invalid')],
                             't4_qed', label_names, 0.05, resamples=0)
    check('MAE with no valid answer is undefined, not 0', empty['mae'], None)
    check('success rate with no valid answer is a real 0',
          empty['success_rate_all_samples'], 0.0)

    # ---- paired comparison ----------------------------------------------
    same = [rec('t1_identity', {'same_molecule': True}, {'same_molecule': True},
                cluster='m%d' % i, clusters=['m%d' % i]) for i in range(12)]
    result = scoring.paired_comparison(same, list(same), 't1_identity',
                                       resamples=300, seed=1)
    check('identical arms give a zero paired difference',
          result['observed_mean_difference'], 0.0)
    check('identical arms give a p-value near 1', result['p_value'] > 0.5, True)

    worse = [rec('t1_identity', {'same_molecule': False}, {'same_molecule': True},
                 cluster='m%d' % i, clusters=['m%d' % i]) for i in range(12)]
    result = scoring.paired_comparison(same, worse, 't1_identity',
                                       resamples=300, seed=1)
    check('a clearly worse arm is detected', result['p_value'] < 0.05, True)
    check('paired test reports how many items were dropped as invalid',
          'n_unpaired_because_invalid' in result, True)

    # example_id embeds the representation, so a naive join across two arms
    # finds ZERO pairs and the comparison silently reports nothing. pair_key
    # must strip it.
    a = rec('t3_atom_counts', {'C': 1, 'N': 0, 'O': 0, 'S': 0},
            {'C': 1, 'N': 0, 'O': 0, 'S': 0}, cluster='123',
            example_id='t3-smiles-123')
    a['representation'] = 'smiles'
    b = rec('t3_atom_counts', {'C': 1, 'N': 0, 'O': 0, 'S': 0},
            {'C': 1, 'N': 0, 'O': 0, 'S': 0}, cluster='123',
            example_id='t3-our_code-123')
    b['representation'] = 'our_code'
    check('two arms answering the same question share a pair key',
          scoring.pair_key(a) == scoring.pair_key(b), True)
    check('the raw example ids really do differ', a['example_id'] != b['example_id'], True)
    dup = scoring.paired_comparison([a], [b], 't3_atom_counts', resamples=0)
    check('the paired test finds the pair rather than reporting none',
          dup['n_pairs'], 1)
    idt = rec('t1_identity', {'same_molecule': True}, {'same_molecule': True},
              example_id='t1-smiles-999-pos')
    idt['representation'] = 'smiles'
    idt2 = rec('t1_identity', {'same_molecule': True}, {'same_molecule': True},
               example_id='t1-our_code-999-pos')
    idt2['representation'] = 'our_code'
    check('a positive pair keeps its -pos suffix in the key',
          scoring.pair_key(idt).endswith('-pos')
          and scoring.pair_key(idt) == scoring.pair_key(idt2), True)

    # ---- budget guard ----------------------------------------------------
    guard = runner.Budget(max_basic_requests=3, max_transport_attempts=4,
                          max_total_tokens=None, in_flight_reserve=0)
    check('budget allows requests below the cap', guard.can_continue(), True)
    guard.note_basic_request()
    guard.record([], {'prompt_tokens': 80, 'completion_tokens': 20, 'total_tokens': 100})
    guard.note_basic_request()
    guard.record([{'attempt': 1}, {'attempt': 2}],   # one request, two transport attempts
                 {'prompt_tokens': 40, 'completion_tokens': 10, 'total_tokens': 50})
    check('basic requests are counted separately from attempts', guard.basic_requests, 2)
    check('a retry costs an extra transport attempt', guard.transport_attempts, 3)
    check('prompt and completion tokens are tracked apart',
          (guard.prompt_tokens, guard.completion_tokens), (120, 30))
    check('total tokens is the sum of the parts', guard.total_tokens, 150)
    guard.note_basic_request()
    guard.record([], {'prompt_tokens': 10, 'completion_tokens': 5, 'total_tokens': 15})
    check('budget refuses once the basic-request cap is reached',
          guard.can_continue(), False)
    check('budget names which cap stopped it',
          'cap reached' in (guard.stopped_reason or ''), True)
    only_basic = runner.Budget(max_basic_requests=1, max_transport_attempts=None,
                               in_flight_reserve=0)
    only_basic.note_basic_request()
    only_basic.record([], {'total_tokens': 1})
    only_basic.can_continue()
    check('with only a basic cap set, that is the reason given',
          'basic request cap' in (only_basic.stopped_reason or ''), True)

    attempt_guard = runner.Budget(max_basic_requests=None, max_transport_attempts=2,
                                  in_flight_reserve=0)
    attempt_guard.record([{'attempt': 1}, {'attempt': 2}], {'total_tokens': 10})
    check('transport-attempt cap stops a retry-heavy run',
          attempt_guard.can_continue(), False)

    # missing usage must never be charged as zero
    quiet = runner.Budget(max_basic_requests=None, in_flight_reserve=0)
    quiet.record([], None)
    check('a response with no usage is charged a reserve, not zero',
          quiet.total_tokens > 0, True)
    check('the unknown-usage count is recorded', quiet.attempts_without_usage, 1)

    # in-flight reserve
    reserved = runner.Budget(max_transport_attempts=2, in_flight_reserve=1)
    reserved.record([], {'total_tokens': 1})
    check('in-flight reserve keeps room for a request already sent',
          reserved.can_continue(), False)

    # RMB cap, priced at the busy-hour rate
    priced = runner.Budget(max_rmb=0.01, in_flight_reserve=0,
                           prices={'input_busy': 2.0, 'output_busy': 8.0})
    priced.record([], {'prompt_tokens': 10000, 'completion_tokens': 1000,
                       'total_tokens': 11000})
    check('RMB cap stops the run', priced.can_continue(), False)
    check('RMB estimate uses the pessimistic (busy) rate',
          round(priced.estimated_rmb(), 4), round(10000 / 1e6 * 2 + 1000 / 1e6 * 8, 4))

    # a subscription provider has no per-token RMB price: the amount is
    # UNKNOWN, which is None, never 0.0 (0.0 would read as "cost nothing")
    sub = runner.Budget(max_rmb=30.0, in_flight_reserve=0, prices=None)
    check('a subscription provider yields no RMB figure', sub.estimated_rmb(), None)
    sub.record([], {'prompt_tokens': 5000, 'completion_tokens': 500, 'total_tokens': 5500})
    check('and its RMB cap cannot stop the run on a None value', sub.can_continue(), True)
    check('no price is assumed for it', not sub.prices, True)

    # ---- endpoint licensing guard ----------------------------------------
    for blocked in ('https://token-plan.cn-beijing.maas.aliyuncs.com/compatible-mode/v1',
                    'https://ws1.cn-beijing.maas.aliyuncs.com/compatible-mode/v1'):
        try:
            client_mod.check_endpoint_permitted(blocked)
            check('blocked endpoint %s is refused' % blocked[:28], 'allowed', 'refused')
        except client_mod.EndpointNotPermitted:
            check('blocked endpoint %s is refused' % blocked[:28], 'refused', 'refused')
    try:
        client_mod.check_endpoint_permitted('https://dashscope.aliyuncs.com/compatible-mode/v1')
        check('the standard DashScope host is allowed', 'allowed', 'allowed')
    except client_mod.EndpointNotPermitted:
        check('the standard DashScope host is allowed', 'refused', 'allowed')
    try:
        client_mod.ModelClient(base_url='https://token-plan.cn-beijing.maas.aliyuncs.com/compatible-mode/v1',
                               model_id='deepseek-v4.1-flash')
        check('the client itself refuses a Token Plan endpoint', 'constructed', 'refused')
    except client_mod.EndpointNotPermitted:
        check('the client itself refuses a Token Plan endpoint', 'refused', 'refused')

    check('experiment credentials are read before any harness variable',
          client_mod.EXPERIMENT_API_KEY_ENVS[0], 'BMB_API_KEY')

    # ---- provider selection ---------------------------------------------
    PLAN = 'https://token-plan.cn-beijing.maas.aliyuncs.com/compatible-mode/v1'
    STD = 'https://dashscope.aliyuncs.com/compatible-mode/v1'

    check('bailian provider still refuses a plan host',
          _refused(lambda: client_mod.check_endpoint_permitted(PLAN, 'bailian')), True)
    check('bailian provider still allows the standard host',
          _refused(lambda: client_mod.check_endpoint_permitted(STD, 'bailian')), False)
    check('token_plan provider allows its own host',
          _refused(lambda: client_mod.check_endpoint_permitted(PLAN, 'token_plan')), False)
    check('token_plan provider refuses a mismatched standard host',
          _refused(lambda: client_mod.check_endpoint_permitted(STD, 'token_plan')), True)
    check('an unknown provider is rejected',
          _refused(lambda: client_mod.check_endpoint_permitted(STD, 'nope')), True)
    check('a plain-HTTP endpoint is refused',
          _refused(lambda: client_mod.check_endpoint_permitted('http://dashscope.aliyuncs.com/compatible-mode/v1')), True)

    check('token_plan reads only its own credential variable',
          client_mod.PROVIDERS['token_plan']['credential_envs'],
          ('BMB_TOKEN_PLAN_API_KEY',))
    check('the plan credential variable is not in the standard provider list',
          'BMB_TOKEN_PLAN_API_KEY' not in client_mod.PROVIDERS['bailian']['credential_envs'], True)
    check('no harness credential variable appears in any provider',
          all('ANTHROPIC' not in name
              for spec in client_mod.PROVIDERS.values()
              for name in spec['credential_envs']), True)

    plan_client = client_mod.ModelClient(base_url=PLAN, model_id='deepseek-v4.1-flash',
                                         provider='token_plan')
    check('a token_plan client resolves the plan credential variable',
          'BMB_TOKEN_PLAN_API_KEY' in plan_client.credential_envs(), True)
    check('the plan provider carries a terms warning to record',
          bool(client_mod.provider_terms_warning('token_plan')), True)
    check('the standard provider carries no terms warning',
          client_mod.provider_terms_warning('bailian'), None)
    check('the warning does not claim the use is permitted',
          'DOES NOT CLAIM' in (client_mod.provider_terms_warning('token_plan') or ''), True)

    # ---- complexity gate -------------------------------------------------
    def mol(smiles):
        return {'cid': smiles, 'smiles': smiles, 'cluster': smiles,
                'description': '', 'formula': '', 'heavy_atoms': 0,
                'source_split': 'test'}

    kept, dropped = dataset.apply_complexity_filter(
        [mol('CCO'), mol('c1ccccc1'), mol('CC(=O)Oc1ccccc1C(=O)O')])
    check('ethanol is dropped by the complexity gate',
          any(r['smiles'] == 'CCO' for r in dropped), True)
    check('benzene is dropped, and for the weight rule first (78 Da)',
          any(r['smiles'] == 'c1ccccc1'
              and r['drop_reason'] == 'below_min_molecular_weight' for r in dropped), True)
    _long, long_dropped = dataset.apply_complexity_filter([mol('CCCCCCCCCCCCCCCCCCCC(=O)O')])
    check('a 312 Da acyclic acid is dropped specifically for the ring rule',
          long_dropped and long_dropped[0]['drop_reason'] == 'below_min_rings', True)
    check('aspirin (180 Da, 1 ring) is also dropped',
          any(r['smiles'] == 'CC(=O)Oc1ccccc1C(=O)O' for r in dropped), True)
    check('the gate reports both drop reasons separately',
          {r['drop_reason'] for r in dropped} <= {'below_min_molecular_weight',
                                                  'below_min_rings'}, True)

    # the frozen main item list must be reproducible from the FROZEN MAIN POOL.
    # Building it from the raw test split instead produces a different item set
    # (they overlapped by only 12 of 200 when this bug was live), which would
    # silently drop the complexity gate from a "main" run.
    frozen = [json.loads(l) for l in
              (paths.prepared('main_items') / 'smiles-t1_identity.jsonl')
              .read_text(encoding='utf-8').splitlines() if l.strip()]
    main_pool = dataset.load_prepared('main_molecules')
    test_pool = dataset.load_prepared('test')
    check('the frozen main pool is smaller than the raw test split',
          len(main_pool) < len(test_pool), True)
    rebuilt_from_main = tasks.get('t1_identity').build_items(
        main_pool, 200, run_config['identity_seed'], 'smiles', 'test')
    rebuilt_from_test = tasks.get('t1_identity').build_items(
        test_pool, 200, run_config['identity_seed'], 'smiles', 'test')
    frozen_ids = {i['example_id'] for i in frozen}
    check('frozen main items reproduce from the frozen main pool',
          {i['example_id'] for i in rebuilt_from_main} == frozen_ids, True)
    check('and NOT from the raw test split (the bug that was fixed)',
          len({i['example_id'] for i in rebuilt_from_test} & frozen_ids) < len(frozen_ids),
          True)
    check('every rebuilt main item uses only gated molecules',
          all(dataset.molecular_weight(m['smiles']) >= 300
              and dataset.ring_count(m['smiles']) >= 2
              for i in rebuilt_from_main
              for m in main_pool
              if m['cluster'] in (i.get('clusters') or (i['cluster'],))), True)

    # the frozen-list verifier must REFUSE a tampered manifest before any request
    import copy as _copy
    _mf = json.loads(paths.DELIVERABLES.joinpath('main_manifest.json')
                     .read_text(encoding='utf-8')) if hasattr(paths, 'DELIVERABLES') else None
    _mf = json.loads((pathlib.Path(__file__).resolve().parent.parent /
                      '08_deliverables' / 'main_manifest.json').read_text(encoding='utf-8'))

    def _frozen_refused(manifest, rep='smiles', task_name='t1_identity'):
        try:
            runner.load_frozen_items(rep, task_name, manifest, run_config=run_config)
            return False
        except runner.FrozenListError:
            return True

    check('the genuine frozen list verifies', _frozen_refused(_mf), False)
    _bad = _copy.deepcopy(_mf)
    _bad['item_list_hashes']['t1_identity']['smiles'] = '0' * 64
    check('a tampered item-list hash is refused', _frozen_refused(_bad), True)
    _bad2 = _copy.deepcopy(_mf)
    _bad2['prompt_hashes']['smiles']['t1_identity'] = ['0' * 64]
    check('a wrong prompt hash set is refused', _frozen_refused(_bad2), True)
    check('a manifest with no recorded list is refused', _frozen_refused({}), True)
    _bad3 = _copy.deepcopy(_mf)
    _bad3['item_list_hashes'].pop('t5_hydrogen_bond_donors', None)
    check('a cell missing from the manifest is refused',
          _frozen_refused(_bad3, task_name='t5_hydrogen_bond_donors'), True)
    check('every frozen list yields exactly the planned item count',
          all(len(runner.load_frozen_items(rep, t, _mf, run_config=run_config)) == 200
              for rep in reps.ROUND1_REPRESENTATIONS for t in tasks.TASK_ORDER), True)

    # T2 must average over the FROZEN label set, not over whatever labels the
    # sample happens to exercise. Re-deriving admits labels below the declared
    # support floor (nitrile at 7 positives) and changes the denominator.
    _frozen10 = labels.frozen_macro_f1_labels()
    check('the frozen T2 label set is recorded in the manifest', bool(_frozen10), True)
    check('and it is the 10-label set, not the 13-label candidate list',
          len(_frozen10 or []), 10)
    check('labels below the declared support floor are excluded',
          all(n not in (_frozen10 or []) for n in ('aldehyde', 'nitrile', 'nitro')), True)
    _src = labels.load_functional_groups()['labels']
    check('the frozen set is a subset of the operative labels',
          set(_frozen10 or []) <= {i['name'] for i in _src if i.get('enabled')}, True)

    big = dataset.load_prepared('main_molecules')
    check('every frozen main molecule passes the gate',
          all(dataset.molecular_weight(r['smiles']) >= 300.0
              and dataset.ring_count(r['smiles']) >= 2 for r in big), True)
    check('the ring threshold is inclusive: 1 ring passes a gate of 1',
          len(dataset.apply_complexity_filter([mol('c1ccccc1')], min_mw=0, min_rings=1)[0]), 1)
    check('and the same molecule fails a gate of 2',
          len(dataset.apply_complexity_filter([mol('c1ccccc1')], min_mw=0, min_rings=2)[1]), 1)

    print('\n%d checks, %d failed' % (len(total), len(failures)))
    return 1 if failures else 0


def cmd_offline_checks(args):
    """Client-side checks that need no network and spend no tokens."""
    model_config, run_config = load_configs()
    failures = []
    total = []

    def _refused_offline(fn):
        try:
            fn()
            return False
        except (client_mod.EndpointNotPermitted, client_mod.ModelConfigMissing):
            return True

    def check(name, got, want):
        total.append(name)
        ok = got == want
        print('  %-58s %s' % (name, 'ok' if ok else 'FAIL got=%r want=%r' % (got, want)))
        if not ok:
            failures.append(name)

    # The client must refuse to exist without an endpoint, rather than guess one
    try:
        client_mod.ModelClient(base_url=None, model_id=None)
        check('client refuses an unconfigured endpoint', 'constructed', 'refused')
    except client_mod.ModelConfigMissing:
        check('client refuses an unconfigured endpoint', 'refused', 'refused')

    # Request shape, built locally
    model = client_mod.ModelClient(base_url='https://example.invalid/compatible-mode/v1',
                                   model_id='example-model', region='beijing')
    payload = model._payload([{'role': 'user', 'content': 'reply with JSON'}])
    check('payload has the model id', payload['model'], 'example-model')
    check('payload is non-streaming', payload['stream'], False)
    check('payload requests json_object', payload['response_format'], {'type': 'json_object'})
    check('payload does not send tools', 'tools' in payload, False)
    check('endpoint joins correctly', model._endpoint(),
          'https://example.invalid/compatible-mode/v1/chat/completions')

    # Redaction must mask a key even when it appears inside a larger string
    import os
    os.environ['BMB_SELFTEST_KEY'] = 'sk-offline-selftest-dummy-000000'
    sample = 'Authorization: Bearer sk-offline-selftest-dummy-000000 (code invalid_api_key)'
    masked = client_mod.redact(sample)
    check('redact masks a bearer token', 'sk-offline-selftest-dummy-000000' in masked, False)
    check('redact keeps the diagnostic tail', 'invalid_api_key' in masked, True)
    del os.environ['BMB_SELFTEST_KEY']

    # The key must never appear in a stored record
    for name in ('raw_responses.jsonl',):
        pass
    check('api key variable name is declared, not its value',
          model_config['api_key_env'], 'BMB_API_KEY')
    check('the experiment does not read the harness credential first',
          client_mod.EXPERIMENT_API_KEY_ENVS[0] == 'BMB_API_KEY'
          and 'ANTHROPIC_AUTH_TOKEN' not in client_mod.EXPERIMENT_API_KEY_ENVS, True)
    active = model_config.get('active') or {}
    active_provider = active.get('provider') or 'bailian'
    check('the active provider is one this client defines',
          active_provider in client_mod.PROVIDERS, True)
    check('the active endpoint matches the active provider',
          not _refused_offline(lambda: client_mod.check_endpoint_permitted(
              active.get('base_url'), active_provider)), True)
    check('the other provider cannot use this endpoint by accident',
          _refused_offline(lambda: client_mod.check_endpoint_permitted(
              active.get('base_url'),
              'bailian' if active_provider == 'token_plan' else 'token_plan')), True)
    check('the active provider reads only its own credential variable',
          client_mod.PROVIDERS[active_provider]['credential_envs'][0]
          in (client_mod.PROVIDERS[active_provider]['credential_envs']), True)

    print('\n%d offline checks, %d failed' % (len(total), len(failures)))
    return 1 if failures else 0


def decode_payload(representation, payload):
    """Decode a representation payload back to a canonical SMILES.

    This is a CONVERSION CONSISTENCY check, not a model result. It proves the
    pipeline is lossless. It is never entered into a benchmark table: the whole
    point of T4 is that the model must produce QED itself rather than having
    the value recomputed externally.
    """
    from rdkit import Chem
    if representation == 'smiles':
        return reps.canonical_smiles(Chem.MolFromSmiles(payload))
    if representation == 'selfies':
        import selfies as sf
        return reps.canonical_smiles(Chem.MolFromSmiles(sf.decoder(payload)))
    if representation == 'our_code':
        from bmb.core_loader import load_core
        core = load_core()
        return reps.canonical_smiles(
            Chem.MolFromSmiles(core.Molecule.from_code(payload).decode_to_smiles()))
    raise ValueError('no decoder for %r' % representation)


def cmd_verify(args):
    """Round-trip every pool molecule through every round-1 representation."""
    from collections import Counter
    pool = dataset.load_prepared('pool')
    out = {}
    for representation in (args.representations.split(',') if args.representations
                           else list(reps.ROUND1_REPRESENTATIONS)):
        ok = 0
        reasons = Counter()
        for entry in pool:
            try:
                payload = reps.render(representation, entry['smiles'])
                decoded = decode_payload(representation, payload)
            except Exception as exc:  # noqa: BLE001
                reasons['%s' % type(exc).__name__] += 1
                continue
            if decoded == entry['smiles']:
                ok += 1
            else:
                reasons['decoded_to_a_different_molecule'] += 1
        out[representation] = {
            'molecules': len(pool),
            'round_trip_ok': ok,
            'coverage': ok / len(pool) if pool else None,
            'failures': dict(reasons),
        }
        print('%-10s round-trip %4d/%4d  coverage=%.3f  %s'
              % (representation, ok, len(pool), ok / len(pool) if pool else 0,
                 dict(reasons) or ''))

    # Identity pairs: the positive must decode to one molecule, the negative to two
    pair_report = {}
    for representation in out:
        positives = negatives = pos_ok = neg_ok = 0
        task = tasks.get('t1_identity')
        items = task.build_items(pool, 64, 999, representation, 'pool')
        for item in items:
            a = decode_payload(representation, item['payload']['a'])
            b = decode_payload(representation, item['payload']['b'])
            if item['label']['same_molecule']:
                positives += 1
                pos_ok += int(a == b)
            else:
                negatives += 1
                neg_ok += int(a != b)
        pair_report[representation] = {
            'positives': positives, 'positives_same': pos_ok,
            'negatives': negatives, 'negatives_different': neg_ok,
            'positive_variant_texts_differ': _variant_texts_differ(pool, representation),
        }
        print('%-10s identity pairs: pos %d/%d same, neg %d/%d different, texts differ=%s'
              % (representation, pos_ok, positives, neg_ok, negatives,
                 pair_report[representation]['positive_variant_texts_differ']))
    out['_identity_pairs'] = pair_report

    paths.DELIV_DIR.mkdir(parents=True, exist_ok=True)
    path = paths.DELIV_DIR / 'conversion_check.json'
    with open(path, 'w', encoding='utf-8') as handle:
        json.dump(out, handle, ensure_ascii=False, indent=2)
    print('wrote', path)
    return 0


def _variant_texts_differ(pool, representation):
    for entry in pool:
        try:
            a = reps.render(representation, entry['smiles'], variant=0)
            b = reps.render(representation, entry['smiles'], variant=1)
        except Exception:  # noqa: BLE001
            continue
        if a != b:
            return True
    return False


def _require_model(model_config):
    active = model_config.get('active') or {}
    env = model_config.get('env_overrides') or {}
    resolved = {
        'region': os.environ.get(env.get('region', ''), '') or active.get('region'),
        'base_url': os.environ.get(env.get('base_url', ''), '') or active.get('base_url'),
        'model_id': os.environ.get(env.get('model_id', ''), '') or active.get('model_id'),
    }
    missing = [k for k, v in resolved.items() if not v]
    if '{WorkspaceId}' in str(resolved.get('base_url') or ''):
        missing.append('base_url (still contains the {WorkspaceId} placeholder)')
    if missing:
        print('This run cannot start. Not configured yet: %s' % ', '.join(missing))
        print('Set them in config/model_config.json under "active", or export '
              'BMB_REGION / BMB_BASE_URL / BMB_MODEL_ID. See README.md.')
        return None
    active = model_config.get('active') or {}
    provider = active.get('provider', 'bailian')
    envs = client_mod.credential_env_names_for(provider)
    if not any(os.environ.get(name) for name in envs):
        print('This run cannot start. No credential for provider %r is set.' % provider)
        print('Set one of: %s' % ', '.join(envs))
        print('Set it in this shell only; do not put it in a file this project tracks '
              'and do not paste it into a chat. See README.md.')
        return None
    return resolved


def cmd_run(args, mode):
    model_config, run_config = load_configs()
    if args.representations:
        representations = args.representations.split(',')
    else:
        representations = list(reps.ROUND1_REPRESENTATIONS)
    task_names = args.tasks.split(',') if args.tasks else list(tasks.TASK_ORDER)
    mode_cfg = run_config[mode]
    run_id = args.run_id or '%s-%s' % (mode, datetime.datetime.now().strftime('%Y%m%d-%H%M%S'))
    manifest = dataset.load_manifest()

    resolved = _require_model(model_config)
    if resolved is None:
        return 2
    model = client_mod.client_from_config(model_config, run_config)
    print('run_id=%s provider=%s model=%s region=%s'
          % (run_id, model.provider, model.model_id, model.region))
    print('base_url=%s' % model.base_url)
    print('credential variable in use: %s' % (model.api_key_var_in_use or 'recognised'))
    warning = client_mod.provider_terms_warning(model.provider)
    if warning:
        print()
        print('!' * 78)
        print('NOTICE -- this project does not claim the selected provider is permitted')
        print('for this use. The platform states, verbatim:')
        print('  %s' % client_mod.TOKEN_PLAN_TERMS_CN)
        print('  source: %s' % client_mod.TOKEN_PLAN_TERMS_SOURCE)
        print('The choice is the user\'s and is recorded in the run manifest.')
        print('!' * 78)
        print()

    budget = mode_cfg.get('per_task', 8)
    if mode == 'smoke':
        budget = min(budget, run_config['request_budget']['smoke_hard_cap'])
    # main MUST run on the frozen complexity-gated molecule list, not the raw split.
    # The two are different lists; using the split silently drops the gate.
    pool_name = 'main_molecules' if mode == 'main' else args.split
    frozen_manifest = None
    if mode == 'main':
        pool_path = paths.prepared('main_molecules.jsonl')
        manifest_path = paths.DELIV_DIR / 'main_manifest.json'
        if not pool_path.is_file() or not manifest_path.is_file():
            print('This run cannot start: the frozen main molecule list or its manifest '
                  'is missing. Run `main-prep` first.')
            print('  expected: %s' % pool_path)
            print('  expected: %s' % manifest_path)
            return 2
        # the frozen item hashes live in main_manifest.json, NOT in the dataset
        # manifest that the rest of the runner uses for counts
        frozen_manifest = load_json(manifest_path)
        frozen_count = len(dataset.load_prepared('main_molecules'))
        print('main pool: frozen complexity-gated set, %d molecules (NOT the raw split)'
              % frozen_count)
        print('frozen manifest: %s' % manifest_path)
    guard = runner.budget_from_config(
        run_config, (model_config.get('active') or {}).get('provider', 'bailian'))
    print('budget guard: basic<=%s attempts<=%s tokens<=%s completion<=%s CNY<=%s'
          % (guard.max_basic_requests, guard.max_transport_attempts,
             guard.max_total_tokens, guard.max_completion_tokens, guard.max_rmb))
    if guard.estimated_rmb() is None:
        print('cost basis: subscription provider -> token usage is recorded, but no RMB '
              'figure is computed; the plan amount is UNKNOWN unless the platform reports it')
    total_planned = 0
    for representation in representations:
        for task_name in task_names:
            count = mode_cfg['per_task']
            if task_name == 't1_identity':
                count = mode_cfg['identity_positives'] + mode_cfg['identity_negatives']
            if args.limit:
                count = min(count, args.limit)
            frozen_items = None
            if mode == 'main':
                # the frozen list is the ONLY source of main items; it is
                # verified against the manifest before a single request is sent
                frozen_items = runner.load_frozen_items(
                    representation, task_name, frozen_manifest, run_config=run_config)[:count]
                print('  verified frozen list: %s/%s, %d items'
                      % (representation, task_name, len(frozen_items)))
            print('--- %s / %s (n=%d)' % (representation, task_name, count))
            runner.run_cell(model, representation, task_name, args.split, count,
                            run_id, run_config, manifest, dry_run=False, budget=guard,
                            pool_name=pool_name, frozen_items=frozen_items)
            total_planned += count
            if guard.stopped_reason:
                print('STOPPED: %s -- remaining cells will be recorded as pending'
                      % guard.stopped_reason)
                break
        if guard.stopped_reason:
            break
    print('budget used:', guard.snapshot())
    metrics = runner.score_run(run_id, run_config, manifest)
    summary, companion = report.write_summary_csv(run_id, metrics, representations, args.split)
    _write_run_manifest(run_id, mode, representations, task_names, model, manifest,
                        run_config, total_planned, budget=guard)
    print('\nwrote', summary)
    print('wrote', companion)
    print('wrote', paths.run_dir(run_id) / 'metrics.json')
    print('wrote', report.write_progress_matrix(
        metrics, report.run_state_map(run_id, representations, task_names)))
    return 0


def cmd_dry_run(args):
    """Build every item and prompt without contacting a model."""
    model_config, run_config = load_configs()
    manifest = dataset.load_manifest()
    representations = (args.representations.split(',') if args.representations
                       else list(reps.ROUND1_REPRESENTATIONS))
    task_names = args.tasks.split(',') if args.tasks else list(tasks.TASK_ORDER)
    run_id = args.run_id or 'dry-run'
    out_dir = paths.run_dir(run_id)
    out_dir.mkdir(parents=True, exist_ok=True)

    summary = []
    for representation in representations:
        for task_name in task_names:
            task = tasks.get(task_name)
            count = run_config['smoke']['per_task']
            if task_name == 't1_identity':
                count = run_config['smoke']['identity_positives'] + run_config['smoke']['identity_negatives']
            pool = dataset.load_prepared(args.split)
            seed = run_config['seed'] if task_name != 't1_identity' else run_config['identity_seed']
            items = task.build_items(pool, count, seed, representation, args.split)
            with open(out_dir / ('items-%s-%s.jsonl' % (representation, task_name)),
                      'w', encoding='utf-8') as handle:
                for item in items:
                    messages = task.messages(representation, item, run_config)
                    handle.write(json.dumps({
                        'example_id': item['example_id'],
                        'prompt_hash': tasks.prompt_hash(representation, item, run_config),
                        'prompt': messages[-1]['content'],
                        'label': item['label'],
                    }, ensure_ascii=False) + '\n')
            summary.append({'representation': representation, 'task': task_name,
                            'items_built': len(items), 'asked': count})
            print('%-10s %-26s items=%d asked=%d' % (representation, task_name, len(items), count))
    with open(out_dir / 'dry_run_report.json', 'w', encoding='utf-8') as handle:
        json.dump({'planned': summary, 'note': 'no model was called'}, handle,
                  ensure_ascii=False, indent=2)
    path = report.write_progress_matrix()
    print('wrote', path)
    return 0


def cmd_extension_prep(args):
    """Freeze the item lists for the five ADDED representations. Calls no model.

    Inherits the frozen main-02 molecules and the same task seeds, so an
    extension item is the SAME question about the SAME molecule as its main-02
    counterpart, with only the payload changed. Images are written to disk and
    referenced by path plus sha256, because an image is not text and must not
    be smuggled through a prompt string.
    """
    import base64
    from bmb import reps_ext
    extension = load_json(paths.CONFIG_DIR / 'extension_config.json')
    main_manifest = load_json(paths.DELIV_DIR / 'main_manifest.json')
    run_config = load_json(paths.RUN_CONFIG)
    pool = dataset.load_prepared('main_molecules')

    out_dir = paths.prepared('extension_items')
    out_dir.mkdir(parents=True, exist_ok=True)
    img_dir = out_dir / 'img'
    img_dir.mkdir(exist_ok=True)

    target = extension['targets']['per_task']
    manifest = {
        'stage': extension['stage'],
        'run_id': extension['run_id'],
        'built_at': datetime.datetime.now().isoformat(timespec='seconds'),
        'inherits': {'molecules': 'main_molecules (431)',
                     'main_molecule_list_sha256': main_manifest['main_molecule_list_sha256'],
                     'dataset_version': main_manifest['dataset_version']},
        'conversion': reps_ext.CONVERSION,
        'identity_views': {
            '_rule': 'slot A uses view 0 and slot B uses view 1 for a positive pair AND for '
                     'a negative pair, so the view a molecule is drawn in carries no '
                     'information about the answer',
            'pair_views': list(EXTENSION_PAIR_VIEWS),
            '_not_applied_to': 'the three inherited representations reuse their frozen '
                               'main-02 responses, which were produced under the earlier '
                               'construction (positives 0/1, negatives 0/0); that is '
                               'recorded as a stimulus-design difference, not hidden',
        },
        'representations': {}, 'item_lists': {}, 'prompt_hashes': {}, 'coverage': {},
    }
    for representation in extension['representations_added']:
        manifest['representations'][representation] = {
            'kind': reps_ext.REPRESENTATION_KIND[representation],
            'display_name': reps_ext.PILOT_DISPLAY_NAME[representation],
        }
        converted = 0
        failed = {}
        for task_name in tasks.TASK_ORDER:
            count = target
            if task_name == 't1_identity':
                count = (extension['targets']['identity_positives']
                         + extension['targets']['identity_negatives'])
            seed = (run_config['identity_seed'] if task_name == 't1_identity'
                    else run_config['seed'])
            extra = {'pair_views': EXTENSION_PAIR_VIEWS} if task_name == 't1_identity' else {}
            try:
                items = tasks.get(task_name).build_items(pool, count, seed,
                                                         representation, 'test', **extra)
            except Exception as exc:  # noqa: BLE001
                failed[task_name] = '%s: %s' % (type(exc).__name__, str(exc)[:120])
                continue
            if task_name == 't1_identity':
                positives = sum(1 for i in items if i['label']['same_molecule'])
                negatives = len(items) - positives
                if positives == 0 or negatives == 0:
                    failed[task_name] = ('degenerate Identity set: %d positive / %d negative '
                                         '-- refused rather than run' % (positives, negatives))
                    continue
            serialisable = []
            for item in items:
                row = {k: v for k, v in item.items() if k != 'payload'}
                row['payload'] = {}
                for key, value in item['payload'].items():
                    if isinstance(value, bytes):
                        name = '%s-%s.png' % (item['example_id'], key)
                        (img_dir / name).write_bytes(value)
                        row['payload'][key] = {
                            'image_path': str((img_dir / name).relative_to(paths.BENCH_DIR)),
                            'sha256': hashlib.sha256(value).hexdigest(),
                            'bytes': len(value)}
                    else:
                        row['payload'][key] = value
                serialisable.append(row)
            path = out_dir / ('%s-%s.jsonl' % (representation, task_name))
            with open(path, 'w', encoding='utf-8') as handle:
                for row in serialisable:
                    handle.write(json.dumps(row, ensure_ascii=False) + '\n')
            blob = '\n'.join(json.dumps(r, ensure_ascii=False, sort_keys=True)
                             for r in serialisable)
            digest = hashlib.sha256(blob.encode('utf-8')).hexdigest()
            manifest['item_lists'].setdefault(task_name, {})[representation] = {
                'n_items': len(serialisable), 'sha256': digest,
                'pool': 'extension_items/%s-%s.jsonl' % (representation, task_name)}
            manifest['prompt_hashes'].setdefault(representation, {})[task_name] = [
                tasks.prompt_hash(representation, item, run_config) for item in items]
            converted += len(serialisable)
            print('  %-9s %-26s %3d items' % (representation, task_name, len(serialisable)))
        manifest['coverage'][representation] = {
            'items_built': converted, 'target': target * len(tasks.TASK_ORDER),
            'failed_tasks': failed}
    path = paths.DELIV_DIR / 'extension_manifest.json'
    path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
    print()
    print('wrote', path)
    return 0


def cmd_extension_checks(args):
    """Offline checks for the five added representations. No network, no tokens.

    These cannot show that a model reads a picture correctly -- only a real call
    can do that. What they CAN show, and what several earlier claims skipped, is
    that the request the client builds really carries an image content block
    rather than a path or a base64 string in the text, that the PNG carries no
    hidden metadata, that the text formats round-trip, and that the frozen items
    are the ones that will be sent.
    """
    _, run_config = load_configs()
    failures, total = [], []

    def check(name, got, want):
        total.append(name)
        ok = got == want
        print('  %-62s %s' % (name, 'ok' if ok else 'FAIL got=%r want=%r' % (got, want)))
        if not ok:
            failures.append(name)

    from bmb import reps_ext

    sample = 'CC(=O)Oc1ccccc1C(=O)O'

    # 1. the image really travels as an image_url content block
    blob = reps_ext.img_png(sample)
    messages = tasks.get('t3_atom_counts').messages(
        'img', {'payload': {'molecule': blob}, 'example_id': 'x', 'task': 't3_atom_counts',
                'rep': 'img', 'cluster': 'x', 'label': {}}, run_config)
    content = messages[-1]['content']
    check('an image request is a block list, not a string', isinstance(content, list), True)
    kinds = [b.get('type') for b in content] if isinstance(content, list) else []
    check('the image request has one text block', kinds.count('text'), 1)
    check('the image request has one image_url block', kinds.count('image_url'), 1)
    url = [b for b in content if b.get('type') == 'image_url'][0]['image_url']['url']
    check('the image is sent as a PNG data URI', url.startswith('data:image/png;base64,'), True)
    check('no file path leaks into the request text',
          'Design Molecule' in json.dumps(content, ensure_ascii=False), False)

    # 2. Identity sends TWO images, one per endpoint
    pair = tasks.get('t1_identity').messages(
        'img', {'payload': {'a': blob, 'b': reps_ext.img_png(sample, variant=1)},
                'example_id': 'x', 'task': 't1_identity', 'rep': 'img',
                'cluster': 'x', 'label': {'same_molecule': True}}, run_config)
    kinds = [b.get('type') for b in pair[-1]['content']]
    check('an Identity pair sends two image blocks', kinds.count('image_url'), 2)
    check('an Identity image pair is marked A and B', 'A:' in pair[-1]['content'][0]['text']
          and 'B:' in pair[-1]['content'][0]['text'], True)

    # 3. the PNG carries pixels only
    check('the PNG has no text/iTXt/zTXt metadata chunk',
          reps_ext.png_has_text_metadata(blob), False)
    check('the two Identity views are different files',
          blob != reps_ext.img_png(sample, variant=1), True)

    # 4. no text representation inlines an image or a path
    for name in ('graph2d', 'graph3d', 'cml', 'moljson'):
        text = reps_ext.render(name, sample)
        check('%s sends a plain text payload' % name, isinstance(text, str), True)
    check('Graph2D is a MOL V2000 block', reps_ext.graph2d_text(sample).splitlines()[3]
          .strip().endswith('V2000'), True)
    check('Graph3D marks the table as 3D',
          '3D' in reps_ext.graph3d_text(sample).splitlines()[1], True)
    check('Graph2D marks the table as 2D',
          '2D' in reps_ext.graph2d_text(sample).splitlines()[1], True)
    check('Graph2D z coordinates are zero',
          all(line.split()[2] == '0.0000'
              for line in reps_ext.graph2d_text(sample).splitlines()[4:4 + 5]), True)
    check('Graph3D z coordinates are not all zero',
          any(float(line.split()[2]) != 0.0
              for line in reps_ext.graph3d_text(sample).splitlines()[4:4 + 5]), True)
    check('CML parses as XML',
          reps_ext.cml_text(sample).startswith('<?xml'), True)
    keys = set(json.loads(reps_ext.moljson_text(sample)))
    check('MolJSON carries the schema-required keys',
          sorted(set(reps_ext.moljson_schema_keys()['required']) - keys), [])
    check('MolJSON ids do not encode element counts',
          any(a['id'].startswith('C') for a in
              json.loads(reps_ext.moljson_text(sample))['atoms']), False)

    # 5. no payload leaks a name, a database id or the answer
    for name in reps_ext.EXTENSION_NAMES:
        payload = reps_ext.render(name, sample)
        blob_text = payload if isinstance(payload, str) else payload.decode('latin-1')
        check('%s carries no example id or answer' % name,
              ('example_id' in blob_text) or ('same_molecule' in blob_text)
              or ('qed' in blob_text.lower()), False)

    print('\n%d offline extension checks, %d failed' % (len(total), len(failures)))
    return 1 if failures else 0


def cmd_extension_dev(args):
    """The dev-side validation of the added representations. Needs a credential.

    Identity is not exercised here: what is being validated is whether the model
    can read a NEW FORMAT at all, and the cheapest way to see that is a question
    whose answer changes with the molecule. For img the same question is asked
    about two different pictures, because HTTP 200 does not show that the image
    was read.
    """
    model_config, run_config = load_configs()
    extension = load_json(paths.CONFIG_DIR / 'extension_config.json')
    resolved = _require_model(model_config)
    if resolved is None:
        return 2
    model = client_mod.client_from_config(model_config, run_config)
    cap = extension['dev_check']['per_representation_max_calls']
    run_id = args.run_id or 'extension-dev'
    out_dir = paths.run_dir(run_id)
    out_dir.mkdir(parents=True, exist_ok=True)

    from bmb import reps_ext
    pool = dataset.load_prepared('dev')
    pool = sorted(pool, key=lambda r: r['cluster'])[:4]
    reps_to_check = (args.representations.split(',') if args.representations
                     else ['img', 'graph2d', 'graph3d', 'cml', 'moljson'])
    task = tasks.get('t3_atom_counts')

    log, calls = [], 0
    for representation in reps_to_check:
        for index, entry in enumerate(pool[:2]):
            if calls >= cap * len(reps_to_check):
                break
            try:
                payload = reps_ext.render(representation, entry['smiles'])
            except Exception as exc:  # noqa: BLE001
                log.append({'representation': representation, 'example_id': entry['cid'],
                            'status': 'conversion_failed', 'error': str(exc)[:200]})
                continue
            item = {'example_id': 'dev-%s-%s' % (representation, entry['cid']),
                    'task': 't3_atom_counts', 'rep': representation, 'split': 'dev',
                    'cluster': entry['cluster'], 'payload': {'molecule': payload},
                    'label': labels.atom_counts(entry['smiles'])}
            messages = task.messages(representation, item, run_config)
            result = model.chat(messages)
            calls += 1
            row = {'representation': representation, 'example_id': item['example_id'],
                   'smiles': entry['smiles'], 'label': item['label'],
                   'status': 'ok' if result['ok'] else 'transport_failed',
                   'model_id': model.model_id,
                   'finish_reason': result.get('finish_reason')}
            if result['ok']:
                parsed, error = task.parse(result['text'])
                row['raw_response'] = result['text']
                row['parsed'] = parsed
                row['parse_error'] = error
                row['exact_match'] = (parsed == item['label']) if parsed else None
                row['image_sent'] = any(b.get('type') == 'image_url'
                                        for b in messages[-1]['content']
                                        ) if isinstance(messages[-1]['content'], list) else False
            else:
                row['transport_error'] = result.get('transport_error')
            log.append(row)
            print('  %-9s %-22s %s' % (representation, item['example_id'], row['status']))
    # the same question about TWO DIFFERENT pictures: a text-only channel would
    # give the same answer twice, or answer without the picture at all
    img_rows = [r for r in log if r['representation'] == 'img' and r.get('parsed')]
    if len(img_rows) >= 2:
        differs = img_rows[0]['parsed'] != img_rows[1]['parsed']
        print('\nimg read two different pictures and %s'
              % ('gave different counts' if differs else 'gave the SAME counts'))
    path = out_dir / 'dev_check.json'
    path.write_text(json.dumps({'run_id': run_id, 'calls': calls, 'cap_per_representation': cap,
                                'rows': log,
                                'note': 'dev only; these responses are NOT scored results'},
                               ensure_ascii=False, indent=2), encoding='utf-8')
    _write_run_manifest(run_id, 'extension-dev', reps_to_check, ['t3_atom_counts'], model,
                        dataset.load_manifest(), run_config, calls)
    print('\nwrote', path)
    return 0


def cmd_extension_run(args):
    """Execute the five ADDED representations against the frozen extension items.

    Same rules as `main`: the frozen list is the only source of items, the model,
    provider, endpoint and generation settings are inherited, only transport
    failures retry, and every response is persisted as it arrives. The stage has
    its own budget block, so this run neither resets nor consumes the previous
    cycle's ledger.
    """
    model_config, run_config = load_configs()
    extension = load_json(paths.CONFIG_DIR / 'extension_config.json')
    manifest = dataset.load_manifest()
    ext_manifest_path = paths.DELIV_DIR / 'extension_manifest.json'
    if not ext_manifest_path.is_file():
        print('This run cannot start: the frozen extension item lists are missing. '
              'Run `extension-prep` first.')
        return 2
    ext_manifest = load_json(ext_manifest_path)

    representations = (args.representations.split(',') if args.representations
                       else list(extension['representations_added']))
    task_names = args.tasks.split(',') if args.tasks else list(tasks.TASK_ORDER)
    run_id = args.run_id or extension['run_id']

    resolved = _require_model(model_config)
    if resolved is None:
        return 2
    model = client_mod.client_from_config(model_config, run_config)
    print('run_id=%s stage=extension provider=%s model=%s region=%s'
          % (run_id, model.provider, model.model_id, model.region))
    print('base_url=%s' % model.base_url)
    print('credential variable in use: %s' % (model.api_key_var_in_use or 'recognised'))
    warning = client_mod.provider_terms_warning(model.provider)
    if warning:
        print()
        print('!' * 78)
        print('NOTICE -- this project does not claim the selected provider is permitted')
        print('for this use. The platform states, verbatim:')
        print('  %s' % client_mod.TOKEN_PLAN_TERMS_CN)
        print('  source: %s' % client_mod.TOKEN_PLAN_TERMS_SOURCE)
        print('The choice is the user\'s and is recorded in the run manifest.')
        print('!' * 78)
        print()

    guard = runner.budget_from_extension(extension, model.provider)
    print('budget guard: basic<=%s attempts<=%s tokens<=%s completion<=%s'
          % (guard.max_basic_requests, guard.max_transport_attempts,
             guard.max_total_tokens, guard.max_completion_tokens))

    # cell order is shuffled with a recorded seed so no added representation is
    # always first or always last; the ITEMS inside a cell keep their frozen order
    cells = [(r, t) for r in representations for t in task_names]
    random.Random(args.order_seed).shuffle(cells)
    print('cell order seed=%d' % args.order_seed)

    total_planned = 0
    for representation, task_name in cells:
        if guard.stopped_reason:
            print('STOPPED: %s -- remaining cells are not started' % guard.stopped_reason)
            break
        try:
            items = runner.load_extension_items(representation, task_name, ext_manifest,
                                                run_config=run_config)
        except runner.FrozenListError as exc:
            print('  SKIP %s/%s: %s' % (representation, task_name, exc))
            continue
        if args.limit:
            items = items[:args.limit]
        print('--- %s / %s (n=%d of frozen %d)'
              % (representation, task_name, len(items),
                 ext_manifest['item_lists'][task_name][representation]['n_items']))
        runner.run_cell(model, representation, task_name, 'test', len(items), run_id,
                        run_config, manifest, dry_run=False, budget=guard,
                        pool_name='main_molecules', frozen_items=items)
        total_planned += len(items)

    print('budget used:', guard.snapshot())
    metrics = runner.score_run(run_id, run_config, manifest)
    summary, companion = report.write_summary_csv(run_id, metrics, representations, 'test')
    path = _write_run_manifest(run_id, 'extension', representations, task_names, model,
                               manifest, run_config, total_planned, budget=guard)
    print('\nwrote', summary)
    print('wrote', companion)
    print('wrote', paths.run_dir(run_id) / 'metrics.json')
    print('wrote', path)
    return 0


def cmd_main_prep(args):
    """Build and freeze the main-scale sample. Calls no model.

    Applies the complexity gate to the PRISTINE test split, checks that every
    round-1 representation can still express the survivors, measures how many
    items each task can actually support, and writes a frozen manifest plus the
    frozen item list with its hash.
    """
    run_config = load_json(paths.RUN_CONFIG)
    manifest = dataset.load_manifest()
    target = run_config['main']['per_task']
    test_rows = dataset.load_prepared('test')
    dev_rows = dataset.load_prepared('dev')

    kept, dropped = dataset.apply_complexity_filter(test_rows)
    kept_dev, dropped_dev = dataset.apply_complexity_filter(dev_rows)
    digest = hashlib.sha256('\n'.join(r['smiles'] for r in kept).encode('utf-8')).hexdigest()

    # every surviving molecule must still be expressible in all three arms
    expressible = 0
    failed = {}
    for row in kept:
        try:
            for rep in reps.ROUND1_REPRESENTATIONS:
                reps.render(rep, row['smiles'])
            expressible += 1
        except Exception as exc:  # noqa: BLE001
            failed[row['cid']] = '%s: %s' % (type(exc).__name__, exc)

    # functional-group support over the filtered evaluation set
    fg = {}
    for name in run_benchmark_labels():
        pos = neg = 0
        for row in kept:
            if labels.functional_group_labels(row['smiles'], only_enabled=True)[name]:
                pos += 1
            else:
                neg += 1
        fg[name] = {'positive': pos, 'negative': neg}

    # how many items can each task actually build at the target size?
    feasibility = {}
    pair_usage = {}
    for task_name in tasks.TASK_ORDER:
        task = tasks.get(task_name)
        counts = {}
        for rep in reps.ROUND1_REPRESENTATIONS:
            count = target
            if task_name == 't1_identity':
                count = run_config['main']['identity_positives'] + run_config['main']['identity_negatives']
            seed = (run_config['identity_seed'] if task_name == 't1_identity'
                    else run_config['seed'])
            items = task.build_items(kept, count, seed, rep, 'test')
            counts[rep] = len(items)
            if task_name == 't1_identity' and rep == 'smiles':
                used = Counter()
                for item in items:
                    for unit in (item.get('clusters') or (item['cluster'],)):
                        used[unit] += 1
                pair_usage = {
                    'n_pairs': len(items),
                    'n_distinct_molecules_used': len(used),
                    'max_pairs_sharing_one_molecule': max(used.values()) if used else 0,
                    'molecules_used_in_more_than_one_pair':
                        sum(1 for v in used.values() if v > 1),
                    'molecules_used_in_more_than_two_pairs':
                        sum(1 for v in used.values() if v > 2),
                }
        feasibility[task_name] = {
            'target': target,
            'built': counts,
            'shortfall': {r: target - c for r, c in counts.items() if c < target},
            'meets_target': all(c >= target for c in counts.values()),
        }

    report = {
        'generated_at': datetime.datetime.now().isoformat(timespec='seconds'),
        'dataset_version': manifest['dataset_version'],
        'complexity_gate': dataset.COMPLEXITY,
        'frozen_items': {},
        'counts': {
            'test_split_before': len(test_rows),
            'test_split_after': len(kept),
            'dropped': len(dropped),
            'drop_reasons': dict(Counter(r['drop_reason'] for r in dropped)),
            'dev_split_before': len(dev_rows),
            'dev_split_after': len(kept_dev),
            'dev_dropped': len(dropped_dev),
        },
        'complexity_before': dataset.complexity_distribution(test_rows),
        'complexity_after': dataset.complexity_distribution(kept),
        'expressible_in_all_three_representations': expressible,
        'expression_failures': failed,
        'functional_group_support': fg,
        'functional_group_labels_used': sorted(fg),
        'identity_pair_usage': pair_usage,
        'task_feasibility': feasibility,
        'note': ('The dev split is filtered with the same gate for any later dev-side '
                 'adaptation check, but dev is NOT part of the main evaluation set.'),
    }

    # ---- freeze the actual item list, and measure the realised label support
    frozen_dir = paths.prepared('main_items')
    frozen_dir.mkdir(parents=True, exist_ok=True)
    rule = run_config['main'].get('label_support_rule', {})
    min_pos = rule.get('min_positives_in_sample', 10)
    min_neg = rule.get('min_negatives_in_sample', 10)
    realised = {}
    for task_name in tasks.TASK_ORDER:
        task = tasks.get(task_name)
        count = target
        if task_name == 't1_identity':
            count = run_config['main']['identity_positives'] + run_config['main']['identity_negatives']
        seed = (run_config['identity_seed'] if task_name == 't1_identity'
                else run_config['seed'])
        per_rep = {}
        for rep in reps.ROUND1_REPRESENTATIONS:
            items = task.build_items(kept, count, seed, rep, 'test')
            path = frozen_dir / ('%s-%s.jsonl' % (rep, task_name))
            with open(path, 'w', encoding='utf-8') as handle:
                for item in items:
                    handle.write(json.dumps(item, ensure_ascii=False) + '\n')
            blob = '\n'.join(json.dumps(i, ensure_ascii=False, sort_keys=True) for i in items)
            per_rep[rep] = {'n_items': len(items),
                            'sha256': hashlib.sha256(blob.encode('utf-8')).hexdigest()}
            if task_name == 't2_functional_groups':
                tally = {}
                for item in items:
                    for name, present in item['label'].items():
                        bucket = tally.setdefault(name, {'positive': 0, 'negative': 0})
                        bucket['positive' if present else 'negative'] += 1
                realised = tally
        report['frozen_items'][task_name] = per_rep

    report['realised_label_support_in_frozen_sample'] = realised
    report['label_support_rule'] = {
        'min_positives_in_sample': min_pos,
        'min_negatives_in_sample': min_neg,
        'declared': ('fixed BEFORE any main model call, on the frozen item list only; a label '
                     'outside this band is excluded from the main macro-F1 and reported as '
                     'under-supported, never scored as 1.0'),
        'per_label': {},
    }
    included = []
    for name, tally in sorted(realised.items()):
        ok = tally['positive'] >= min_pos and tally['negative'] >= min_neg
        report['label_support_rule']['per_label'][name] = dict(
            tally, included_in_main_macro_f1=ok,
            reason='' if ok else ('too few positives in the frozen sample'
                                  if tally['positive'] < min_pos
                                  else 'too few negatives in the frozen sample'))
        if ok:
            included.append(name)
    report['labels_in_main_macro_f1'] = sorted(included)
    report['n_labels_in_main_macro_f1'] = len(included)
    report['n_labels_defined'] = len(realised)
    # record the NAMES, not just the threshold: scoring must average over this
    # exact set rather than re-deriving 'exercised' labels from the sample
    report['label_support_rule']['label_set'] = sorted(included)

    # ---- Identity endpoint containment and separation --------------------
    # Each pair owns TWO molecules. Both must be inside the frozen main set,
    # inside the complexity gate, and absent from dev and from the quarantine.
    dev_keys = {row['cluster'] for row in dev_rows}
    quarantined_path = paths.prepared('quarantined.jsonl')
    quarantine_keys = ({row['cluster'] for row in dataset.load_prepared('quarantined')}
                       if quarantined_path.is_file() else set())
    main_keys = {row['cluster'] for row in kept}
    endpoints = set()
    offenders = []
    for item in tasks.get('t1_identity').build_items(
            kept, target, run_config['identity_seed'], 'smiles', 'test'):
        owned = item.get('clusters') or (item['cluster'],)
        endpoints.update(owned)
        for unit in owned:
            if unit not in main_keys:
                offenders.append({'example_id': item['example_id'], 'reason': 'not in frozen main set'})
            elif unit in dev_keys:
                offenders.append({'example_id': item['example_id'], 'reason': 'also in dev'})
            elif unit in quarantine_keys:
                offenders.append({'example_id': item['example_id'], 'reason': 'in quarantine'})
    report['identity_endpoint_check'] = {
        'n_endpoints': len(endpoints),
        'all_endpoints_inside_frozen_main_set': not offenders,
        'any_endpoint_in_dev': bool([o for o in offenders if o['reason'] == 'also in dev']),
        'any_endpoint_in_quarantine': bool([o for o in offenders if 'quarantine' in o['reason']]),
        'violations': offenders[:20],
        'rule': ('both endpoints of every Identity pair come from the same frozen main set; '
                 'no endpoint may lie in dev or in the debug quarantine'),
        'gate_note': ('the main set already passed the complexity gate, so containment implies '
                      'both endpoints satisfy MW >= 300 Da and rings >= 2'),
    }

    # ---- dev-side offline adaptation check -------------------------------
    # dev is filtered by the same gate and must be able to build the same item
    # shapes, so any later dev-side adaptation has the same structure as main.
    # dev is NOT part of the evaluation set.
    dev_feasibility = {}
    for task_name in tasks.TASK_ORDER:
        task = tasks.get(task_name)
        count = target
        if task_name == 't1_identity':
            count = run_config['main']['identity_positives'] + run_config['main']['identity_negatives']
        seed = (run_config['identity_seed'] if task_name == 't1_identity'
                else run_config['seed'])
        built = len(task.build_items(kept_dev, count, seed, 'smiles', 'dev'))
        dev_feasibility[task_name] = {'target': count, 'built_on_dev': built,
                                      'meets_target': built >= count}
    report['dev_adaptation_check'] = dev_feasibility

    # ---- freeze ----------------------------------------------------------
    model_config = load_json(paths.MODEL_CONFIG)
    prompt_hashes = {}
    for task_name in tasks.TASK_ORDER:
        task = tasks.get(task_name)
        count = target
        if task_name == 't1_identity':
            count = run_config['main']['identity_positives'] + run_config['main']['identity_negatives']
        seed = (run_config['identity_seed'] if task_name == 't1_identity'
                else run_config['seed'])
        for rep in reps.ROUND1_REPRESENTATIONS:
            items = task.build_items(kept, count, seed, rep, 'test')
            prompt_hashes.setdefault(rep, {})[task_name] = sorted(
                tasks.prompt_hash(rep, item, run_config) for item in items)
    report['prompt_hashes'] = prompt_hashes
    report['frozen'] = {
        'frozen_at': datetime.datetime.now().isoformat(timespec='seconds'),
        'dataset_version': manifest['dataset_version'],
        'main_molecule_list_sha256': digest,
        'n_main_molecules': len(kept),
        'item_list_hashes': {task: {rep: v['sha256'] for rep, v in per.items()}
                             for task, per in report['frozen_items'].items()},
        # per-item prompt hashes, so a run can verify that what it is about to
        # send is what was frozen rather than trusting the list by construction
        'prompt_hashes': prompt_hashes,
        'provider': (model_config.get('active') or {}).get('provider', 'bailian'),
        'provider_selected_by': (model_config.get('active') or {}).get('provider_selected_by'),
        'provider_terms_warning': client_mod.provider_terms_warning(
            (model_config.get('active') or {}).get('provider', 'bailian')),
        'model_id': (model_config.get('active') or {}).get('model_id'),
        'region': (model_config.get('active') or {}).get('region'),
        'base_url': (model_config.get('active') or {}).get('base_url'),
        'credential_variable': client_mod.PROVIDERS.get(
            (model_config.get('active') or {}).get('provider', 'bailian'),
            {}).get('credential_envs'),
        'generation': {k: v for k, v in run_config['generation'].items()
                       if not k.startswith('_')},
        'scoring': {
            'qed_success_threshold': run_config['qed_success_threshold'],
            'bootstrap': run_config['bootstrap'],
            'interval_methods': ['wilson', 'clopper-pearson', 'percentile cluster bootstrap'],
            'paired_test': 'cluster sign-flip permutation',
            'multiple_comparison': scoring.DECLARED_COMPARISON_FAMILY,
            'label_rule': run_config['main']['label_support_rule'],
        },
        'seed': run_config['seed'],
        'identity_seed': run_config['identity_seed'],
        'budget': run_config['budget'],
        'resume_and_cache': run_config['resume_and_cache'],
        'statement': ('frozen before any main model call. The formal test list is never used '
                      'to tune prompts, thresholds or scoring rules.'),
    }
    freeze_path = paths.DELIV_DIR / 'main_manifest.json'
    freeze_path.write_text(json.dumps(report['frozen'], ensure_ascii=False, indent=2),
                           encoding='utf-8')

    paths.DELIV_DIR.mkdir(parents=True, exist_ok=True)
    out = paths.DELIV_DIR / 'main_feasibility.json'
    out.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')

    frozen = paths.prepared('main_molecules.jsonl')
    with open(frozen, 'w', encoding='utf-8') as handle:
        for row in kept:
            handle.write(json.dumps(row, ensure_ascii=False) + '\n')
    (paths.prepared('main_molecules.canonical_smiles.sha256')).write_text(digest + '\n',
                                                                         encoding='utf-8')

    print('test split %d -> %d after the complexity gate (dropped %d)'
          % (len(test_rows), len(kept), len(dropped)))
    print('expressible in all three representations: %d/%d' % (expressible, len(kept)))
    print('frozen main molecule list hash:', digest[:32])
    print()
    print('%-26s %6s %8s %8s' % ('task', 'target', 'built', 'ok'))
    for name, entry in feasibility.items():
        print('%-26s %6d %8d %8s' % (name, entry['target'],
                                     min(entry['built'].values()),
                                     'yes' if entry['meets_target'] else 'NO'))
    print()
    print('%-16s %8s %8s' % ('functional group', 'positive', 'negative'))
    for name, entry in sorted(fg.items()):
        print('%-16s %8d %8d' % (name, entry['positive'], entry['negative']))
    print()
    print('realised support in the FROZEN 200-item T2 sample '
          '(rule: >=%d positives and >=%d negatives)'
          % (min_pos, min_neg))
    print('%-16s %8s %8s   %s' % ('label', 'positive', 'negative', 'in main macro-F1'))
    for name, entry in sorted(report['label_support_rule']['per_label'].items()):
        print('%-16s %8d %8d   %s' % (name, entry['positive'], entry['negative'],
                                      'yes' if entry['included_in_main_macro_f1']
                                      else 'NO -- ' + entry['reason']))
    print()
    print('labels entering the main macro-F1: %d of %d'
          % (report['n_labels_in_main_macro_f1'], report['n_labels_defined']))
    print('wrote', out)
    return 0


def run_benchmark_labels():
    return [item['name'] for item in labels.enabled_labels()]


def cmd_compare(args):
    """The pre-declared paired comparisons, with Holm-Bonferroni correction.

    Run only at main scale. Two arms are compared on the SAME items (paired),
    by a cluster sign-flip permutation, and the three pairs within a task form
    one declared family. No comparison is run on smoke or on the connection
    check: those are pipeline checks, not results.
    """
    _, run_config = load_configs()
    manifest = dataset.load_manifest()
    resamples = run_config.get('bootstrap', {}).get('resamples', 2000)
    threshold = run_config.get('qed_success_threshold', 0.05)
    label_names = [item['name'] for item in labels.enabled_labels()]
    records = runner.load_run_records(args.run)
    if not records:
        print('no records in run %r' % args.run)
        return 2

    splits = sorted({r.get('split') for r in records})
    out = {'run_id': args.run, 'family': scoring.DECLARED_COMPARISON_FAMILY,
           'per_task': {}}
    for split_name in splits:
        for task_name in tasks.TASK_ORDER:
            cells = {}
            for rep in reps.ROUND1_REPRESENTATIONS:
                rows = [r for r in records
                        if r['representation'] == rep and r['task'] == task_name
                        and r.get('split') == split_name]
                if rows:
                    cells[rep] = rows
            if len(cells) < 2:
                continue
            raw_p = {}
            details = {}
            metric_details = {}
            for a, b in scoring.DECLARED_COMPARISON_FAMILY['comparisons']:
                if a not in cells or b not in cells:
                    continue
                key = '%s_vs_%s' % (a, b)
                # PRIMARY: the difference in the metric the main table prints,
                # resampled in pairs so both arms see the same draw
                metric_details[key] = scoring.paired_metric_bootstrap(
                    cells[a], cells[b], task_name, label_names, threshold,
                    resamples=resamples, seed=run_config['seed'])
                # SECONDARY: a per-item permutation test. Its per-item score is
                # NOT the main-table metric for t2, so its p-value is reported
                # separately and must not be read as a Macro-F1 result.
                result = scoring.paired_comparison(cells[a], cells[b], task_name,
                                                   resamples=resamples,
                                                   seed=run_config['seed'])
                raw_p[key] = result.get('p_value')
                details[key] = result
            adjusted = scoring.holm_adjust(raw_p)
            out['per_task']['%s|%s' % (task_name, split_name)] = {
                'metric_comparisons': metric_details,
                'per_item_permutation': details, 'p_raw': raw_p, 'p_holm': adjusted}

    path = paths.run_dir(args.run) / 'comparison.json'
    path.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding='utf-8')
    print('wrote', path)
    print()
    print('PRIMARY: paired difference in the MAIN-TABLE metric '
          '(cluster resampled in pairs; CI excluding 0 = real difference)')
    for cell, block in out['per_task'].items():
        print(cell)
        for key, detail in block['metric_comparisons'].items():
            ci = detail.get('ci') or {}
            print('  %-22s n=%-4s %s=%.3f  diff=%+.3f  CI=[%.3f, %.3f]  excludes0=%s'
                  % (key, detail.get('n_pairs'), detail.get('metric'),
                     detail.get('observed_difference') or 0.0,
                     detail.get('observed_difference') or 0.0,
                     ci.get('low', 0), ci.get('high', 0), ci.get('excludes_zero')))
    print()
    print('SECONDARY: per-item permutation test (different per-item score; '
          'for t2 this is NOT Macro-F1)')
    for cell, block in out['per_task'].items():
        for key, detail in block['per_item_permutation'].items():
            print('  %-24s %-20s diff=%+.3f p_raw=%-8s p_holm=%s'
                  % (cell.split('|')[0], key, detail.get('observed_mean_difference') or 0.0,
                     ('%.4f' % block['p_raw'][key]) if block['p_raw'].get(key) is not None else 'n/a',
                     ('%.4f' % block['p_holm'][key]) if block['p_holm'].get(key) is not None else 'n/a'))
    return 0


def cmd_ext_compare(args):
    """The pre-declared nine-pair family, per task, across all eight rows.

    The five added representations come from the extension run; the three
    inherited ones are reused from main-02 ONLY when their model, provider,
    endpoint and generation match what the extension run actually used. If they
    do not match, the old rows are reported as a separate historical table and
    are NOT placed in the same ranking -- a silent model swap would make the
    whole comparison meaningless.
    """
    model_config, run_config = load_configs()
    extension = load_json(paths.CONFIG_DIR / 'extension_config.json')
    resamples = run_config.get('bootstrap', {}).get('resamples', 2000)
    threshold = run_config.get('qed_success_threshold', 0.05)
    label_names = labels.frozen_macro_f1_labels() or [
        item['name'] for item in labels.enabled_labels()]

    added = list(extension['representations_added'])
    inherited = list(extension['representations_inherited'])

    new_records = runner.load_run_records(args.run)
    if not new_records:
        print('no records in extension run %r' % args.run)
        return 2
    old_records = runner.load_run_records(args.reuse_run) if args.reuse_run else []

    new_manifest_path = paths.run_dir(args.run) / 'experiment_manifest.json'
    new_model = {}
    if new_manifest_path.is_file():
        new_model = (load_json(new_manifest_path).get('model') or {})

    def condition_of(records, manifest_model=None):
        for record in records:
            row = {k: record.get(k) for k in ('model_id', 'provider', 'base_url')}
            if manifest_model:
                row['generation'] = manifest_model.get('generation')
            return row
        return None

    reuse_ok, reuse_note = False, 'no reuse run was named'
    old_condition = condition_of(old_records)
    new_condition = condition_of(new_records)
    if old_records:
        old_manifest_path = paths.run_dir(args.reuse_run) / 'experiment_manifest.json'
        old_manifest_model = {}
        if old_manifest_path.is_file():
            old_manifest_model = (load_json(old_manifest_path).get('model') or {})
        old_condition = condition_of(old_records, old_manifest_model)
        new_condition = condition_of(new_records, new_model)
        comparable = ('model_id', 'provider', 'base_url')
        mismatched = [k for k in comparable if old_condition.get(k) != new_condition.get(k)]
        if mismatched:
            reuse_note = ('the reused run was produced under a different %s, so the '
                          'inherited rows are NOT pooled with the added ones'
                          % ', '.join(mismatched))
        else:
            reuse_ok = True
            reuse_note = ('model, provider, endpoint and generation match between %s and %s; '
                          'the inherited responses carry the same prompt for the same item'
                          % (args.reuse_run, args.run))

    pooled = list(new_records) + (list(old_records) if reuse_ok else [])
    representations = added + (inherited if reuse_ok else [])
    family = ([tuple(p) for p in scoring.EXTENSION_COMPARISON_FAMILY['within_each_task']]
              + [tuple(p) for p in scoring.EXTENSION_COMPARISON_FAMILY['extra_pairs']])
    family_size = scoring.EXTENSION_COMPARISON_FAMILY['size_per_task']

    out = {
        'extension_run': args.run, 'reuse_run': args.reuse_run,
        'reuse_verified': reuse_ok, 'reuse_note': reuse_note,
        'condition_added': new_condition, 'condition_inherited': old_condition,
        'representations_in_this_table': representations,
        'representations_reported_separately': [] if reuse_ok else inherited,
        'family': {'pairs': [list(p) for p in family], 'size_per_task': family_size,
                   'correction': scoring.EXTENSION_COMPARISON_FAMILY['correction'],
                   '_untested_handling': ('an untested pair keeps its place in the family; '
                                          'the divisor is not reduced to the pairs that '
                                          'happened to be computable')},
        'inference': extension['inference_rules'],
        'per_task': {},
    }

    for task_name in tasks.TASK_ORDER:
        cells = {}
        for rep in representations:
            rows = [r for r in pooled
                    if r['representation'] == rep and r['task'] == task_name
                    and r.get('split') == 'test']
            if rows:
                cells[rep] = rows
        if len(cells) < 2:
            continue
        details, raw_p = {}, {}
        for a, b in family:
            key = '%s_vs_%s' % (a, b)
            if a not in cells or b not in cells:
                raw_p[key] = None
                details[key] = {'status': 'not tested',
                                'why': 'one or both arms have no valid responses for this task'}
                continue
            # EVIDENCE 1: the paired bootstrap interval for the main-table
            # metric difference. Reported UNCORRECTED, and labelled as such.
            detail = scoring.paired_metric_bootstrap(
                cells[a], cells[b], task_name, label_names, threshold,
                resamples=resamples, seed=run_config['seed'])
            detail['interval_note'] = ('95% percentile bootstrap over clusters, UNCORRECTED; '
                                       'Holm correction applies to the p-values, not to this '
                                       'interval')
            # EVIDENCE 2: a permutation p-value for the SAME metric, with whole
            # item prediction vectors swapped, so the family can be corrected
            perm = scoring.paired_metric_permutation(
                cells[a], cells[b], task_name, label_names, threshold,
                resamples=resamples, seed=run_config['seed'])
            detail['permutation'] = perm
            detail['status'] = 'tested'
            details[key] = detail
            raw_p[key] = perm.get('p_value')
        out['per_task'][task_name] = {
            'n_cells': len(cells),
            'coverage': {rep: {'n_records': len(rows),
                               'n_ok': sum(1 for r in rows if r.get('status') == 'ok')}
                         for rep, rows in sorted(cells.items())},
            'comparisons': details,
            'p_raw': raw_p,
            'p_holm': scoring.holm_adjust_family(raw_p, family_size),
        }

    path = paths.DELIV_DIR / 'extension_comparison.json'
    path.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding='utf-8')
    print('reuse verified: %s -- %s' % (reuse_ok, reuse_note))
    print('representations pooled in one table: %s' % ', '.join(representations))
    print()
    for task_name, block in out['per_task'].items():
        print(task_name)
        for key, detail in block['comparisons'].items():
            if detail.get('status') != 'tested':
                print('  %-22s NOT TESTED' % key)
                continue
            ci = detail.get('ci') or {}
            print('  %-22s n=%-4s %s: %.3f vs %.3f  diff=%+.3f  CI=[%.3f, %.3f]  '
                  'p_holm=%s'
                  % (key, detail.get('n_pairs'), detail.get('metric'),
                     detail.get('value_a'), detail.get('value_b'),
                     detail.get('observed_difference') or 0.0,
                     ci.get('low', float('nan')), ci.get('high', float('nan')),
                     _fmt_p(block['p_holm'].get(key))))
    print()
    print('wrote', path)
    return 0


def _fmt_p(value):
    return 'n/a' if value is None else '%.4f' % value


def cmd_score(args):
    _, run_config = load_configs()
    manifest = dataset.load_manifest()
    metrics = runner.score_run(args.run, run_config, manifest)
    print('scored %d cells from %s' % (len(metrics), args.run))
    return 0


def cmd_tables(args):
    _, run_config = load_configs()
    manifest = dataset.load_manifest()
    metrics = runner.score_run(args.run, run_config, manifest)
    representations = (args.representations.split(',') if args.representations
                       else list(reps.ROUND1_REPRESENTATIONS))
    summary, companion = report.write_summary_csv(args.run, metrics, representations, args.split)
    task_names = args.tasks.split(',') if args.tasks else list(tasks.TASK_ORDER)
    print('wrote', summary)
    print('wrote', companion)
    print('wrote', report.write_progress_matrix(
        metrics, report.run_state_map(args.run, representations, task_names)))
    return 0


def _write_run_manifest(run_id, mode, representations, task_names, model, manifest,
                        run_config, planned, budget=None):
    payload = {
        'run_id': run_id,
        'mode': mode,
        'created_at': datetime.datetime.now().isoformat(timespec='seconds'),
        'code': {
            'core_sha256': core_sha256(),
            'core_matches_baseline': core_matches_baseline(),
            'benchmark_version': __import__('bmb').__version__,
        },
        'model': {'provider': model.provider, 'model_id': model.model_id, 'region': model.region,
                  'base_url': model.base_url, 'generation': run_config['generation'],
                  'credential_variable': model.api_key_var_in_use,
                  'provider_terms_warning': client_mod.provider_terms_warning(model.provider),
                  'provider_terms_verbatim_cn': (
                      client_mod.TOKEN_PLAN_TERMS_CN if model.provider == 'token_plan' else None),
                  'provider_terms_source': (
                      client_mod.TOKEN_PLAN_TERMS_SOURCE if model.provider == 'token_plan'
                      else None)},
        'dataset': {'version': manifest['dataset_version'],
                    'counts': manifest['counts']},
        'representations': representations,
        'tasks': task_names,
        'planned_requests': planned,
        'seed': run_config['seed'],
        'budget': budget.snapshot() if budget else None,
        'scoring': {'qed_success_threshold': run_config['qed_success_threshold'],
                    'bootstrap': run_config['bootstrap']},
        'environment': {'python': sys.version.split()[0], 'rdkit': labels.rdkit_version()},
    }
    path = paths.run_dir(run_id) / 'experiment_manifest.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
    return path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest='command', required=True)

    p = sub.add_parser('prepare', help='build pool, splits and manifest')
    p.add_argument('--quarantine-runs', default=None,
                   help='comma-separated run ids whose molecules are removed from dev/test '
                        '(default: %s). Use an empty string to keep everything.'
                        % ','.join(DEBUG_RUNS))
    p.set_defaults(func=cmd_prepare)

    p = sub.add_parser('selftest', help='scoring checks with known/wrong answers')
    p.set_defaults(func=cmd_selftest)

    p = sub.add_parser('offline-checks', help='client checks that need no network')
    p.set_defaults(func=cmd_offline_checks)

    p = sub.add_parser('verify', help='round-trip the pool through every representation')
    p.add_argument('--representations', default=None)
    p.set_defaults(func=cmd_verify)

    p = sub.add_parser('dry-run', help='build prompts without calling a model')
    p.add_argument('--representations', default=None)
    p.add_argument('--tasks', default=None)
    p.add_argument('--split', default='test', choices=['dev', 'test'])
    p.add_argument('--run-id', default=None)
    p.set_defaults(func=cmd_dry_run)

    for mode in ('smoke', 'main'):
        p = sub.add_parser(mode, help='real model calls, %s size' % mode)
        p.add_argument('--representations', default=None)
        p.add_argument('--tasks', default=None)
        p.add_argument('--split', default='test', choices=['dev', 'test'])
        p.add_argument('--run-id', default=None)
        p.add_argument('--limit', type=int, default=None,
                       help='override the per-task item count for this run (e.g. 5 for a '
                            'connection check). It only sets how many items are attempted; it '
                            'does not change the frozen main plan.')
        p.set_defaults(func=lambda a, m=mode: cmd_run(a, m))

    p = sub.add_parser('extension-prep',
                       help='freeze item lists for the five added representations (no model)')
    p.set_defaults(func=cmd_extension_prep)

    p = sub.add_parser('extension-checks',
                       help='offline checks for the five added representations (no model)')
    p.set_defaults(func=cmd_extension_checks)

    p = sub.add_parser('extension-dev',
                       help='dev-side validation of the five added representations (real calls)')
    p.add_argument('--representations', default=None)
    p.add_argument('--run-id', default=None)
    p.set_defaults(func=cmd_extension_dev)

    p = sub.add_parser('extension-run',
                       help='real model calls for the five added representations')
    p.add_argument('--representations', default=None)
    p.add_argument('--tasks', default=None)
    p.add_argument('--run-id', default=None)
    p.add_argument('--limit', type=int, default=None,
                   help='attempt only the first N frozen items of each cell (a dev or '
                        'staged run); it never changes what was frozen')
    p.add_argument('--order-seed', type=int, default=20261006,
                   help='seed for the cell order; the items inside a cell stay frozen')
    p.set_defaults(func=cmd_extension_run)

    p = sub.add_parser('main-prep',
                       help='apply the complexity gate and freeze the main sample (no model)')
    p.set_defaults(func=cmd_main_prep)

    p = sub.add_parser('score', help='re-score a saved run')
    p.add_argument('--run', required=True)
    p.set_defaults(func=cmd_score)

    p = sub.add_parser('compare',
                       help='pre-declared paired comparisons with Holm correction (main only)')
    p.add_argument('--run', required=True)
    p.set_defaults(func=cmd_compare)

    p = sub.add_parser('ext-compare',
                       help='the nine-pair extension family across all eight rows')
    p.add_argument('--run', required=True, help='the extension run')
    p.add_argument('--reuse-run', default='main-02',
                   help='the run whose three inherited representations are reused, and only '
                        'when model, provider, endpoint and generation all match')
    p.set_defaults(func=cmd_ext_compare)

    p = sub.add_parser('tables', help='rebuild tables from a saved run')
    p.add_argument('--run', required=True)
    p.add_argument('--representations', default=None)
    p.add_argument('--tasks', default=None)
    p.add_argument('--split', default='test', choices=['dev', 'test'])
    p.set_defaults(func=cmd_tables)

    args = parser.parse_args(argv)
    load_local_env()
    return args.func(args)


if __name__ == '__main__':
    sys.exit(main())