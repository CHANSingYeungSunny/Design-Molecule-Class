"""The run loop: one cell = (representation x task x split x model x mode).

Guarantees
----------
* Every scheduled item ends with a status. Nothing is dropped silently.
* Resumable: raw responses are appended to JSONL as they arrive, so an
  interrupted run continues instead of restarting.
* Scoring reads only saved responses, so tables can be rebuilt without spending
  another token.
* Exactly one session per item: the model never sees another representation's
  answer. Call order is shuffled with a recorded seed so no representation
  systematically benefits from a warm cache or a lucky batch.
* A transport failure that survives the declared retries is 'pending'; the cell
  is marked incomplete rather than reported with a smaller denominator.
* A wrong answer is never retried and never replaced.
"""
from __future__ import annotations

import pathlib
import hashlib
import json
import random
import time

from . import client as client_mod
from . import paths, scoring, tasks


def generation_config(run_config):
    gen = dict(run_config.get('generation', {}))
    gen.pop('_comment', None)
    return gen


def cache_key(record_meta, run_config):
    blob = json.dumps({
        'dataset_version': record_meta['dataset_version'],
        'example_id': record_meta['example_id'],
        'representation': record_meta['representation'],
        'task': record_meta['task'],
        'prompt_hash': record_meta['prompt_hash'],
        # an image or a conformer is not captured by prompt text alone, so the
        # payload digest and the pixel / coordinate digests are part of the key
        'payload_sha256': record_meta.get('payload_sha256'),
        'image_sha256': record_meta.get('image_sha256'),
        'coords_sha256': record_meta.get('coords_sha256'),
        # the CHANNEL belongs in the key: the same question answered over a
        # different provider or endpoint is a different observation
        'provider': record_meta.get('provider'),
        'base_url': record_meta.get('base_url'),
        'model_id': record_meta['model_id'],
        'generation': generation_config(run_config),
    }, sort_keys=True)
    return hashlib.sha256(blob.encode('utf-8')).hexdigest()


def _read_jsonl(path):
    if not path.is_file():
        return []
    with open(path, encoding='utf-8') as handle:
        return [json.loads(line) for line in handle if line.strip()]


def _append_jsonl(path, record):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'a', encoding='utf-8') as handle:
        handle.write(json.dumps(record, ensure_ascii=False) + '\n')


TERMINAL = {'ok', 'invalid'}


class Budget:
    """Hard stop on requests, tokens and money.

    Three request counts are kept apart, because they answer different questions:

      basic_requests     planned items we actually attempted
      transport_attempts every HTTP attempt, so a retry costs one
      completed_items    items that came back with a terminal answer

    A run that retries heavily therefore cannot hide behind the basic count.

    Uncertainty is handled conservatively in one direction: when a response
    carries no usage figure, the attempt is charged a RESERVE rather than
    zero, because the provider did perform the work and will bill for it.
    """

    def __init__(self, max_basic_requests=None, max_transport_attempts=None,
                 max_total_tokens=None, max_completion_tokens=None, max_rmb=None,
                 in_flight_reserve=1, prices=None, unknown_usage_token_reserve=8000):
        self.max_basic_requests = max_basic_requests
        self.max_transport_attempts = max_transport_attempts
        self.max_total_tokens = max_total_tokens
        self.max_completion_tokens = max_completion_tokens
        self.max_rmb = max_rmb
        self.in_flight_reserve = in_flight_reserve
        self.prices = prices or {}
        self.unknown_usage_token_reserve = unknown_usage_token_reserve
        self.basic_requests = 0
        self.transport_attempts = 0
        self.completed_items = 0
        self.prompt_tokens = 0
        self.completion_tokens = 0
        self.other_tokens = 0
        self.attempts_without_usage = 0
        self.stopped_reason = None

    # -- accounting --------------------------------------------------------
    @property
    def total_tokens(self):
        return self.prompt_tokens + self.completion_tokens + self.other_tokens

    def note_basic_request(self):
        self.basic_requests += 1

    def note_completed(self):
        self.completed_items += 1

    def record(self, attempts, usage):
        """Charge one logical request. attempts include retries; usage may be absent."""
        made = max(1, len(attempts or []))
        self.transport_attempts += made
        if usage:
            prompt = usage.get('prompt_tokens') or 0
            completion = usage.get('completion_tokens') or 0
            total = usage.get('total_tokens') or (prompt + completion)
            self.prompt_tokens += prompt
            self.completion_tokens += completion
            self.other_tokens += max(0, total - prompt - completion)
        else:
            # no usage figure: charge a conservative reserve, never zero
            self.attempts_without_usage += made
            self.other_tokens += self.unknown_usage_token_reserve * made

    # -- gating ------------------------------------------------------------
    def _reserve(self):
        """Room kept for requests already in flight when concurrency is used."""
        return max(0, self.in_flight_reserve or 0)

    def can_continue(self):
        checks = [
            (self.max_transport_attempts,
             self.transport_attempts + self._reserve(), 'transport attempt cap'),
            (self.max_basic_requests,
             self.basic_requests + self._reserve(), 'basic request cap'),
            (self.max_total_tokens, self.total_tokens, 'total token cap'),
            (self.max_completion_tokens, self.completion_tokens, 'completion token cap'),
        ]
        for cap, used, label in checks:
            if cap is not None and used >= cap:
                self.stopped_reason = '%s reached (%s of %s)' % (label, used, cap)
                return False
        rmb = self.estimated_rmb()
        if self.max_rmb is not None and rmb is not None and rmb >= self.max_rmb:
            self.stopped_reason = ('RMB cap reached (%.2f of %.2f CNY, busy-hour pricing)'
                                   % (rmb, self.max_rmb))
            return False
        return True

    def estimated_rmb(self, busy=True):
        """RMB at list price, or None when no price applies to this provider.

        Returns None (not 0.0) for a subscription provider: zero would read as
        'this run cost nothing', which is a claim we cannot support.
        """
        p = self.prices
        if not p:
            return None
        rate_in = p.get('input_busy' if busy else 'input_idle')
        rate_out = p.get('output_busy' if busy else 'output_idle')
        if rate_in is None or rate_out is None:
            return None
        return ((self.prompt_tokens + self.other_tokens) / 1e6 * rate_in
                + self.completion_tokens / 1e6 * rate_out)

    def snapshot(self):
        return {
            'basic_requests': self.basic_requests,
            'max_basic_requests': self.max_basic_requests,
            'transport_attempts': self.transport_attempts,
            'max_transport_attempts': self.max_transport_attempts,
            'completed_items': self.completed_items,
            'prompt_tokens': self.prompt_tokens,
            'completion_tokens': self.completion_tokens,
            'other_tokens': self.other_tokens,
            'total_tokens': self.total_tokens,
            'max_total_tokens': self.max_total_tokens,
            'max_completion_tokens': self.max_completion_tokens,
            'attempts_without_usage': self.attempts_without_usage,
            'estimated_rmb_busy_hour': self.estimated_rmb(busy=True),
            'max_rmb': self.max_rmb,
            'stopped_reason': self.stopped_reason,
        }


def budget_from_extension(extension_config, provider='bailian'):
    """The EXTENSION stage has its own bound.

    The previous cycle's ceiling is neither reset nor reused as this stage's
    limit, so the two ledgers cannot be confused with one another.
    """
    budget_cfg = extension_config.get('budget') or {}
    return Budget(
        max_basic_requests=budget_cfg.get('max_basic_requests'),
        max_transport_attempts=budget_cfg.get('max_transport_attempts'),
        max_total_tokens=budget_cfg.get('max_total_tokens'),
        max_completion_tokens=budget_cfg.get('max_completion_tokens'),
        max_rmb=budget_cfg.get('max_rmb'),
        in_flight_reserve=budget_cfg.get('in_flight_reserve', 0),
        prices=None,          # subscription provider: no RMB figure is computed
    )


def budget_from_config(run_config, provider='bailian'):
    budget_cfg = run_config.get('budget') or {}
    # The published per-token RMB rates describe the pay-as-you-go model API.
    # A subscription provider consumes plan credits instead, so no RMB figure
    # is computed for it: the amount is UNKNOWN unless the platform reports
    # credits, and inventing a conversion would be a fabricated number.
    prices = ((run_config.get('pricing') or {}).get('beijing')
              if provider == 'bailian' else None)
    return Budget(
        max_basic_requests=budget_cfg.get('max_basic_requests'),
        max_transport_attempts=budget_cfg.get('max_transport_attempts'),
        max_total_tokens=budget_cfg.get('max_total_tokens'),
        max_completion_tokens=budget_cfg.get('max_completion_tokens'),
        max_rmb=budget_cfg.get('max_rmb'),
        in_flight_reserve=budget_cfg.get('in_flight_reserve', 1),
        prices=prices,
    )


class FrozenListError(RuntimeError):
    """The frozen item list cannot be trusted, so no request may be sent."""


def load_frozen_items(representation, task_name, manifest, verify_prompts=True,
                      run_config=None):
    """Load a frozen main item list and verify it before it is used.

    main executes the FROZEN LIST, not items rebuilt from a pool. Rebuilding is
    what silently dropped the complexity gate before; reading the file the
    manifest points at removes that whole class of error.

    Checks, all of which must pass:
      1. the file exists and is named for exactly this representation x task;
      2. its sha256 equals manifest.frozen.item_list_hashes[task][rep];
      3. every item's task/rep/split match this cell;
      4. every item's prompt hash is one the manifest recorded for this cell;
      5. the molecules the items name are the ones the frozen pool holds.
    """
    path = paths.prepared('main_items') / ('%s-%s.jsonl' % (representation, task_name))
    if not path.is_file():
        raise FrozenListError('frozen item list missing: %s' % path)
    blob = path.read_text(encoding='utf-8')
    digest = hashlib.sha256('\n'.join(
        json.dumps(json.loads(l), ensure_ascii=False, sort_keys=True)
        for l in blob.splitlines() if l.strip()).encode('utf-8')).hexdigest()

    # main_manifest.json IS the frozen block; a run manifest nests it under
    # 'frozen'. Accept either so both callers work.
    frozen = manifest.get('frozen') if 'frozen' in manifest else manifest
    frozen = frozen or {}
    recorded = ((frozen.get('item_list_hashes') or {}).get(task_name) or {}).get(representation)
    if recorded is None:
        raise FrozenListError(
            'the manifest records no item list for %s/%s; run main-prep first'
            % (representation, task_name))
    if digest != recorded:
        raise FrozenListError(
            'frozen item list hash mismatch for %s/%s\n  file    : %s\n  manifest: %s'
            % (representation, task_name, digest, recorded))

    items = [json.loads(l) for l in blob.splitlines() if l.strip()]
    allowed_hashes = set(
        ((frozen.get('prompt_hashes') or {}).get(representation) or {}).get(task_name) or [])
    problems = []
    for item in items:
        if item.get('task') != task_name or item.get('rep') != representation:
            problems.append('%s belongs to %s/%s' % (item.get('example_id'),
                                                     item.get('rep'), item.get('task')))
            continue
        if verify_prompts and allowed_hashes:
            if run_config is None:
                raise FrozenListError('run_config is required to verify prompt hashes')
            current = tasks.prompt_hash(representation, item, run_config)
            if current not in allowed_hashes:
                problems.append('%s prompt hash not among the frozen hashes'
                                % item.get('example_id'))
    if problems:
        raise FrozenListError(
            'the frozen list for %s/%s does not match the manifest:\n  %s'
            % (representation, task_name, '\n  '.join(problems[:10])))
    return items


def load_extension_items(representation, task_name, manifest, run_config=None):
    """Load a frozen EXTENSION item list and verify it before it is used.

    Same discipline as ``load_frozen_items``: the run executes the file the
    manifest points at, never a rebuild, so a conversion that fails on one
    molecule cannot be quietly re-drawn later. An image payload is stored as a
    path plus sha256 rather than inline; it is re-read here and the digest is
    re-checked, so a replaced or corrupted PNG stops the run instead of being
    sent.
    """
    recorded = ((manifest.get('item_lists') or {}).get(task_name) or {}).get(representation)
    if recorded is None:
        raise FrozenListError(
            'the extension manifest records no item list for %s/%s; run extension-prep first'
            % (representation, task_name))
    path = paths.prepared(recorded['pool'])
    if not path.is_file():
        raise FrozenListError('frozen extension item list missing: %s' % path)
    items = [json.loads(l) for l in path.read_text(encoding='utf-8').splitlines() if l.strip()]
    digest = hashlib.sha256('\n'.join(
        json.dumps(i, ensure_ascii=False, sort_keys=True) for i in items).encode('utf-8')).hexdigest()
    if digest != recorded['sha256']:
        raise FrozenListError(
            'frozen extension item list hash mismatch for %s/%s\n  file    : %s\n  manifest: %s'
            % (representation, task_name, digest, recorded['sha256']))

    problems = []
    for item in items:
        if item.get('task') != task_name or item.get('rep') != representation:
            problems.append('%s belongs to %s/%s'
                            % (item.get('example_id'), item.get('rep'), item.get('task')))
            continue
        for key, value in list(item['payload'].items()):
            if isinstance(value, dict) and 'image_path' in value:
                blob_path = paths.BENCH_DIR.joinpath(*pathlib.PureWindowsPath(value['image_path']).parts)
                if not blob_path.is_file():
                    problems.append('%s image missing: %s' % (item['example_id'], blob_path))
                    continue
                data = blob_path.read_bytes()
                if hashlib.sha256(data).hexdigest() != value['sha256']:
                    problems.append('%s image digest mismatch' % item['example_id'])
                    continue
                item['payload'][key] = data
    # the prompt each item produces TODAY must be one of the prompts the manifest
    # recorded for this cell; a change to the note, the instruction or the
    # renderer therefore stops the run instead of silently sending a new stimulus
    allowed = set(((manifest.get('prompt_hashes') or {}).get(representation) or {})
                  .get(task_name) or [])
    if allowed:
        if run_config is None:
            raise FrozenListError('run_config is required to verify extension prompt hashes')
        for item in items:
            if tasks.prompt_hash(representation, item, run_config) not in allowed:
                problems.append('%s prompt hash not among the frozen hashes'
                                % item.get('example_id'))
    if problems:
        raise FrozenListError('the frozen extension list for %s/%s is not usable:\n  %s'
                              % (representation, task_name, '\n  '.join(problems[:10])))
    return items


def payload_digests(representation, item):
    """Digests that identify the exact bytes a request carries."""
    from . import reps_ext
    out = {'payload_sha256': None, 'image_sha256': None, 'coords_sha256': None}
    payload = item.get('payload') or {}
    parts = []
    for key in sorted(payload):
        value = payload[key]
        if isinstance(value, bytes):
            digest = hashlib.sha256(value).hexdigest()
            if representation == 'img':
                out['image_sha256'] = out['image_sha256'] or digest
            parts.append(digest)
        else:
            parts.append(hashlib.sha256(str(value).encode('utf-8')).hexdigest())
    out['payload_sha256'] = hashlib.sha256(''.join(parts).encode('utf-8')).hexdigest()
    if representation == 'graph3d' and isinstance(payload.get('molecule'), str):
        # the coordinates ARE the payload for this representation, so the same
        # digest answers the 'did the conformer change' question directly
        out['coords_sha256'] = out['payload_sha256']
    return out


def run_cell(client, representation, task_name, split_name, count, run_id,
             run_config, manifest, dry_run=False, budget_left=None, budget=None,
             pool_name=None, frozen_items=None):
    """Execute one cell and return the records (also persisted).

    ``pool_name`` selects WHICH molecule list the items are drawn from. `main`
    must pass the frozen complexity-gated set, not the raw split: the two are
    different molecule lists and the frozen one is the whole point of the
    complexity gate. Defaulting to the split name keeps smoke behaviour as it
    was.
    """
    task = tasks.get(task_name)
    if frozen_items is not None:
        # main executes the VERIFIED FROZEN LIST verbatim; nothing is rebuilt
        # and nothing is reshuffled, so the executed items are the frozen ones
        items = list(frozen_items)
    else:
        pool = _load_split(pool_name or split_name)
        seed = (run_config['seed'] if task_name != 't1_identity'
                else run_config['identity_seed'])
        items = task.build_items(pool, count, seed, representation, split_name)
        random.Random(seed).shuffle(items)

    raw_path = paths.run_dir(run_id) / 'raw_responses.jsonl'
    # Resume must not silently reuse answers produced by a DIFFERENT model,
    # endpoint or prompt. A stored record only counts as done when its model id
    # and prompt hash match what this run would send now.
    existing = set()
    ignored = 0
    for r in _read_jsonl(raw_path):
        if (r.get('representation') != representation or r.get('task') != task_name
                or r.get('status') not in TERMINAL):
            continue
        if r.get('model_id') != client.model_id:
            ignored += 1
            continue
        if r.get('provider') not in (None, getattr(client, 'provider', None)):
            ignored += 1
            continue
        if r.get('base_url') not in (None, getattr(client, 'base_url', None)):
            ignored += 1
            continue
        existing.add(r['example_id'])
    if ignored:
        print('  note: %d stored response(s) for %s/%s were produced by a different model, '
              'provider or endpoint and are NOT reused'
              % (ignored, representation, task_name))
    if existing:
        # prompt hashes are recomputed per item below, but a mismatch at any
        # item means the prompt changed; flag it loudly rather than reuse it
        stale = 0
        for r in _read_jsonl(raw_path):
            if r.get('example_id') not in existing:
                continue
            item = next((i for i in items if i['example_id'] == r['example_id']), None)
            if item is None:
                continue
            if tasks.prompt_hash(representation, item, run_config) != r.get('prompt_hash'):
                existing.discard(r['example_id'])
                stale += 1
        if stale:
            print('  note: %d stored response(s) were produced by a DIFFERENT prompt and are '
                  'NOT reused' % stale)

    records = []
    requests_made = 0
    for item in items:
        if item['example_id'] in existing:
            continue
        messages = task.messages(representation, item, run_config)
        meta = {
            'run_id': run_id,
            'dataset_version': manifest['dataset_version'],
            'example_id': item['example_id'],
            'split': split_name,
            'representation': representation,
            'task': task_name,
            'model_id': client.model_id,
            'provider': getattr(client, 'provider', None),
            'base_url': getattr(client, 'base_url', None),
            'prompt_hash': tasks.prompt_hash(representation, item, run_config),
            'cluster': item['cluster'],
            # every dependence unit this item belongs to. An Identity pair owns
            # both endpoint molecules; a single-molecule item owns one.
            'clusters': list(item.get('clusters') or (item['cluster'],)),
            # which molecule list the item was drawn from, so a run can be
            # audited against the frozen set rather than trusted by construction
            'pool': pool_name or split_name,
        }
        if item.get('payload'):
            meta.update(payload_digests(representation, item))
        meta['cache_key'] = cache_key(meta, run_config)

        if dry_run:
            record = dict(meta, status='dry_run', raw_response=None, parsed=None,
                          label=item['label'], prompt_text=messages[-1]['content'],
                          usage=None, latency_s=None, attempts=[], transport_error=None)
            records.append(record)
            continue

        if budget is not None and not budget.can_continue():
            record = dict(meta, status='pending', raw_response=None, parsed=None,
                          label=item['label'], prompt_text=messages[-1]['content'],
                          usage=None, latency_s=None, attempts=[],
                          transport_error='run stopped: %s' % budget.stopped_reason)
            records.append(record)
            _append_jsonl(raw_path, record)
            continue

        if budget_left is not None and requests_made >= budget_left:
            record = dict(meta, status='pending', raw_response=None, parsed=None,
                          label=item['label'], prompt_text=messages[-1]['content'],
                          usage=None, latency_s=None, attempts=[],
                          transport_error='run stopped: request budget reached')
            records.append(record)
            _append_jsonl(raw_path, record)
            continue

        if budget is not None:
            budget.note_basic_request()
        result = client.chat(messages)
        requests_made += 1
        if budget is not None:
            # a retry is a request: count every transport attempt and the
            # provider's own token figure. A missing usage figure is charged a
            # reserve inside record(), never zero.
            budget.record(result.get('attempts'), result.get('usage'))
            if result['ok'] or result.get('transport_error'):
                budget.note_completed()
        if not result['ok']:
            record = dict(meta, status='pending', raw_response=None, parsed=None,
                          label=item['label'], prompt_text=messages[-1]['content'],
                          usage=None, latency_s=result.get('latency_s'),
                          attempts=result.get('attempts', []),
                          transport_error=result.get('transport_error'))
        else:
            parsed, error = task.parse(result['text'])
            finish = result.get('finish_reason')
            # A length-stopped response is a configuration problem, not a wrong
            # answer. Flagging it separately keeps it diagnosable instead of
            # hiding inside a generic 'invalid'.
            truncated = finish == 'length'
            if truncated and error is not None:
                error = 'response truncated (finish_reason=length); %s' % error
            record = dict(meta,
                          status='ok' if error is None else 'invalid',
                          raw_response=result['text'], parsed=parsed, label=item['label'],
                          prompt_text=messages[-1]['content'],
                          usage=result.get('usage'), latency_s=result.get('latency_s'),
                          attempts=result.get('attempts', []),
                          finish_reason=finish, truncated=truncated,
                          parse_error=error, transport_error=None)
        _append_jsonl(raw_path, record)
        records.append(record)
        print('  %-58s %s' % (item['example_id'], record['status']))
    return records


def _load_split(split_name):
    from . import dataset
    return dataset.load_prepared(split_name)


def load_run_records(run_id, representation=None, task_name=None):
    records = _read_jsonl(paths.run_dir(run_id) / 'raw_responses.jsonl')
    if representation:
        records = [r for r in records if r.get('representation') == representation]
    if task_name:
        records = [r for r in records if r.get('task') == task_name]
    return records


def rebuild_item_index(run_id, run_config, split_name, pool_name=None):
    """Re-derive every scheduled item and index it by example_id.

    Item construction is seeded and the pool is frozen by hash, so this
    reproduces exactly the items a run was built from. It is used to attach
    dependence units to responses recorded before those units were stored.
    The raw responses themselves are never rewritten.
    """
    from . import dataset
    mode = run_id.split('-')[0]
    mode_cfg = run_config.get(mode) or run_config.get('smoke')
    if mode not in ('smoke', 'main', 'probe'):
        mode_cfg = run_config.get('smoke')
    index = {}
    try:
        pool = dataset.load_prepared(pool_name or split_name)
    except FileNotFoundError:
        return index
    for representation in ('our_code', 'smiles', 'selfies'):
        for task_name in tasks.TASK_ORDER:
            count = mode_cfg.get('per_task', 8)
            if task_name == 't1_identity':
                count = mode_cfg.get('identity_positives', 4) + mode_cfg.get('identity_negatives', 4)
            seed = (run_config.get('identity_seed', run_config['seed'])
                    if task_name == 't1_identity' else run_config['seed'])
            try:
                items = tasks.get(task_name).build_items(pool, count, seed,
                                                         representation, split_name)
            except Exception:  # noqa: BLE001
                continue
            for item in items:
                index[item['example_id']] = item
    return index


def debug_used_clusters(run_ids, run_config, splits=('dev', 'test')):
    """Every molecule any of the named runs touched, BOTH endpoints included.

    Records saved before `clusters` was written only carry one endpoint, so the
    deterministic item index is consulted to recover the other one. Quarantining
    on the single stored field would silently leave pair partners in the pool.
    """
    used = set()
    for run_id in run_ids:
        raw = paths.run_dir(run_id) / 'raw_responses.jsonl'
        if not raw.is_file():
            continue
        rows = [json.loads(l) for l in raw.read_text(encoding='utf-8').splitlines() if l.strip()]
        for split_name in splits:
            index = rebuild_item_index(run_id, run_config, split_name)
            for row in rows:
                if row.get('split') != split_name:
                    continue
                owned = row.get('clusters')
                if not owned:
                    item = index.get(row['example_id']) or {}
                    owned = item.get('clusters') or ([row['cluster']] if row.get('cluster') else [])
                used.update(unit for unit in owned if unit)
    return used


def score_run(run_id, run_config, manifest):
    """Rebuild every metric block from stored responses; no model calls."""
    from . import labels
    # Prefer the FROZEN macro-F1 label set. Re-deriving the labels from the
    # sample's realised positives would admit a label below the declared
    # support floor and change the mean's denominator between samples.
    label_names = labels.frozen_macro_f1_labels() or [
        item['name'] for item in labels.enabled_labels()]
    threshold = run_config.get('qed_success_threshold', 0.05)
    resamples = run_config.get('bootstrap', {}).get('resamples', 2000)

    records = load_run_records(run_id)

    # backfill dependence units for responses saved before they were recorded
    pools = {(r.get('split'), r.get('pool')) for r in records}
    index = {}
    for split_name, pool_name in pools:
        index.update(rebuild_item_index(run_id, run_config, split_name, pool_name))
    backfilled = 0
    for record in records:
        if record.get('clusters'):
            continue
        item = index.get(record['example_id'])
        if item:
            record['clusters'] = list(item.get('clusters') or (item['cluster'],))
            backfilled += 1
    if backfilled:
        print('  backfilled dependence units for %d stored responses '
              '(derived, not from the raw file)' % backfilled)
    cells = {}
    for record in records:
        cells.setdefault((record['representation'], record['task'], record['split']), []).append(record)

    out = {}
    for (representation, task_name, split_name), rows in sorted(cells.items()):
        block = scoring.evaluate(rows, task_name, label_names, threshold,
                                 resamples=resamples, seed=run_config['seed'])
        out['%s|%s|%s' % (representation, task_name, split_name)] = block

    scored_path = paths.run_dir(run_id) / 'scored.jsonl'
    with open(scored_path, 'w', encoding='utf-8') as handle:
        for record in records:
            detail = scoring.score_item(record)
            handle.write(json.dumps(dict(record, per_item_score=detail),
                                    ensure_ascii=False) + '\n')
    with open(paths.run_dir(run_id) / 'metrics.json', 'w', encoding='utf-8') as handle:
        json.dump(out, handle, ensure_ascii=False, indent=2)
    return out