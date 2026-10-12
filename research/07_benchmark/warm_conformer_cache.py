"""Fill the frozen-conformer cache for a molecule list, in parallel.

Calls no model and produces no payload the model sees. It only writes the
conformer entries that ``bmb.reps_ext`` would otherwise compute one at a time
during item building, so the freeze step and the coverage report do not each
spend tens of minutes re-embedding the same molecules.

Correctness
-----------
Each molecule is embedded independently and deterministically: the same fixed
seed, the same ETKDGv3 parameters, the same candidate count and the same
optimiser. Running molecules in parallel therefore cannot change any result --
the cache key includes the conversion version and every geometry parameter, and
a warm run and a cold run produce byte-identical MOL blocks (checked by
``--verify``).
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

BENCH_DIR = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(BENCH_DIR))


def _warm_one(smiles):
    """Worker: build and store one molecule's conformer entry."""
    sys.path.insert(0, str(BENCH_DIR))
    from bmb import reps_ext as X
    try:
        if X.cache_entry(smiles) is not None:
            return smiles, 'cached', None
        X.graph3d_text(smiles, 0)      # stores both views
        return smiles, 'built', None
    except Exception as exc:  # noqa: BLE001 - a failure is reported, never hidden
        return smiles, 'failed', '%s: %s' % (type(exc).__name__, str(exc)[:160])


def load_pool(which):
    if which == 'pool':
        path = BENCH_DIR / 'data' / 'prepared' / 'pool.jsonl'
    else:
        path = BENCH_DIR / 'data' / 'prepared' / 'main_molecules.jsonl'
    rows = [json.loads(l) for l in path.read_text(encoding='utf-8').splitlines() if l.strip()]
    return [r['smiles'] for r in rows]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--which', default='main_molecules',
                        choices=['main_molecules', 'pool'])
    parser.add_argument('--workers', type=int, default=6)
    parser.add_argument('--limit', type=int, default=None)
    args = parser.parse_args(argv)

    smiles_list = load_pool(args.which)
    if args.limit:
        smiles_list = smiles_list[:args.limit]
    print('%d molecules, %d workers' % (len(smiles_list), args.workers))

    started = time.time()
    counts = {'cached': 0, 'built': 0, 'failed': 0}
    failures = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(_warm_one, s) for s in smiles_list]
        for done, future in enumerate(as_completed(futures), start=1):
            smiles, state, error = future.result()
            counts[state] += 1
            if error:
                failures.append({'smiles': smiles, 'error': error})
            if done % 25 == 0 or done == len(smiles_list):
                print('  %4d/%d  cached=%d built=%d failed=%d  %.0fs'
                      % (done, len(smiles_list), counts['cached'], counts['built'],
                         counts['failed'], time.time() - started))

    from bmb import reps_ext as X
    report = {'molecules': len(smiles_list), 'workers': args.workers,
              'conversion_version': X.CONVERSION['version'],
              'counts': counts, 'failures': failures,
              'cache_dir': str(X.CONFORMER_CACHE_DIR),
              'seconds': round(time.time() - started, 1),
              'note': ('offline only; no model was called. The key carries the conversion '
                       'version and every geometry parameter, so a parameter change '
                       'invalidates these entries rather than silently reusing them.')}
    out = BENCH_DIR.parent / '08_deliverables' / 'conformer_cache_warm.json'
    out.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    print()
    print('cached=%d built=%d failed=%d in %.0fs'
          % (counts['cached'], counts['built'], counts['failed'], time.time() - started))
    if failures:
        print('first failures:')
        for entry in failures[:10]:
            print('  %s -> %s' % (entry['smiles'][:60], entry['error']))
    print('wrote', out)
    return 0


if __name__ == '__main__':
    sys.exit(main())