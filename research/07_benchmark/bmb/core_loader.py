"""Load the deliverable core (Atom / Bond / Molecule) for the benchmark.

The core is the professor's prototype, unchanged. The benchmark needs it for
exactly one thing: turning a reference molecule into the "Our Code" text that
the model reads. Scoring never goes through here -- reference labels come from
RDKit (see labels.py) so that Our Code cannot mark its own homework.

Python 3.12+ note: the module must be registered in sys.modules *before*
exec_module, otherwise @dataclass cannot resolve cls.__module__ for the
string annotations under `from __future__ import annotations`.
"""
from __future__ import annotations

import functools
import hashlib
import importlib.util
import sys

from . import paths


def core_sha256() -> str:
    return hashlib.sha256(paths.CORE_PATH.read_bytes()).hexdigest()


def core_matches_baseline() -> bool:
    return paths.CORE_PATH.read_bytes() == paths.BASELINE_PATH.read_bytes()


@functools.lru_cache(maxsize=1)
def load_core(module_name: str = 'bmb_our_code_core'):
    """Import the core under a private module name and verify it is the baseline."""
    if not core_matches_baseline():
        raise RuntimeError(
            '01_implementation/molecule_core.py differs from the professor '
            'baseline. The benchmark must not run against a modified core.')
    if module_name in sys.modules:
        return sys.modules[module_name]
    spec = importlib.util.spec_from_file_location(module_name, paths.CORE_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module