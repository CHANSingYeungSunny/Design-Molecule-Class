"""Shared loader for the professor's baseline prototype (read-only baseline copy).

Only the standard library is used. The prototype is loaded by file path so that
the baseline copy in 00_baseline/ is the single object under test.
"""
from __future__ import annotations

import hashlib
import importlib.util
import pathlib
import sys

BASELINE_DIR = pathlib.Path(__file__).resolve().parent.parent / '00_baseline'
PROTOTYPE_PATH = BASELINE_DIR / 'molecule_standalone.py'
CORE_PATH = BASELINE_DIR.parent / '01_implementation' / 'molecule_core.py'
REPO_ROOT = BASELINE_DIR.parent.parent
EXPECTED_BASELINE_SHA256 = (
    '55d314aa312707b498536d339ff59dc89bbbe4956cab54a74a4746cd385af57b'
)


def load_prototype(module_name: str = 'professor_baseline_molecule_standalone'):
    """Import the baseline prototype under a private module name."""
    return _load(PROTOTYPE_PATH, module_name)


def load_core(module_name: str = 'our_code_core'):
    """Import the deliverable core (Atom/Bond/Molecule) under a private name.

    The core is the professor original; the benchmark must reach the molecule
    classes through here rather than growing private methods on them.
    """
    return _load(CORE_PATH, module_name)


def _load(path: pathlib.Path, module_name: str):
    if module_name in sys.modules:
        return sys.modules[module_name]
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def baseline_sha256() -> str:
    return hashlib.sha256(PROTOTYPE_PATH.read_bytes()).hexdigest()
