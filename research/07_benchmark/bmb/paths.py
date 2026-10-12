"""Canonical paths for the benchmark module.

Everything is derived from this file's location so the module can be run from
any working directory.
"""
from __future__ import annotations

import pathlib

BENCH_DIR = pathlib.Path(__file__).resolve().parent.parent      # 07_benchmark/
RESEARCH_DIR = BENCH_DIR.parent                                  # research/
PROJECT_ROOT = RESEARCH_DIR.parent                               # project root

# --- the deliverable core (professor original, byte-identical to the baseline)
CORE_PATH = RESEARCH_DIR / '01_implementation' / 'molecule_core.py'
BASELINE_PATH = RESEARCH_DIR / '00_baseline' / 'molecule_standalone.py'
EXPECTED_CORE_SHA256 = '55d314aa312707b498536d339ff59dc89bbbe4956cab54a74a4746cd385af57b'

# --- benchmark data
DATA_DIR = BENCH_DIR / 'data'
RAW_DIR = DATA_DIR / 'raw'
PREPARED_DIR = DATA_DIR / 'prepared'

# --- configuration (contains NO secrets)
CONFIG_DIR = BENCH_DIR / 'config'
MODEL_CONFIG = CONFIG_DIR / 'model_config.json'
FUNCTIONAL_GROUPS = CONFIG_DIR / 'functional_groups.json'
# Written by `prepare` from measured support. This is the OPERATIVE label set:
# tasks and scoring read it in preference to the candidate file above.
FUNCTIONAL_GROUPS_FROZEN = CONFIG_DIR / 'functional_groups.frozen.json'
RUN_CONFIG = CONFIG_DIR / 'run_config.json'

# --- outputs
RUNS_DIR = BENCH_DIR / 'runs'
PROMPT_DIR = BENCH_DIR / 'prompts'

# --- deliverables
DELIV_DIR = RESEARCH_DIR / '08_deliverables'


def prepared(name: str) -> pathlib.Path:
    return PREPARED_DIR / name


def run_dir(run_id: str) -> pathlib.Path:
    return RUNS_DIR / run_id