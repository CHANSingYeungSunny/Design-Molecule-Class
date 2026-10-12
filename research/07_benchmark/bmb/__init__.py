"""Molecular representation benchmark (MolRepBench-inspired extension).

Round 1 compares three representations -- Our Code, SMILES, SELFIES -- on one
model, one molecule list and one set of scoring rules. Image, Graph2D, Graph3D
and MolJSON are declared but not implemented this round; see reps.py.

The core molecule classes are never modified here. Scoring never runs through
them: RDKit produces the reference labels, so Our Code cannot mark its own work.
"""
from __future__ import annotations

__all__ = ['client', 'core_loader', 'dataset', 'labels', 'paths', 'report',
           'reps', 'runner', 'scoring', 'tasks']

__version__ = '1.0.0-round1'