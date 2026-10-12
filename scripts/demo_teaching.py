"""A small teaching example: one molecule, its representations and its RDKit reference.

TEACHING ONLY. Ethanol (CCO) is far below the complexity gate of the published
benchmark (MW >= 300 Da and at least 2 rings), so nothing printed here enters the
formal 8 x 5 table, and no model is asked anything.

Run with the project environment:
    .venv\\Scripts\\python.exe scripts/demo_teaching.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / 'research/07_benchmark'
sys.path.insert(0, str(BENCH))

from bmb import labels, reps  # noqa: E402

TEACHING_SMILES = 'CCO'


def main():
    smiles = reps.canonical_smiles(reps.reference_mol(TEACHING_SMILES))
    reference = {
        'atom_counts': labels.atom_counts(smiles),
        'qed': round(labels.qed_reference(smiles), 6),
        'hbd': labels.descriptor_reference(smiles)['hbd'],
        'formula': labels.molecular_formula(smiles),
        'functional_groups': labels.functional_group_labels(smiles, only_enabled=True),
    }
    payload = {
        'status': 'teaching example; excluded from the published scores',
        'molecule': 'ethanol',
        'canonical_smiles': smiles,
        'our_code': reps.our_code_text(smiles),
        'smiles': reps.smiles_text(smiles),
        'selfies': reps.selfies_text(smiles),
        'rdkit_reference': reference,
        'rdkit_version': labels.rdkit_version(),
    }
    text = json.dumps(payload, ensure_ascii=False, indent=2)
    print(text)
    out = ROOT / 'reproduced'
    out.mkdir(parents=True, exist_ok=True)
    (out / 'teaching_example.json').write_text(text + '\n', encoding='utf-8')
    print('\nThis molecule is too small for the formal gate (MW >= 300 Da, >= 2 rings).')
    print('Written to', out / 'teaching_example.json')


if __name__ == '__main__':
    main()
