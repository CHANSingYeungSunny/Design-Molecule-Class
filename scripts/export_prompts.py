"""Export the original prompt builder without calling a model."""
from __future__ import annotations
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / 'research/07_benchmark'
sys.path.insert(0, str(BENCH))
from bmb import tasks, labels  # noqa: E402


def write(path, text):
    # newline='' keeps '\n' as written, so the export is byte-identical on
    # every platform instead of picking up Windows CRLF.
    with path.open('w', encoding='utf-8', newline='') as stream:
        stream.write(text)


def main():
    destination = ROOT / 'prompts'
    destination.mkdir(exist_ok=True)
    system = tasks.system_message()['content']
    write(destination / 'system.txt', system + '\n')
    export = {
        'source': 'research/07_benchmark/bmb/tasks.py',
        'status': 'historical prompts preserved verbatim; T2/T3 example double braces retained',
        'system': system, 'representation_notes': tasks.REP_NOTES,
        'asked_functional_groups': [entry['name'] for entry in labels.enabled_labels()],
        'scored_functional_groups': labels.frozen_macro_f1_labels(),
        'templates': {},
    }
    for name in tasks.TASK_ORDER:
        payload = ({'a': '{molecule_A}', 'b': '{molecule_B}'} if name == 't1_identity'
                   else {'molecule': '{molecule_payload}'})
        item = {'task': name, 'payload': payload}
        full = tasks.get(name).messages('our_code', item, {})[-1]['content']
        export['templates'][name] = {'user_with_our_code_note': full,
                                    'task_instruction': full.split('\n\n')[-1]}
    write(destination / 'benchmark_templates.json',
          json.dumps(export, ensure_ascii=False, indent=2))
    examples = []
    raw = BENCH / 'runs/main-02/raw_responses.jsonl'
    for line in raw.read_text(encoding='utf-8').splitlines():
        record = json.loads(line)
        if record['representation'] == 'our_code' and record['task'] not in {x['task'] for x in examples}:
            examples.append({'task': record['task'], 'example_id': record['example_id'],
                             'prompt_hash': record['prompt_hash'],
                             'messages': [tasks.system_message(), {'role': 'user', 'content': record['prompt_text']}],
                             'model_raw_response': record['raw_response'], 'reference': record['label'],
                             'warning': 'Reference is separate evidence and was never included in messages.'})
    write(destination / 'historical_our_code_examples.json',
          json.dumps(examples, ensure_ascii=False, indent=2))
    print('Exported the system prompt, eight format notes, five task templates, and five historical examples.')


if __name__ == '__main__':
    main()
