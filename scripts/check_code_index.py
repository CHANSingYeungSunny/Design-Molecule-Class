"""Check that every row of docs/CORE_CODE_INDEX.md still points at real code.

The index is documentation about source locations, so it can silently rot when a
file is edited. This reads the table back and confirms that each listed line
range really contains the definition it claims. Standard library only.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INDEX = ROOT / 'docs/CORE_CODE_INDEX.md'

ROW = re.compile(r'^\|\s*`([^`]+)`\s*\|\s*`([^`]+)`\s*\|\s*(\d+)[-\u2013](\d+)\s*\|', re.M)


def main():
    if not INDEX.is_file():
        print('missing', INDEX)
        return 1
    rows = ROW.findall(INDEX.read_text(encoding='utf-8'))
    if not rows:
        print('no index rows found')
        return 1
    problems = []
    for path, symbol, start, end in rows:
        target = ROOT / path
        if not target.is_file():
            problems.append('%s: file missing' % path)
            continue
        lines = target.read_text(encoding='utf-8').splitlines()
        segment = '\n'.join(lines[int(start) - 1:int(end)])
        name = symbol.split('.')[-1]
        if ('def %s' % name) not in segment and ('class %s' % name) not in segment:
            problems.append('%s %s: lines %s-%s do not define %s'
                            % (path, symbol, start, end, name))
    print('index rows checked:', len(rows))
    for problem in problems:
        print('  PROBLEM', problem)
    print('core code index OK' if not problems else 'core code index OUT OF DATE')
    return 1 if problems else 0


if __name__ == '__main__':
    sys.exit(main())
