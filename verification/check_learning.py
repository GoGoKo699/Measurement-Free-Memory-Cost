#!/usr/bin/env python3
"""Check the selected source, thirty-milestone map and local navigation.

These source-structure checks do not certify hosted rendering, scientific
truth or a reader's understanding. Run from any directory; --output saves a
fresh report without modifying a committed reference.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import unicodedata
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def without_code(text: str) -> str:
    return re.sub(r'^```[^\n]*\n.*?^```[ \t]*$', '', text, flags=re.M | re.S)


def heading_ids(text: str) -> set[str]:
    seen: dict[str, int] = {}
    result = set()
    for title in re.findall(r'^#{1,6}\s+(.+?)\s*#*$', without_code(text), re.M):
        title = re.sub(r'[\*_`]', '', title.strip()).lower()
        slug = ''.join(ch for ch in title
                       if ch in '-_ ' or unicodedata.category(ch)[0] in 'LN').replace(' ', '-')
        count = seen.get(slug, 0)
        seen[slug] = count + 1
        result.add(slug if count == 0 else f'{slug}-{count}')
    return result


def milestone_checks() -> dict:
    text = (ROOT / 'learning/MILESTONES.md').read_text(encoding='utf-8')
    matches = list(re.finditer(r'^#{2,3} M(\d{2})\s*(?:—|–|:|-)\s+(.+)$', text, re.M))
    require([int(match[1]) for match in matches] == list(range(1, 31)),
            'M01–M30 must occur once each and in order')
    edges = 0
    phases = {'background': 0, 'narrative': 0}
    for i, match in enumerate(matches):
        number = int(match[1])
        end = matches[i + 1].start() if i + 1 < len(matches) else len(text)
        body = text[match.end():end]
        for field in ('Phase', 'Central question', 'Understanding target',
                      'Conceptual connection', 'Source'):
            require(f'**{field}' in body, f'M{number:02}: missing {field}')
        require('Dependencies' in body and ('Story role' in body or 'story role' in body),
                f'M{number:02}: dependencies and story role required')
        phase = re.search(r'\*\*Phase:\*\*\s*([^.*]+)', body)
        require(phase is not None, f'M{number:02}: missing phase')
        phase_text = phase[1].lower().replace(' ', '')
        key = 'background' if 'background' in phase_text else 'narrative'
        require(key == 'background' or 'narrative' in phase_text,
                f'M{number:02}: unknown phase')
        phases[key] += 1
        dep_line = re.search(r'\*\*Dependencies(?: and story role)?:\*\*\s*([^*]+)', body)
        require(dep_line is not None, f'M{number:02}: missing dependencies')
        deps = [int(dep) for dep in re.findall(r'M(\d{2})', dep_line[1].split(';')[0])]
        require(all(1 <= dep < number for dep in deps), f'M{number:02}: forward dependency')
        require(len(deps) == len(set(deps)), f'M{number:02}: repeated dependency')
        require(bool(deps) or number == 1, f'M{number:02}: disconnected milestone')
        edges += len(deps)
    require(phases['background'] > 0 and phases['narrative'] > 0, 'Both learning phases required')
    bridge = (ROOT / 'learning/WATROUS_BRIDGE.md').read_text(encoding='utf-8')
    require('Watrous' in bridge and '2018' in bridge, 'Selected Watrous2018 bridge missing')
    require('Watrous' in text and 'GOTTESMAN_BRIDGE.md' not in text,
            'Milestones must use the selected common source')
    require('learning/WATROUS_BRIDGE.md' in (ROOT / 'README.md').read_text(),
            'Root selected-source route missing')
    require('learning/WATROUS_BRIDGE.md' in (ROOT / 'llms.txt').read_text(),
            'LLM selected-source route missing')
    return {'milestones': 30, 'dependency_edges': edges, **phases}


def navigation_checks() -> dict:
    files = sorted(path for path in ROOT.rglob('*.md')
                   if not any(part.startswith('.') for part in path.relative_to(ROOT).parts))
    links, fragments = 0, 0
    for path in files:
        text = path.read_text(encoding='utf-8')
        require(len(re.findall(r'^```', text, re.M)) % 2 == 0,
                f'{path.relative_to(ROOT)}: unbalanced fences')
        for target in re.findall(r'\]\(([^\s)]+)\)', text):
            parsed = urlsplit(target)
            if parsed.scheme or parsed.netloc:
                continue
            candidate = (path.parent / unquote(parsed.path)).resolve() if parsed.path else path
            require(candidate.is_relative_to(ROOT), f'{path}: link leaves repository: {target}')
            require(candidate.exists(), f'{path}: missing destination: {target}')
            links += 1
            if parsed.fragment:
                require(candidate.suffix == '.md', f'{path}: unsupported fragment: {target}')
                require(unquote(parsed.fragment) in heading_ids(candidate.read_text()),
                        f'{path}: missing heading: {target}')
                fragments += 1
    return {'markdown_files': len(files), 'local_links': links, 'heading_fragments': fragments}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = {'status': 'passed',
              'anchor': 'Watrous, The Theory of Quantum Information (2018)',
              'milestone_map': milestone_checks(), 'navigation': navigation_checks(),
              'hosted_rendering': 'not_tested_by_this_script'}
    report = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(report, encoding='utf-8')
    print(report, end='')


if __name__ == '__main__':
    main()
