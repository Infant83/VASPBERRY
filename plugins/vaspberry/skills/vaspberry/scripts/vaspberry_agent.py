#!/usr/bin/env python3
"""Local VASPBERRY diagnostics, exact-source credit, and released analytic demo."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys

REPOSITORY = 'https://github.com/Infant83/VASPBERRY'
TARGET_VERSION = '1.6.6'
TARGET_COMMIT = '2e7067b8b8b6d3ab444eef5bc025f4ce5720ca4d'
HISTORICAL_PRB_DOI = '10.1103/PhysRevB.93.041404'
HISTORICAL_PRB_SOURCE = (REPOSITORY + '/blob/'
                         '134b71ee69c9fa74215ce3b812d1b561070c7b8d/README.md#L158')
HISTORICAL_PRB_BIB = '''@article{PhysRevB.93.041404,
  title = {Competing magnetic orderings and tunable topological states in two-dimensional hexagonal organometallic lattices},
  author = {Kim, Hyun-Jung and Li, Chaokai and Feng, Ji and Cho, Jun-Hyung and Zhang, Zhenyu},
  journal = {Phys. Rev. B},
  volume = {93},
  pages = {041404},
  year = {2016},
  doi = {10.1103/PhysRevB.93.041404},
  url = {https://doi.org/10.1103/PhysRevB.93.041404}
}'''
METHODS = {
    'fhs': ('Fukui2005', 'Fukui, Takahiro and Hatsugai, Yasuhiro and Suzuki, Hiroshi',
            'Chern Numbers in Discretized Brillouin Zone: Efficient Method of Computing (Spin) Hall Conductances',
            '2005', '10.1143/JPSJ.74.1674',
            {'journal': 'J. Phys. Soc. Jpn.', 'volume': '74', 'pages': '1674--1677'}),
    'fh-z2': ('Fukui2007', 'Fukui, Takahiro and Hatsugai, Yasuhiro',
              'Quantum Spin Hall Effect in Three Dimensional Materials: Lattice Computation of Z2 Topological Invariants and Its Application to Bi and Sb',
              '2007', '10.1143/JPSJ.76.053702',
              {'journal': 'J. Phys. Soc. Jpn.', 'volume': '76', 'pages': '053702'}),
    'berry-review': ('Xiao2010', 'Xiao, Di and Chang, Ming-Che and Niu, Qian',
                     'Berry Phase Effects on Electronic Properties', '2010', '10.1103/RevModPhys.82.1959',
                     {'journal': 'Rev. Mod. Phys.', 'volume': '82', 'pages': '1959'}),
}


def utc():
    return datetime.now(timezone.utc).isoformat()


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def git(source, *args):
    if not shutil.which('git'):
        return None
    p = subprocess.run(['git', '-C', str(source), *args], capture_output=True, text=True, check=False)
    return p.stdout.strip() if p.returncode == 0 else None


def source_path(value):
    source = Path(value or os.environ.get('VASPBERRY_ROOT') or Path.cwd()).expanduser().resolve()
    for name in ('VERSION', 'CITATION.cff', 'tools/vaspberry_post.py'):
        if not (source / name).is_file():
            raise ValueError(f'Not a supported VASPBERRY source directory: missing {source / name}')
    return source


def scalar(text, key):
    """Read a top-level single-line scalar from the release's CFF (not general YAML)."""
    match = re.search(r'^' + re.escape(key) + r':\s*([^\n]+)$', text, re.MULTILINE)
    return match.group(1).strip().strip('\"\'') if match else None


def readme_citation_section(text):
    """Select the README Citation heading and its subsections, excluding later peers."""
    heading = re.search(r'^(#{1,6})[ \t]+Citation(?:s| of the code)?\b[^\n]*\n',
                        text, re.MULTILINE | re.IGNORECASE)
    if not heading:
        return None
    remainder = text[heading.end():]
    following = re.search(r'^#{1,' + str(len(heading[1])) + r'}[ \t]+',
                          remainder, re.MULTILINE)
    return remainder[:following.start()] if following else remainder


def article_entries(text):
    """Extract balanced article entries verbatim; report malformed input without guessing."""
    entries, warnings, consumed = [], [], 0
    for match in re.finditer(r'@article\b\s*([{(])', text, re.IGNORECASE):
        if match.start() < consumed:
            continue
        delimiter = match[1]
        base = 1 if delimiter == '{' else 0
        depth, quoted, escaped, comment = base, False, False, False
        end = None
        for pos in range(match.end(), len(text)):
            char = text[pos]
            if comment:
                comment = char != '\n'
                continue
            if escaped:
                escaped = False
                continue
            if char == '\\':
                escaped = True
            elif char == '%' and depth == base and not quoted:
                comment = True
            elif char == '"' and depth == base:
                quoted = not quoted
            elif char == '{':
                depth += 1
            elif char == '}':
                depth -= 1
                if depth < 0 or (delimiter == '{' and depth == 0 and quoted):
                    break
                if delimiter == '{' and depth == 0:
                    end = pos + 1
                    break
            elif char == ')' and delimiter == '(' and depth == 0 and not quoted:
                end = pos + 1
                break
        if end is None:
            warnings.append('Malformed README @article entry omitted; inspect the Citation section.')
            continue
        block = text[match.start():end]
        key = re.match(r'@article\s*[{(]\s*([^\s,{}()]+)\s*,', block, re.IGNORECASE)
        if not key or not all(re.search(r'\b' + field + r'\s*=', block, re.IGNORECASE)
                              for field in ('author', 'title', 'journal', 'year')):
            warnings.append('Incomplete README @article entry omitted; inspect the Citation section.')
            continue
        consumed = end
        entries.append((key[1], block))
    return entries, warnings


def bibtex_doi(block):
    match = re.search(r'\bdoi\s*=\s*[{\"]([^}\"]+)[}\"]', block, re.IGNORECASE)
    return match[1].strip() if match else None


def recommended_references(source, meta):
    """Preserve README article credit, plus the verified historical README PRB reference."""
    readme = source / 'README.md'
    blocks, references, warnings = [], [], []
    readme_url = (f"{REPOSITORY}/blob/{meta['commit']}/README.md#citation" if meta['commit']
                  else REPOSITORY + '#citation')
    readme_meta = {'path': 'README.md', 'sha256': None, 'citation_section_sha256': None,
                   'source_url': readme_url}
    if readme.is_file():
        readme_text = readme.read_text(encoding='utf-8')
        readme_meta['sha256'] = digest(readme)
        section = readme_citation_section(readme_text)
        if section is None:
            warnings.append('README has no Citation section; no current README article entries were inferred.')
        else:
            readme_meta['citation_section_sha256'] = hashlib.sha256(section.encode('utf-8')).hexdigest()
            entries, parse_warnings = article_entries(section)
            warnings.extend(parse_warnings)
            if not entries:
                warnings.append('No complete @article entries found in the README Citation section.')
            for key, block in entries:
                if any(item['key'] == key for item in references):
                    warnings.append(f'Duplicate README article key {key} omitted; inspect the Citation section.')
                    continue
                blocks.append(block)
                references.append({'key': key, 'doi': bibtex_doi(block), 'source_kind': 'engine_readme',
                                   'source_url': readme_url,
                                   'bibtex_sha256': hashlib.sha256(block.encode('utf-8')).hexdigest()})
    else:
        warnings.append('README.md is absent; no current README article entries were inferred.')
    if not any((item['doi'] or '').lower() == HISTORICAL_PRB_DOI.lower()
               or item['key'].lower() == 'physrevb.93.041404' for item in references):
        blocks.append(HISTORICAL_PRB_BIB)
        references.append({'key': 'PhysRevB.93.041404', 'doi': HISTORICAL_PRB_DOI,
                           'source_kind': 'historical_readme',
                           'source_url': HISTORICAL_PRB_SOURCE,
                           'metadata_url': 'https://journals.aps.org/prb/abstract/' + HISTORICAL_PRB_DOI,
                           'bibtex_sha256': hashlib.sha256(HISTORICAL_PRB_BIB.encode('utf-8')).hexdigest()})
    return blocks, references, readme_meta, warnings


def identity(source):
    version = (source / 'VERSION').read_text().strip()
    if not re.fullmatch(r'[0-9]+\.[0-9]+\.[0-9]+(?:[-+][A-Za-z0-9.-]+)?', version):
        raise ValueError('VERSION is not a supported semantic version')
    cff = (source / 'CITATION.cff').read_text(encoding='utf-8')
    if scalar(cff, 'version') != version:
        raise ValueError('CITATION.cff version differs from VERSION; reconcile metadata before citation')
    # Explicit support for this source's author record; do not guess on a new schema.
    if not re.search(r'family-names:\s*[\"\']?Kim[\"\']?\s*$', cff, re.MULTILINE) or not re.search(
            r'given-names:\s*[\"\']?Hyun-Jung[\"\']?\s*$', cff, re.MULTILINE):
        raise ValueError('CITATION.cff author format changed; read and cite its authors directly')
    git_root = git(source, 'rev-parse', '--show-toplevel')
    commit = git(source, 'rev-parse', 'HEAD') if git_root and Path(git_root).resolve() == source else None
    status = git(source, 'status', '--porcelain', '--untracked-files=no') if commit else None
    doi = scalar(cff, 'doi')
    warnings = []
    if doi and doi in ('10.5281/zenodo.1402593', 'https://doi.org/10.5281/zenodo.1402593') and version not in ('1.0', '1.0.0'):
        warnings.append('Historical VASPBERRY 1.0 DOI omitted for this later version.')
        doi = None
    if doi and not re.fullmatch(r'10\.[0-9]{4,9}/[^\s{}]+', doi):
        raise ValueError('CFF DOI must be a plain DOI identifier; inspect metadata')
    files = ('VERSION', 'CITATION.cff', 'README.md', 'vaspberry.f', 'tools/vaspberry_post.py', 'tools/wavecar_fukui.py',
             'tools/vaspberry_transport.py', 'validation/models/fukui-chern/run.py')
    return {'software': 'VASPBERRY', 'author': 'Hyun-Jung Kim', 'version': version,
            'commit': commit, 'tracked_worktree_modified': bool(status) if status is not None else None,
            'repository': REPOSITORY, 'source_url': f'{REPOSITORY}/tree/{commit}' if commit else REPOSITORY,
            'release_date': scalar(cff, 'date-released'), 'doi': doi,
            'is_target_release_commit': commit == TARGET_COMMIT,
            'source_sha256': {name: digest(source / name) for name in files if (source / name).is_file()},
            'source_hash_scope': 'Selected source files only; preserve the full source snapshot or patch, '
                                 'relevant untracked files, native executable and build record for a reproducible modified run.',
            'software_citation_policy': 'GitHub source URL and exact version/commit are sufficient; a matching DOI is optional.',
            'warnings': warnings}


def environment():
    packages = {}
    for name in ('numpy', 'matplotlib'):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    return {'python': sys.executable, 'python_version': platform.python_version(),
            'platform': platform.platform(), 'packages': packages,
            'commands': {name: shutil.which(name) for name in ('git', 'make', 'gfortran', 'ifx', 'mpiifx', 'mpifort', 'mpiexec')}}


def citations(source, output, selected):
    meta = identity(source)
    recommended_bib, recommended, readme_meta, recommended_warnings = recommended_references(source, meta)
    meta['warnings'].extend(recommended_warnings)
    report = source / 'docs/TECHNICAL_REPORT.md'
    report_text = report.read_text(encoding='utf-8') if report.exists() else ''
    selected = list(dict.fromkeys(selected))
    for key in selected:
        if METHODS[key][4] not in report_text:
            raise ValueError(f'{key} reference was not found in this checkout technical report; inspect manually')
    output.mkdir(parents=True, exist_ok=False)
    cff = (source / 'CITATION.cff').read_text(encoding='utf-8')
    if scalar(cff, 'doi') and meta['doi'] is None:
        cff = re.sub(r'^doi:[^\n]*\n?', '', cff, flags=re.MULTILINE)
    (output / 'CITATION.cff').write_text(cff, encoding='utf-8')
    key = 'KimVASPBERRY' + re.sub('[^A-Za-z0-9]', '', meta['version'])
    note = f"Source commit {meta['commit']}" if meta['commit'] else 'Source archive; Git commit unavailable, see provenance hashes'
    if meta['tracked_worktree_modified']:
        note += '; modified working tree, see provenance hashes'
    fields = {'author': 'Kim, Hyun-Jung', 'title': 'VASPBERRY', 'version': meta['version'],
              'url': meta['source_url'], 'note': note}
    if meta['release_date'] and re.fullmatch(r'[0-9]{4}-[0-9]{2}-[0-9]{2}', meta['release_date']):
        fields['year'] = meta['release_date'][:4]
    if meta['doi']:
        fields['doi'] = meta['doi']
    bib = '@software{' + key + ',\n' + ',\n'.join(f'  {k} = {{{v}}}' for k, v in fields.items()) + '\n}\n'
    credit = (f"Software reference: VASPBERRY {meta['version']} by Hyun-Jung Kim "
              f"([source]({meta['source_url']})). {note}.\n")
    credit += '\nGitHub source and exact version/commit provide software credit; a matching Zenodo DOI is optional.\n'
    if recommended:
        credit += ('\nRecommended research papers (from current and historical README instructions; '
                   'these are not universal algorithm references):\n')
        for ref in recommended:
            link = 'https://doi.org/' + ref['doi'] if ref['doi'] else ref['source_url']
            credit += f"\n- [{ref['key']}]({link}); [citation source]({ref['source_url']}).\n"
    for method in selected:
        ref, author, title, year, doi, publication = METHODS[method]
        publication_fields = ''.join(f'  {key} = {{{value}}},\n' for key, value in publication.items())
        bib += (f'\n@article{{{ref},\n  author = {{{author}}},\n  title = {{{title}}},\n'
                f'{publication_fields}  year = {{{year}}},\n  doi = {{{doi}}}\n}}\n')
        credit += f'\nSelected reference ({method}): [{title}](https://doi.org/{doi}).\n'
    credit += '\n' + '\n'.join(meta['warnings']) + '\n'
    (output / 'citations.bib').write_text(bib, encoding='utf-8')
    (output / 'recommended-references.bib').write_text('\n\n'.join(recommended_bib) + '\n', encoding='utf-8')
    (output / 'credit.md').write_text(credit, encoding='utf-8')
    meta.update(created_utc=utc(), selected_references=selected, environment=environment(),
                recommended_references=recommended, readme_citation_provenance=readme_meta,
                citation_scope='Exact software source and selected calculation methods in citations.bib; '
                               'author-recommended research papers in recommended-references.bib.')
    write_json(output / 'provenance.json', meta)
    return meta


def demo(source, output):
    meta = identity(source)
    runner = source / 'validation/models/fukui-chern/run.py'
    model_input = runner.parent / 'input.json'
    if not runner.is_file() or not model_input.is_file():
        raise ValueError('This checkout has no released Fukui model demo')
    output.mkdir(parents=True, exist_ok=False)
    command = [sys.executable, str(runner), '--output-dir', str(output / 'calculation')]
    receipt = {'schema_version': 1, 'status': 'RUNNING', 'started_utc': utc(),
               'command': command, 'cwd': str(source), 'software': meta, 'environment': environment(),
               'inputs_sha256': {'input.json': digest(model_input), 'runner': digest(runner), 'helper': digest(__file__)},
               'scope': 'Analytic QWZ model, not a VASP or real-material calculation.'}
    write_json(output / 'receipt.json', receipt)
    try:
        with (output / 'stdout.log').open('w') as out, (output / 'stderr.log').open('w') as err:
            process = subprocess.run(command, cwd=source, stdout=out, stderr=err, check=False)
        receipt['exit_code'] = process.returncode
        if process.returncode:
            raise RuntimeError(f'Engine demo exited {process.returncode}; see {output / "stderr.log"}')
        result = json.loads((output / 'calculation/result.json').read_text())
        if result.get('status') != 'PASS':
            raise RuntimeError('Engine result does not declare PASS')
        receipt['scientific_validation'] = result['numerical_checks']
        receipt['chern'] = [case['chern'] for case in result['cases']]
        citations(source, output / 'citations', ['fhs'])
        receipt['status'] = 'PASS'
    except BaseException as exc:
        receipt.update(status='FAIL', error=str(exc))
        raise
    finally:
        receipt['finished_utc'] = utc()
        receipt['outputs_sha256'] = {str(p.relative_to(output)): digest(p)
                                     for p in sorted(output.rglob('*')) if p.is_file() and p != output / 'receipt.json'}
        write_json(output / 'receipt.json', receipt)
    return {'status': 'PASS', 'output': str(output), 'chern': receipt['chern'], 'scope': receipt['scope']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    for action in ('doctor', 'cite', 'demo'):
        command = sub.add_parser(action)
        command.add_argument('--source', help='engine checkout; defaults to VASPBERRY_ROOT or working directory')
        if action != 'doctor':
            command.add_argument('--output', type=Path, required=True, help='new output directory; never overwritten')
        if action == 'cite':
            command.add_argument('--method', choices=sorted(METHODS), action='append', default=[])
    args = parser.parse_args()
    try:
        source = source_path(args.source)
        if args.action == 'doctor':
            result = {'status': 'INSPECTED', 'source': str(source), 'software': identity(source),
                      'environment': environment(),
                      'native_binaries': [str(p) for p in sorted((source / 'build').glob('vaspberry*'))
                                          if p.is_file() and os.access(p, os.X_OK)],
                      'note': 'Discovery only; dependency importability, native linkage and scientific inputs require their own checks.'}
        elif args.action == 'cite':
            result = citations(source, args.output.expanduser().resolve(), args.method)
        else:
            result = demo(source, args.output.expanduser().resolve())
        print(json.dumps(result, indent=2, allow_nan=False))
        return 0
    except (ValueError, OSError, RuntimeError) as exc:
        print(json.dumps({'status': 'FAIL', 'error': str(exc)}), file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
