"""Publish version 1.4.0 only after validation of this exact source commit."""
import json
import os
from pathlib import Path
import re
import subprocess
import time
import urllib.request

REPO = 'Infant83/VASPBERRY'
VERSION = '1.4.0'
TAG = 'v' + VERSION
PREVIOUS_TAG = 'v1.3.0'
PREVIOUS_COMMIT = '85150154fc9ab1f881742516bbc9faae482e9e51'
CHECKS = {
    'ci.yml': frozenset({
        'python-tests (3.10)', 'python-tests (3.12)', 'feature-examples',
        'fortran-z2-regression',
        'fortran-portability (ubuntu-22.04, gfortran-11)',
        'fortran-portability (ubuntu-24.04, gfortran-13)',
        'intel-serial-portability (intel, 2025.0, ifx, ifx, build/vaspberry-ifx)',
        'intel-serial-portability (intel-classic, 2021.10, ifort, ifort, build/vaspberry-ifort)',
    }),
    'bi-z2-validation.yml': frozenset({'validate'}),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def api(path, payload=None):
    request = urllib.request.Request(
        f'https://api.github.com/repos/{REPO}/{path}',
        data=None if payload is None else json.dumps(payload).encode(),
        headers={'Authorization': 'Bearer ' + os.environ['GH_TOKEN'],
                 'Accept': 'application/vnd.github+json',
                 'Content-Type': 'application/json',
                 'X-GitHub-Api-Version': '2022-11-28'})
    with urllib.request.urlopen(request, timeout=45) as response:
        return json.load(response)


def git(*args):
    return subprocess.check_output(['git', *args], text=True).strip()


def commit_of(obj):
    seen = set()
    while obj['type'] == 'tag':
        require(obj['sha'] not in seen and len(seen) < 16,
                'Cyclic or excessively nested tag object')
        seen.add(obj['sha'])
        obj = api(f'git/tags/{obj["sha"]}')['object']
    require(obj['type'] == 'commit', 'Tag does not resolve to a commit')
    return obj['sha']


def exact_tag(tag):
    refs = api(f'git/matching-refs/tags/{tag}')
    matches = [ref for ref in refs if ref['ref'] == f'refs/tags/{tag}']
    require(len(matches) <= 1, 'Duplicate exact tag')
    return commit_of(matches[0]['object']) if matches else None


def validate_run(run, target, workflow):
    require(run['head_sha'] == target and run['head_branch'] == 'master'
            and run['event'] == 'push', 'Validation run has a different source')
    require(run['status'] == 'completed' and run['conclusion'] == 'success',
            'Validation run is not successful')
    require(run['path'] == f'.github/workflows/{workflow}'
            and type(run['run_attempt']) is int and run['run_attempt'] >= 1,
            'Unexpected workflow or run attempt')
    data = api(f'actions/runs/{run["id"]}/attempts/{run["run_attempt"]}/jobs?per_page=100')
    jobs = data['jobs']
    require(len(jobs) == data['total_count'] and len(jobs) == len(CHECKS[workflow])
            and {job['name'] for job in jobs} == CHECKS[workflow],
            'Incomplete validation jobs')
    require(all(job['status'] == 'completed' and job['conclusion'] == 'success'
                for job in jobs), 'A validation job failed or was skipped')
    require(all(job['head_sha'] == target and job['run_id'] == run['id']
                and job['run_attempt'] == run['run_attempt'] for job in jobs),
            'Validation jobs do not belong to the checked source and attempt')
    return {'workflow': run['path'], 'sha': target, 'url': run['html_url'],
            'conclusion': run['conclusion'], 'attempt': run['run_attempt'],
            'jobs': [{'name': j['name'], 'url': j['html_url'], 'conclusion': j['conclusion']}
                     for j in jobs]}


def checked_workflow(name, target, timeout=2400):
    deadline = time.monotonic() + timeout
    while True:
        runs = api(f'actions/workflows/{name}/runs?head_sha={target}'
                   '&event=push&per_page=100')['workflow_runs']
        if runs:
            run = max(runs, key=lambda row: row['id'])
            require(run['head_sha'] == target and run['head_branch'] == 'master',
                    'Unexpected workflow source')
            if run['status'] == 'completed':
                return validate_run(run, target, name)
        if time.monotonic() >= deadline:
            raise RuntimeError(f'Missing completed validation: {name} at {target}')
        time.sleep(20)


def release_body(path, target, checks):
    root = Path.cwd().resolve()
    text = path.read_text(encoding='utf-8')

    def link(match):
        label, dest = match.groups()
        if re.match(r'^[a-z]+:', dest) or dest.startswith('#'):
            return match.group(0)
        local, sep, fragment = dest.partition('#')
        resolved = (path.parent / local).resolve()
        require(resolved.exists() and resolved.is_relative_to(root),
                f'Invalid release link: {dest}')
        kind = 'tree' if resolved.is_dir() else 'blob'
        url = f'https://github.com/{REPO}/{kind}/{TAG}/{resolved.relative_to(root).as_posix()}'
        return f'[{label}]({url}{sep}{fragment})'

    text = re.sub(r'\[([^]\n]+)\]\(([^)]+)\)', link, text)
    text += f'\n\n## Published source and checks\n\nSource commit: [{target[:7]}](https://github.com/{REPO}/commit/{target}).\n\n'
    for check in checks:
        text += f'- [{check["workflow"]}: {len(check["jobs"])} jobs passed]({check["url"]})\n'
    text += '\nBoth workflows checked this exact release commit. Historical material-convergence studies retain their original provenance.\n'
    return text


def validate_metadata():
    require(Path('VERSION').read_text().strip() == VERSION, 'Incorrect source version')
    require(f'## [{VERSION}] - 2026-09-25' in Path('CHANGELOG.md').read_text(),
            'Missing versioned changelog')
    citation = Path('CITATION.cff').read_text()
    require(re.search(rf'(?m)^version:\s*{re.escape(VERSION)}\s*$', citation)
            and f'/releases/tag/{TAG}' in citation, 'Citation version mismatch')
    for source in ('vaspberry.f', 'vaspberry_gfortran_serial.f'):
        text = Path(source).read_text()
        require(f'Version {VERSION}' in text and f'VASPBERRY (Ver {VERSION})' in text,
                'Fortran version mismatch')
    for path in ('docs/VALIDATION_1.4.0.md', 'docs/releases/v1.4.0.md',
                 'examples/features/procar-character/README.md'):
        require(Path(path).is_file(), f'Missing release documentation: {path}')


def publish(target, checks):
    require(api('git/ref/heads/master')['object']['sha'] == target,
            'Master advanced while validation was running')
    require(exact_tag(PREVIOUS_TAG) == PREVIOUS_COMMIT, 'Previous release tag changed')
    existing_tag = exact_tag(TAG)
    require(existing_tag in (None, target), 'Existing release tag points elsewhere')
    releases = [r for r in api('releases?per_page=100') if r['tag_name'] == TAG]
    require(len(releases) <= 1, 'Duplicate release')
    if releases:
        release = releases[0]
        require(existing_tag == target and not release['draft'] and not release['prerelease'],
                'Existing release is not the expected published source')
        return 'ALREADY_PUBLISHED_UNCHANGED', release
    payload = {
        'tag_name': TAG, 'target_commitish': target,
        'name': 'VASPBERRY 1.4.0 — native workflows and projected transport',
        'body': release_body(Path('docs/releases/v1.4.0.md'), target, checks),
        'draft': False, 'prerelease': False, 'make_latest': 'true'}
    require(api('git/ref/heads/master')['object']['sha'] == target,
            'Master advanced before publication')
    release = api('releases', payload)
    require(exact_tag(TAG) == target, 'Published tag does not match the checked commit')
    require(release['tag_name'] == TAG and not release['draft'] and not release['prerelease'],
            'Publication did not produce the requested final release')
    return 'PUBLISHED', release


def main():
    target = os.environ['GITHUB_SHA']
    require(os.environ['GITHUB_REPOSITORY'] == REPO and
            os.environ['GITHUB_REF'] == 'refs/heads/master', 'Unexpected repository/ref')
    require(re.fullmatch(r'[0-9a-f]{40}', target) and git('rev-parse', 'HEAD') == target,
            'Checkout does not match publication target')
    require(git('merge-base', PREVIOUS_COMMIT, target) == PREVIOUS_COMMIT,
            'Release is not descended from the previous source')
    validate_metadata()
    paths = git('ls-tree', '-r', '--name-only', target).splitlines()
    require(not any(Path(p).name.upper().startswith('POTCAR') for p in paths),
            'Licensed input accidentally included')
    checks = [checked_workflow(name, target) for name in CHECKS]
    status, release = publish(target, checks)
    record = dict(status=status, target=target, tag=TAG, url=release['html_url'], checks=checks)
    Path('release-validation.json').write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps(record))
    with open(os.environ['GITHUB_STEP_SUMMARY'], 'a') as summary:
        summary.write(f'[{TAG}]({release["html_url"]}) from `{target}`: {status}.\n')


if __name__ == '__main__':
    main()
