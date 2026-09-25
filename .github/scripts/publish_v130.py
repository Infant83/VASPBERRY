"""Publish the checked 1.3.0 source once; never move a tag or edit a release."""
import json
import os
from pathlib import Path
import re
import subprocess
import time
import urllib.request

REPO = 'Infant83/VASPBERRY'
BASELINE = '9c92e45aca15ede6b6438304805f506ee5ddd49c'
TAG = 'v1.3.0'
ALLOWED = {
    'README.md', 'CITATION.cff', 'CHANGELOG.md',
    'docs/BUILD.md', 'docs/MIGRATION.md', 'docs/RELEASING.md',
    'docs/VALIDATION_1.3.0.md', 'docs/releases/v1.3.0.md',
    'docs/TECHNICAL_REPORT.pdf', 'tests/test_version_metadata.py',
    '.github/workflows/publish-v1.3.0.yml', '.github/scripts/publish_v130.py',
}


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


def checked_workflow(name, target, wait=False):
    deadline = time.monotonic() + (900 if wait else 0)
    while True:
        runs = api(f'actions/workflows/{name}/runs?head_sha={target}'
                   '&event=push&per_page=100')['workflow_runs']
        if runs:
            run = max(runs, key=lambda row: row['id'])
            assert run['head_sha'] == target and run['head_branch'] == 'master'
            if run['status'] == 'completed':
                assert run['conclusion'] == 'success', f'{name}: {run["conclusion"]}'
                return {'workflow': name, 'sha': target, 'url': run['html_url'],
                        'conclusion': run['conclusion'], 'attempt': run['run_attempt']}
        if not wait or time.monotonic() >= deadline:
            raise RuntimeError(f'Missing successful completed validation: {name} at {target}')
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
        assert resolved.exists() and resolved.is_relative_to(root), dest
        kind = 'tree' if resolved.is_dir() else 'blob'
        url = f'https://github.com/{REPO}/{kind}/{TAG}/{resolved.relative_to(root).as_posix()}'
        return f'[{label}]({url}{sep}{fragment})'

    text = re.sub(r'\[([^]\n]+)\]\(([^)]+)\)', link, text)
    text += f'\n\n## Published source and checks\n\nSource commit: [{target[:7]}](https://github.com/{REPO}/commit/{target}).\n\n'
    for check in checks:
        text += f'- [{check["workflow"]}: {check["sha"][:7]}]({check["url"]})\n'
    text += '\nThe release commit has its own CI; the two baseline checks validate the unchanged scientific source and actual Bi serial/MPI calculation.\n'
    return text


def main():
    target = os.environ['GITHUB_SHA']
    assert os.environ['GITHUB_REPOSITORY'] == REPO
    assert os.environ['GITHUB_REF'] == 'refs/heads/master'
    assert re.fullmatch(r'[0-9a-f]{40}', target)
    assert git('rev-parse', 'HEAD') == target
    assert git('merge-base', BASELINE, target) == BASELINE
    changed = set(git('diff', '--name-only', BASELINE, target).splitlines())
    assert changed and changed <= ALLOWED, f'Unvalidated scientific changes: {sorted(changed - ALLOWED)}'
    assert Path('VERSION').read_text().strip() == '1.3.0'
    assert '## [1.3.0] - 2026-09-25' in Path('CHANGELOG.md').read_text()
    paths = git('ls-tree', '-r', '--name-only', target).splitlines()
    assert not any(Path(p).name.upper().startswith('POTCAR') for p in paths)
    checks = [checked_workflow('ci.yml', BASELINE),
              checked_workflow('bi-z2-validation.yml', BASELINE),
              checked_workflow('ci.yml', target, wait=True)]
    refs = api(f'git/matching-refs/tags/{TAG}')
    exact = [ref for ref in refs if ref['ref'] == f'refs/tags/{TAG}']
    assert len(exact) <= 1
    if exact:
        obj = exact[0]['object']
        while obj['type'] == 'tag':
            obj = api(f'git/tags/{obj["sha"]}')['object']
        assert obj['type'] == 'commit' and obj['sha'] == target, 'Existing tag points elsewhere'
    existing = [r for r in api('releases?per_page=100') if r['tag_name'] == TAG]
    assert len(existing) <= 1
    if existing:
        release = existing[0]
        assert exact and not release['draft'] and not release['prerelease']
        result = 'ALREADY_PUBLISHED_UNCHANGED'
    else:
        release = api('releases', {
            'tag_name': TAG, 'target_commitish': target,
            'name': 'VASPBERRY 1.3.0 — Kubo transport, spin Hall and optical workflows',
            'body': release_body(Path('docs/releases/v1.3.0.md'), target, checks),
            'draft': False, 'prerelease': False, 'make_latest': 'true'})
        result = 'PUBLISHED'
    record = dict(status=result, target=target, tag=TAG, url=release['html_url'],
                  scientific_baseline=BASELINE, changed_paths=sorted(changed), checks=checks)
    Path('release-validation.json').write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps(record))
    with open(os.environ['GITHUB_STEP_SUMMARY'], 'a') as summary:
        summary.write(f'[{TAG}]({release["html_url"]}) from `{target}`: {result}.\n')


if __name__ == '__main__':
    main()
