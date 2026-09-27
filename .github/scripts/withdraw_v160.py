"""Withdraw only the duplicate 1.6.0 release/tag at the maintainer's request."""
import json
import os
from pathlib import Path
import subprocess
import urllib.error
import urllib.request

REPO = 'Infant83/VASPBERRY'
KEEP_TAG = 'v1.5.1'
KEEP_COMMIT = 'a9b2aa15f9b1d89d7da1fc10889a62bbf026059a'
KEEP_ID = 397441924
WITHDRAW_TAG = 'v1.6.0'
WITHDRAW_COMMIT = 'eb0b468aa0674f5ec579e75a014fe4591ec4cdfe'
WITHDRAW_ID = 397438702
OLD_INTRO = ('This patch release supersedes the earlier 1.6.0 numbering with the same\n'
             'automatic-spinor functionality. Use **1.5.1** for new installations. The\n'
             'published 1.6.0 tag and historical reference results remain unchanged.')
NEW_INTRO = ('This patch release contains the automatic-spinor update. The duplicate 1.6.0\n'
             'release and tag have been withdrawn. Use **1.5.1** for new installations;\n'
             'its source tag and scientific results are unchanged.')


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def api(path, payload=None, method=None, missing_ok=False):
    request = urllib.request.Request(
        f'https://api.github.com/repos/{REPO}/{path}',
        data=None if payload is None else json.dumps(payload).encode(), method=method,
        headers={'Authorization': 'Bearer ' + os.environ['GH_TOKEN'],
                 'Accept': 'application/vnd.github+json',
                 'Content-Type': 'application/json', 'X-GitHub-Api-Version': '2022-11-28'})
    try:
        with urllib.request.urlopen(request, timeout=45) as response:
            data = response.read()
            return json.loads(data) if data else None
    except urllib.error.HTTPError as exc:
        if missing_ok and exc.code == 404:
            return None
        raise


def exact_tag(tag):
    refs = api(f'git/matching-refs/tags/{tag}')
    refs = [r for r in refs if r['ref'] == f'refs/tags/{tag}']
    require(len(refs) <= 1, 'Duplicate exact tag')
    if not refs:
        return None
    obj = refs[0]['object']
    require(obj['type'] == 'commit', 'Expected the original lightweight tag')
    return obj['sha']


def verify_keep():
    require(exact_tag(KEEP_TAG) == KEEP_COMMIT, 'Retained source tag changed')
    release = api(f'releases/tags/{KEEP_TAG}')
    require(release['id'] == KEEP_ID and release['tag_name'] == KEEP_TAG
            and release['target_commitish'] == KEEP_COMMIT
            and not release['draft'] and not release['prerelease'],
            'Unexpected retained release')
    require(api('releases/latest')['id'] == KEEP_ID, 'Retained release is not latest')
    return release


def withdraw(target, record_path=Path('withdrawal-validation.json')):
    require(api('git/ref/heads/master')['object']['sha'] == target, 'Master advanced')
    kept = verify_keep()
    old_tag = exact_tag(WITHDRAW_TAG)
    require(old_tag in (None, WITHDRAW_COMMIT), 'Withdrawal tag points elsewhere')
    old = api(f'releases/tags/{WITHDRAW_TAG}', missing_ok=True)
    if old is not None:
        require(old['id'] == WITHDRAW_ID and old['tag_name'] == WITHDRAW_TAG
                and old['target_commitish'] == WITHDRAW_COMMIT
                and not old['draft'] and not old['prerelease'] and not old['assets'],
                'Unexpected withdrawal release or assets')
    body = kept.get('body') or ''
    require((body.count(OLD_INTRO) == 1 and NEW_INTRO not in body)
            or (body.count(NEW_INTRO) == 1 and OLD_INTRO not in body),
            'Retained release intro changed; review before modifying it')
    corrected_body = body.replace(OLD_INTRO, NEW_INTRO)
    record = {'status': 'PREFLIGHT_PASS', 'target': target,
              'keep_before': kept, 'withdraw_before': old, 'withdraw_tag_before': old_tag}
    record_path.write_text(json.dumps(record, indent=2) + '\n')
    require(api('git/ref/heads/master')['object']['sha'] == target, 'Master advanced')
    if old is not None:
        require(api(f'releases/tags/{WITHDRAW_TAG}') == old,
                'Withdrawal release changed after review')
        api(f'releases/{WITHDRAW_ID}', method='DELETE')
    # A retry after partial completion is safe; unrelated refs/releases never match.
    if old_tag is not None:
        require(exact_tag(WITHDRAW_TAG) == WITHDRAW_COMMIT, 'Withdrawal tag changed')
        api(f'git/refs/tags/{WITHDRAW_TAG}', method='DELETE')
    require(api(f'releases/tags/{WITHDRAW_TAG}', missing_ok=True) is None,
            'Withdrawn release still exists')
    require(exact_tag(WITHDRAW_TAG) is None, 'Withdrawn tag still exists')
    if corrected_body != body:
        require(verify_keep().get('body') == body, 'Retained release changed after review')
        api(f'releases/{KEEP_ID}', {'body': corrected_body, 'make_latest': 'true'},
            method='PATCH')
    kept_after = verify_keep()
    require(kept_after.get('body') == corrected_body, 'Corrected release note not saved')
    require(kept_after['assets'] == kept['assets']
            and kept_after['target_commitish'] == kept['target_commitish'],
            'Retained release source or assets changed')
    record.update(status='WITHDRAWN', keep_after=kept_after,
                  withdrawn_release_id=WITHDRAW_ID, withdrawn_tag=WITHDRAW_TAG,
                  retained_tag=KEEP_TAG, retained_commit=KEEP_COMMIT)
    record_path.write_text(json.dumps(record, indent=2) + '\n')
    return record


def main():
    require(os.environ['GITHUB_REPOSITORY'] == REPO
            and os.environ['GITHUB_REF'] == 'refs/heads/master', 'Unexpected repository/ref')
    target = os.environ['GITHUB_SHA']
    require(subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == target,
            'Checkout does not match administrative commit')
    require(Path('VERSION').read_text().strip() == '1.5.1', 'Current version changed')
    result = withdraw(target)
    with open(os.environ['GITHUB_STEP_SUMMARY'], 'a') as stream:
        stream.write('Withdrawn duplicate v1.6.0 release and tag; v1.5.1 remains latest.\n')
    print(json.dumps({'status': result['status'], 'retained_commit': KEEP_COMMIT}))


if __name__ == '__main__':
    main()
