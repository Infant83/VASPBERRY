"""The explicitly requested withdrawal must never change the retained source."""
import copy
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT/'.github/scripts'/f'{name}.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


w = load('withdraw_v160')
TARGET = 'a' * 40


class WithdrawalTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.record = Path(self.temp.name)/'record.json'
        self.state = {
            'master': TARGET, 'latest': w.KEEP_ID,
            'tags': {w.KEEP_TAG: w.KEEP_COMMIT, w.WITHDRAW_TAG: w.WITHDRAW_COMMIT},
            'keep': dict(id=w.KEEP_ID, tag_name=w.KEEP_TAG, draft=False, prerelease=False,
                         target_commitish=w.KEEP_COMMIT, assets=[],
                         body='Title\n\n'+w.OLD_INTRO+'\n\nOriginal scientific notes and CI links.'),
            'old': dict(id=w.WITHDRAW_ID, tag_name=w.WITHDRAW_TAG, draft=False,
                        prerelease=False, target_commitish=w.WITHDRAW_COMMIT, assets=[],
                        body='Original withdrawn notes'),
            'mutations': [],
        }

    def api(self, path, payload=None, method=None, missing_ok=False):
        s = self.state
        if method in ('DELETE', 'PATCH'):
            self.assertTrue(self.record.exists(), 'Pre-mutation record must exist')
            s['mutations'].append((method, path, copy.deepcopy(payload)))
            if method == 'DELETE':
                self.assertIsNone(payload)
                if path == f'releases/{w.WITHDRAW_ID}': s['old'] = None
                elif path == f'git/refs/tags/{w.WITHDRAW_TAG}': s['tags'].pop(w.WITHDRAW_TAG)
                else: self.fail('Unexpected deletion '+path)
                return None
            self.assertEqual(path, f'releases/{w.KEEP_ID}')
            self.assertEqual(set(payload), {'body', 'make_latest'})
            self.assertEqual(payload['make_latest'], 'true')
            s['keep']['body'] = payload['body']
            return copy.deepcopy(s['keep'])
        self.assertIsNone(method); self.assertIsNone(payload)
        if path == 'git/ref/heads/master': return {'object': {'sha': s['master']}}
        if path == 'releases/latest': return {'id': s['latest']}
        if path == f'releases/tags/{w.KEEP_TAG}': return copy.deepcopy(s['keep'])
        if path == f'releases/tags/{w.WITHDRAW_TAG}':
            if s['old'] is None: self.assertTrue(missing_ok)
            return copy.deepcopy(s['old'])
        if path.startswith('git/matching-refs/tags/'):
            tag = path.removeprefix('git/matching-refs/tags/')
            if tag not in s['tags']: return []
            return [{'ref': 'refs/tags/'+tag,
                     'object': {'type': 'commit', 'sha': s['tags'][tag]}}]
        self.fail('Unexpected API operation '+path)

    def call(self, api=None):
        with mock.patch.object(w, 'api', side_effect=api or self.api):
            return w.withdraw(TARGET, self.record)

    def test_deletes_only_duplicate_and_preserves_retained_source_and_notes(self):
        before = copy.deepcopy(self.state['keep'])
        result = self.call()
        self.assertEqual(result['status'], 'WITHDRAWN')
        self.assertIsNone(self.state['old'])
        self.assertEqual(self.state['tags'], {w.KEEP_TAG: w.KEEP_COMMIT})
        self.assertEqual(self.state['keep'], dict(before, body=before['body'].replace(w.OLD_INTRO,w.NEW_INTRO)))
        self.assertEqual([(m,p) for m,p,_ in self.state['mutations']],
                         [('DELETE',f'releases/{w.WITHDRAW_ID}'),
                          ('DELETE',f'git/refs/tags/{w.WITHDRAW_TAG}'),
                          ('PATCH',f'releases/{w.KEEP_ID}')])

    def test_completed_run_is_idempotent(self):
        self.call(); before = copy.deepcopy(self.state)
        self.call()
        self.assertEqual(self.state, before)

    def test_retry_after_release_deletion_only_finishes_remaining_actions(self):
        self.state['old'] = None
        self.call()
        self.assertEqual(len(self.state['mutations']), 2)

    def test_retry_after_both_deletions_only_corrects_note(self):
        self.state['old'] = None; del self.state['tags'][w.WITHDRAW_TAG]
        self.call()
        self.assertEqual(len(self.state['mutations']), 1)
        self.assertEqual(self.state['mutations'][0][0], 'PATCH')

    def test_changed_master_or_latest_prevents_mutation(self):
        for key, value in [('master','b'*40),('latest',0)]:
            with self.subTest(key=key):
                self.setUp(); self.state[key] = value
                with self.assertRaises(RuntimeError): self.call()
                self.assertFalse(self.state['mutations'])

    def test_changed_tag_prevents_mutation(self):
        for tag in [w.KEEP_TAG,w.WITHDRAW_TAG]:
            with self.subTest(tag=tag):
                self.setUp(); self.state['tags'][tag] = 'b'*40
                with self.assertRaises(RuntimeError): self.call()
                self.assertFalse(self.state['mutations'])

    def test_unexpected_release_or_assets_prevents_mutation(self):
        for release,key,value in [('keep','id',0),('keep','draft',True),('keep','prerelease',True),
                                   ('keep','target_commitish','b'*40),
                                   ('old','id',0),('old','tag_name','v1.6.1'),
                                   ('old','target_commitish','b'*40),('old','assets',[{'id':1}]),
                                   ('old','draft',True),('old','prerelease',True)]:
            with self.subTest(release=release,key=key):
                self.setUp(); self.state[release][key] = value
                with self.assertRaises(RuntimeError): self.call()
                self.assertFalse(self.state['mutations'])

    def test_changed_intro_prevents_any_deletion(self):
        self.state['keep']['body'] = 'Unexpected manually edited introduction'
        with self.assertRaises(RuntimeError): self.call()
        self.assertFalse(self.state['mutations'])

    def test_racing_old_release_edit_is_not_deleted(self):
        count = 0
        def changed(path, **kwargs):
            nonlocal count
            if path == f'releases/tags/{w.WITHDRAW_TAG}':
                count += 1
                if count == 2: self.state['old']['body'] = 'New edit'
            return self.api(path, **kwargs)
        with self.assertRaisesRegex(RuntimeError, 'changed after review'): self.call(changed)
        self.assertFalse(self.state['mutations'])

    def test_failed_deletion_is_reported_without_patching_retained_release(self):
        def ignored(path, **kwargs):
            if path == f'releases/{w.WITHDRAW_ID}' and kwargs.get('method') == 'DELETE': return None
            return self.api(path, **kwargs)
        with self.assertRaisesRegex(RuntimeError, 'still exists'): self.call(ignored)
        self.assertNotIn('PATCH', [m for m,_,_ in self.state['mutations']])

    def test_tag_prefix_collision_is_not_deleted(self):
        refs = [{'ref':'refs/tags/v1.6.0-rc1','object':{'type':'commit','sha':'b'*40}}]
        with mock.patch.object(w,'api',return_value=refs):
            self.assertIsNone(w.exact_tag(w.WITHDRAW_TAG))

    def test_old_publishers_are_retired_before_any_api_call(self):
        for name in ['publish_v151','publish_v160']:
            publisher = load(name)
            with mock.patch.object(publisher,'api') as api:
                with self.assertRaisesRegex(RuntimeError,'Retired'): publisher.main()
                api.assert_not_called()
        self.assertFalse((ROOT/'.github/workflows/publish-v1.5.1.yml').exists())
        self.assertFalse((ROOT/'.github/workflows/publish-v1.6.0.yml').exists())


if __name__ == '__main__':
    unittest.main()
