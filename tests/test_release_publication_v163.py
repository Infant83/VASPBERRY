"""Publication must consume exact-source CI and never move or rewrite releases."""
import copy
import importlib.util
import json
import os
from pathlib import Path
import re
import tempfile
import unittest
import urllib.error
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('publish_v163', ROOT/'.github/scripts/publish_v163.py')
publisher = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(publisher)
TARGET = 'a'*40


class ReleaseChecks(unittest.TestCase):
    def test_automatic_trigger_is_limited_to_this_workflow_and_receipts_survive_failure(self):
        text = (ROOT/'.github/workflows/publish-v1.6.3.yml').read_text()
        self.assertIn('  push:', text)
        self.assertIn("    paths: ['.github/workflows/publish-v1.6.3.yml']", text)
        self.assertEqual(text.count('paths:'), 1)
        self.assertIn('  workflow_dispatch:', text)
        self.assertIn('  cancel-in-progress: false', text)
        self.assertIn('        if: always()', text)
        self.assertIn('          path: release-validation.json', text)

    def test_completed_publishers_have_no_automatic_trigger(self):
        for version in ('1.4.2', '1.6.1'):
            with self.subTest(version=version):
                workflow = (ROOT/f'.github/workflows/publish-v{version}.yml').read_text()
                self.assertIn('  workflow_dispatch:', workflow)
                self.assertNotIn('\n  push:', workflow)

    def test_archive_installation_jobs_cannot_be_omitted_or_skipped(self):
        workflow = 'installation-validation.yml'
        self.assertEqual(publisher.CHECKS[workflow], frozenset({
            'Archive install (macOS ARM64, GNU/OpenMPI)',
            'Archive install (Linux, GNU/MPICH)',
        }))
        for variant in ('missing', 'skipped', 'wrong_source'):
            run, jobs = self.fixture(workflow)
            if variant == 'missing':
                jobs['jobs'].pop(); jobs['total_count'] -= 1
            elif variant == 'skipped':
                jobs['jobs'][0]['conclusion'] = 'skipped'
            else:
                jobs['jobs'][0]['head_sha'] = 'b' * 40
            with mock.patch.object(publisher, 'api', return_value=jobs):
                with self.assertRaises(RuntimeError):
                    publisher.validate_run(run, TARGET, workflow)

    def test_intel_mpi_both_compilers_are_required(self):
        self.assertEqual(publisher.CHECKS['intel-mpi-validation.yml'],
                         frozenset({'Intel MPI (ifx)', 'Intel MPI (ifort)'}))
        self.assertEqual(sum(map(len, publisher.CHECKS.values())), 13)
        for conclusion in ('skipped', 'failure'):
            run, jobs = self.fixture('intel-mpi-validation.yml')
            jobs['jobs'][0]['conclusion'] = conclusion
            with mock.patch.object(publisher, 'api', return_value=jobs):
                with self.assertRaises(RuntimeError):
                    publisher.validate_run(run, TARGET, 'intel-mpi-validation.yml')

    def fixture(self, workflow='ci.yml'):
        run = dict(head_sha=TARGET, head_branch='master', event='push', status='completed',
                   conclusion='success', id=123, run_attempt=2, path='.github/workflows/'+workflow,
                   html_url='https://example.test/run/123')
        jobs = [dict(name=name, head_sha=TARGET, run_id=123, run_attempt=2, status='completed',
                     conclusion='success', html_url='https://example.test/job/'+str(i))
                for i,name in enumerate(sorted(publisher.CHECKS[workflow]))]
        return run, dict(jobs=jobs,total_count=len(jobs))

    def test_success_uses_exact_named_jobs_and_attempt_endpoint(self):
        for workflow in publisher.CHECKS:
            with self.subTest(workflow=workflow):
                run, jobs = self.fixture(workflow)
                with mock.patch.object(publisher,'api',return_value=jobs) as api:
                    result = publisher.validate_run(run,TARGET,workflow)
                api.assert_called_once_with('actions/runs/123/attempts/2/jobs?per_page=100')
                self.assertEqual(result['sha'],TARGET)
                self.assertEqual(result['attempt'],2)
                self.assertEqual(len(result['jobs']),len(publisher.CHECKS[workflow]))

    def test_workflow_wrong_source_or_failed_conclusion_fails_before_jobs(self):
        for key,value in [('head_sha','b'*40),('head_branch','feature/test'),('event','pull_request'),
                          ('status','in_progress'),('conclusion','failure'),('conclusion','cancelled'),
                          ('path','.github/workflows/unrelated.yml'),('run_attempt',0)]:
            with self.subTest(key=key,value=value):
                run,jobs = self.fixture(); run[key]=value
                with mock.patch.object(publisher,'api',return_value=jobs) as api:
                    with self.assertRaises(RuntimeError): publisher.validate_run(run,TARGET,'ci.yml')
                    api.assert_not_called()

    def test_successful_workflow_with_skipped_or_failed_job_is_rejected(self):
        for conclusion in ('failure','cancelled','skipped','neutral',None):
            with self.subTest(conclusion=conclusion):
                run,jobs = self.fixture(); jobs['jobs'][0]['conclusion']=conclusion
                with mock.patch.object(publisher,'api',return_value=jobs):
                    with self.assertRaisesRegex(RuntimeError,'failed or was skipped'):
                        publisher.validate_run(run,TARGET,'ci.yml')

    def test_incomplete_duplicate_or_replaced_job_set_is_rejected(self):
        for variant in ('missing','duplicate','unrelated','truncated_page'):
            with self.subTest(variant=variant):
                run,jobs = self.fixture()
                if variant == 'missing': jobs['jobs'].pop(); jobs['total_count']-=1
                if variant == 'duplicate': jobs['jobs'][0]['name']=jobs['jobs'][1]['name']
                if variant == 'unrelated': jobs['jobs'][0]['name']='unrelated successful job'
                if variant == 'truncated_page': jobs['total_count']+=1
                with mock.patch.object(publisher,'api',return_value=jobs):
                    with self.assertRaisesRegex(RuntimeError,'Incomplete validation jobs'):
                        publisher.validate_run(run,TARGET,'ci.yml')

    def test_jobs_from_wrong_source_run_or_attempt_are_rejected(self):
        for key,value in [('head_sha','b'*40),('run_id',124),('run_attempt',1)]:
            with self.subTest(key=key):
                run,jobs=self.fixture(); jobs['jobs'][0][key]=value
                with mock.patch.object(publisher,'api',return_value=jobs):
                    with self.assertRaisesRegex(RuntimeError,'checked source and attempt'):
                        publisher.validate_run(run,TARGET,'ci.yml')

    def test_newest_failed_run_does_not_fall_back_to_old_success(self):
        good,_ = self.fixture()
        bad = dict(good,id=124,conclusion='failure')
        with mock.patch.object(publisher,'api',return_value={'workflow_runs':[good,bad]}), \
                mock.patch.object(publisher.time,'sleep') as sleep:
            with self.assertRaisesRegex(RuntimeError,'not successful'):
                publisher.checked_workflow('ci.yml',TARGET,timeout=0)
            sleep.assert_not_called()

    def test_missing_workflow_times_out_without_publication(self):
        with mock.patch.object(publisher,'api',return_value={'workflow_runs':[]}), \
                mock.patch.object(publisher.time,'sleep') as sleep:
            with self.assertRaisesRegex(RuntimeError,'Missing completed validation'):
                publisher.checked_workflow('ci.yml',TARGET,timeout=0)
            sleep.assert_not_called()

    def test_newer_running_attempt_blocks_final_revalidation(self):
        run, _ = self.fixture()
        newer = dict(run, id=124, run_attempt=3, status='in_progress', conclusion=None)
        with mock.patch.object(publisher, 'api', return_value={'workflow_runs': [run, newer]}):
            with self.assertRaisesRegex(RuntimeError, 'Missing completed validation'):
                publisher.checked_workflow('ci.yml', TARGET, timeout=0)


class ImmutablePublication(unittest.TestCase):
    def fixture(self, *, master=TARGET, previous=None, current=None, existing=False):
        state = dict(master=master,previous=publisher.PREVIOUS_COMMIT if previous is None else previous,
                     current=current,posts=[],master_reads=0)
        release = dict(tag_name=publisher.TAG,draft=False,prerelease=False,
                       html_url='https://example.test/release')
        def tag(name):
            return state['previous'] if name == publisher.PREVIOUS_TAG else state['current']
        def api(path,payload=None):
            if payload is not None:
                self.assertEqual(path,'releases')
                state['posts'].append(copy.deepcopy(payload)); state['current']=payload['target_commitish']
                return dict(release)
            if path == 'git/ref/heads/master':
                state['master_reads']+=1
                value = state['master']
                if isinstance(value,list): value=value[min(state['master_reads']-1,len(value)-1)]
                return {'object':{'sha':value}}
            if path == f'releases/tags/{publisher.TAG}':
                if existing or state['posts']:
                    return release
                raise urllib.error.HTTPError('https://example.test/release', 404, 'missing', {}, None)
            raise AssertionError('Unexpected API read '+path)
        return state,release,tag,api

    def call(self,state,tag,api):
        with mock.patch.object(publisher,'exact_tag',side_effect=tag), \
                mock.patch.object(publisher,'api',side_effect=api), \
                mock.patch.object(publisher,'release_body',return_value='reviewed release notes'):
            return publisher.publish(TARGET,[])

    def test_first_publication_pins_exact_target_and_verifies_new_tag(self):
        state,release,tag,api=self.fixture()
        status,actual=self.call(state,tag,api)
        self.assertEqual(status,'PUBLISHED'); self.assertEqual(actual,release)
        self.assertEqual(len(state['posts']),1)
        self.assertEqual(state['posts'][0]['target_commitish'],TARGET)
        self.assertFalse(state['posts'][0]['draft']); self.assertFalse(state['posts'][0]['prerelease'])
        self.assertEqual(state['master_reads'],2)

    def test_existing_final_release_is_returned_without_mutation(self):
        state,release,tag,api=self.fixture(current=TARGET,existing=True)
        self.assertEqual(self.call(state,tag,api),('ALREADY_PUBLISHED_UNCHANGED',release))
        self.assertEqual(state['posts'],[])

    def test_repeated_publication_creates_only_once(self):
        state, release, tag, api = self.fixture()
        self.assertEqual(self.call(state, tag, api)[0], 'PUBLISHED')
        self.assertEqual(self.call(state, tag, api), ('ALREADY_PUBLISHED_UNCHANGED', release))
        self.assertEqual(len(state['posts']), 1)

    def test_only_direct_release_lookup_404_means_absent(self):
        for code in (404, 401, 403, 500):
            with self.subTest(code=code):
                error = urllib.error.HTTPError('https://example.test/release', code, 'error', {}, None)
                with mock.patch.object(publisher, 'api', side_effect=error) as api:
                    if code == 404:
                        self.assertIsNone(publisher.existing_release())
                    else:
                        with self.assertRaises(urllib.error.HTTPError):
                            publisher.existing_release()
                    api.assert_called_once_with(f'releases/tags/{publisher.TAG}')

    def test_tag_changes_during_release_body_preparation_block_post(self):
        for name in (publisher.PREVIOUS_TAG, publisher.TAG):
            with self.subTest(tag=name):
                state, _, tag, api = self.fixture()
                reads = 0
                def changing_tag(query):
                    nonlocal reads
                    if query == name:
                        reads += 1
                        if reads == 2:
                            return 'b' * 40
                    return tag(query)
                with self.assertRaises(RuntimeError):
                    self.call(state, changing_tag, api)
                self.assertEqual(state['posts'], [])

    def test_concurrent_creation_or_lost_post_response_is_reconciled_without_retry(self):
        for error in (TimeoutError('lost response'), urllib.error.HTTPError(
                'https://example.test/release', 422, 'already exists', {}, None)):
            with self.subTest(error=type(error).__name__):
                state, release, tag, api = self.fixture()
                def ambiguous_api(path, payload=None):
                    result = api(path, payload)
                    if payload is not None:
                        raise error
                    return result
                self.assertEqual(self.call(state, tag, ambiguous_api),
                                 ('ALREADY_PUBLISHED_UNCHANGED', release))
                self.assertEqual(len(state['posts']), 1)

    def test_ambiguous_post_does_not_certify_missing_conflicting_or_draft_release(self):
        for variant in ('missing', 'conflicting_tag', 'draft', 'prerelease', 'wrong_name'):
            with self.subTest(variant=variant):
                state, release, tag, api = self.fixture()
                posts = 0
                def ambiguous_api(path, payload=None):
                    nonlocal posts
                    if payload is None:
                        return api(path)
                    posts += 1
                    if variant != 'missing':
                        api(path, payload)
                    if variant == 'conflicting_tag': state['current'] = 'b' * 40
                    if variant in ('draft', 'prerelease'): release[variant] = True
                    if variant == 'wrong_name': release['tag_name'] = 'v9.9.9'
                    raise TimeoutError('lost response')
                with self.assertRaises(TimeoutError):
                    self.call(state, tag, ambiguous_api)
                self.assertEqual(posts, 1)

    def test_changed_master_old_tag_or_new_tag_blocks_post(self):
        for overrides in [dict(master='b'*40),dict(master=[TARGET,'b'*40]),
                          dict(previous='b'*40),dict(current='b'*40)]:
            with self.subTest(overrides=overrides):
                state,_,tag,api=self.fixture(**overrides)
                with self.assertRaises(RuntimeError): self.call(state,tag,api)
                self.assertEqual(state['posts'],[])

    def test_existing_draft_prerelease_or_missing_tag_is_not_modified(self):
        for variant in ('draft','prerelease','missing_tag'):
            with self.subTest(variant=variant):
                state,release,tag,api=self.fixture(current=TARGET,existing=True)
                if variant == 'missing_tag': state['current']=None
                else: release[variant]=True
                with self.assertRaises(RuntimeError): self.call(state,tag,api)
                self.assertEqual(state['posts'],[])

    def test_annotation_resolution_and_prefix_collision(self):
        refs=[{'ref':'refs/tags/'+publisher.TAG+'-rc1','object':{'type':'commit','sha':'b'*40}},
              {'ref':'refs/tags/'+publisher.TAG,'object':{'type':'tag','sha':'c'*40}}]
        with mock.patch.object(publisher,'api',side_effect=[refs,{'object':{'type':'commit','sha':TARGET}}]):
            self.assertEqual(publisher.exact_tag(publisher.TAG),TARGET)
        with mock.patch.object(publisher,'api',return_value=refs[:1]):
            self.assertIsNone(publisher.exact_tag(publisher.TAG))

    def test_cyclic_or_noncommit_tag_is_rejected(self):
        obj={'type':'tag','sha':'c'*40}
        with mock.patch.object(publisher,'api',return_value={'object':obj}):
            with self.assertRaisesRegex(RuntimeError,'Cyclic'): publisher.commit_of(obj)
        with self.assertRaisesRegex(RuntimeError,'does not resolve'):
            publisher.commit_of({'type':'tree','sha':TARGET})

    def test_post_creation_tag_race_is_reported(self):
        state,_,tag,api=self.fixture()
        original_api=api
        def changed_tag_after_post(path,payload=None):
            result=original_api(path,payload)
            if payload is not None: state['current']='b'*40
            return result
        with self.assertRaisesRegex(RuntimeError,'Published tag does not match'):
            self.call(state,tag,changed_tag_after_post)
        self.assertEqual(len(state['posts']),1)


class PublicationMetadata(unittest.TestCase):
    def setUp(self):
        self.temporary=tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.before=Path.cwd()
        os.chdir(self.temporary.name)
        self.addCleanup(os.chdir,self.before)
        Path('VERSION').write_text('1.6.3\n')
        Path('CHANGELOG.md').write_text('## [1.6.3] - 2026-09-30\n')
        Path('CITATION.cff').write_text('version: 1.6.3\nurl: https://example.test/releases/tag/v1.6.3\n')
        for name in ('vaspberry.f','vaspberry_gfortran_serial.f'):
            Path(name).write_text('! PROGRAM VASPBERRY Version 1.6.3\nVASPBERRY (Ver 1.6.3)\n')
        for name in ('docs/VALIDATION_1.6.3.md','docs/releases/v1.6.3.md',
                     'tools/native_kubo_csv.py', 'tools/vaspberry_post.py','tools/postprocess_config.py',
                     'docs/POSTPROCESSING.md',
                     'examples/features/simple-postprocess/README.md',
                     'examples/features/simple-postprocess/bi.ini',
                     'examples/features/simple-postprocess/bi-rescan.ini',
                     'examples/features/procar-character/README.md',
                     'vaspberry_spin_chern.inc', 'vaspberry_spin_kubo.inc',
                     'vaspberry_help.inc', 'tests/test_native_help.py',
                     'docs/NATIVE_COMMANDS.md', 'docs/BUILD.md',
                     '.github/scripts/validate_native_spin.py',
                     'docs/SPIN_CHERN.md', 'docs/SPIN_KUBO.md',
                     'examples/materials/graphene-spin-chern/README.md',
                     'examples/materials/bi-spin-hall/spin-chern-kubo/README.md'):
            path=Path(name); path.parent.mkdir(parents=True,exist_ok=True); path.write_text('release\n')

    def test_release_and_previous_source_are_pinned(self):
        self.assertEqual(publisher.VERSION, '1.6.3')
        self.assertEqual(publisher.TAG, 'v1.6.3')
        self.assertEqual(publisher.PREVIOUS_TAG, 'v1.6.2')
        self.assertEqual(publisher.PREVIOUS_COMMIT, '411d0ebf09ad900af3e8661809d8deb20abd9c21')

    def test_frontend_guide_and_runnable_examples_must_be_packaged(self):
        for name in ('tools/native_kubo_csv.py', 'tools/vaspberry_post.py', 'tools/postprocess_config.py',
                     'docs/POSTPROCESSING.md',
                     'examples/features/simple-postprocess/README.md',
                     'examples/features/simple-postprocess/bi.ini',
                     'examples/features/simple-postprocess/bi-rescan.ini',
                     'vaspberry_spin_chern.inc', 'vaspberry_spin_kubo.inc',
                     'vaspberry_help.inc', 'tests/test_native_help.py',
                     'docs/NATIVE_COMMANDS.md', 'docs/BUILD.md',
                     '.github/scripts/validate_native_spin.py',
                     'docs/SPIN_CHERN.md', 'docs/SPIN_KUBO.md',
                     'examples/materials/graphene-spin-chern/README.md',
                     'examples/materials/bi-spin-hall/spin-chern-kubo/README.md'):
            with self.subTest(path=name):
                path = Path(name)
                original = path.read_text()
                path.unlink()
                with self.assertRaisesRegex(RuntimeError, 'Missing release documentation'):
                    publisher.validate_metadata()
                path.write_text(original)

    def test_consistent_metadata(self):
        publisher.validate_metadata()

    def test_version_runtime_changelog_citation_or_missing_example_rejected(self):
        for path,replacement in [('VERSION','1.3.0'),('CHANGELOG.md','## [Unreleased]'),
                                 ('CITATION.cff','version: 1.6.3rc1'),
                                 ('vaspberry.f','! Version 1.6.3\nVASPBERRY (Ver 1.3.0)'),
                                 ('examples/features/procar-character/README.md',None)]:
            with self.subTest(path=path):
                target=Path(path); previous=target.read_text()
                if replacement is None: target.unlink()
                else: target.write_text(replacement)
                with self.assertRaises(RuntimeError): publisher.validate_metadata()
                target.write_text(previous)

    def test_release_links_pin_tag_and_reject_escaping_repo(self):
        notes=Path('docs/releases/v1.6.3.md')
        notes.write_text('[example](../../examples/features/procar-character/README.md#commands)\n')
        result=publisher.release_body(notes,TARGET,[])
        self.assertIn('/blob/v1.6.3/examples/features/procar-character/README.md#commands',result)
        notes.write_text('[bad](../../../outside.md)\n')
        with self.assertRaisesRegex(RuntimeError,'Invalid release link'):
            publisher.release_body(notes,TARGET,[])


class PublicationLifecycle(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        before = Path.cwd()
        os.chdir(temporary.name)
        self.addCleanup(os.chdir, before)
        env = mock.patch.dict(os.environ, GITHUB_SHA=TARGET,
                              GITHUB_REPOSITORY=publisher.REPO,
                              GITHUB_REF='refs/heads/master',
                              GITHUB_STEP_SUMMARY=str(Path('step-summary.md').resolve()))
        env.start(); self.addCleanup(env.stop)
        def git(*args):
            if args == ('rev-parse', 'HEAD'): return TARGET
            if args[0] == 'merge-base': return publisher.PREVIOUS_COMMIT
            if args[0] == 'ls-tree': return 'VERSION\nvaspberry.f'
            raise AssertionError(args)
        self.git = self.patch('git', side_effect=git)
        self.api = self.patch('api', return_value={'object': {'sha': TARGET}})
        self.metadata = self.patch('validate_metadata')
        self.checked = self.patch('checked_workflow', side_effect=lambda name, target, **kw:
                                  dict(workflow=name, sha=target, run_id=123, attempt=1, jobs=[]))
        self.publish = self.patch('publish', return_value=('PUBLISHED', {'html_url': 'https://example.test/release'}))

    def patch(self, name, **kwargs):
        patch = mock.patch.object(publisher, name, **kwargs)
        value = patch.start(); self.addCleanup(patch.stop)
        return value

    def receipt(self):
        return json.loads(Path('release-validation.json').read_text())

    def test_final_checks_and_metadata_are_refreshed_before_publication(self):
        publisher.main()
        names = list(publisher.CHECKS)
        expected = [mock.call(name, TARGET) for name in names]
        expected += [mock.call(name, TARGET, timeout=0) for name in names]
        self.assertEqual(self.checked.call_args_list, expected)
        self.assertEqual(self.metadata.call_count, 2)
        self.publish.assert_called_once()
        self.assertEqual(self.receipt()['status'], 'PUBLISHED')
        self.assertEqual(self.receipt()['phase'], 'complete')

    def test_obsolete_queued_target_fails_before_waiting(self):
        self.api.return_value = {'object': {'sha': 'b' * 40}}
        with self.assertRaisesRegex(RuntimeError, 'no longer master'):
            publisher.main()
        self.checked.assert_not_called(); self.publish.assert_not_called()
        self.assertEqual(self.receipt()['status'], 'FAILED')
        self.assertFalse(self.receipt()['publication_may_have_happened'])

    def test_failed_new_attempt_is_retained_and_prevents_publication(self):
        initial = self.checked.side_effect
        def check(name, target, **kwargs):
            if 'timeout' in kwargs:
                raise RuntimeError('latest attempt failed')
            return initial(name, target)
        self.checked.side_effect = check
        with self.assertRaisesRegex(RuntimeError, 'latest attempt failed'):
            publisher.main()
        self.publish.assert_not_called()
        receipt = self.receipt()
        self.assertEqual(receipt['status'], 'FAILED')
        self.assertEqual(receipt['phase'], 'revalidation: ci.yml')
        self.assertEqual(len(receipt['checks']), 4)
        self.assertFalse(receipt['publication_may_have_happened'])

    def test_initial_validation_failure_preserves_completed_checks(self):
        initial = self.checked.side_effect
        def check(name, target, **kwargs):
            if name == 'intel-mpi-validation.yml':
                raise RuntimeError('Intel MPI failed')
            return initial(name, target)
        self.checked.side_effect = check
        with self.assertRaisesRegex(RuntimeError, 'Intel MPI failed'):
            publisher.main()
        self.publish.assert_not_called()
        self.assertEqual(len(self.receipt()['checks']), 2)
        self.assertEqual(self.receipt()['phase'], 'validation: intel-mpi-validation.yml')

    def test_final_metadata_failure_preserves_checks_without_post(self):
        self.metadata.side_effect = [None, RuntimeError('Incorrect source version')]
        with self.assertRaisesRegex(RuntimeError, 'Incorrect source version'):
            publisher.main()
        self.publish.assert_not_called()
        self.assertEqual(self.receipt()['phase'], 'final metadata')
        self.assertFalse(self.receipt()['publication_may_have_happened'])

    def test_publication_exception_preserves_uncertainty_and_original_exception(self):
        self.publish.side_effect = TimeoutError('lost response')
        with self.assertRaisesRegex(TimeoutError, 'lost response'):
            publisher.main()
        receipt = self.receipt()
        self.assertEqual(receipt['status'], 'FAILED')
        self.assertEqual(receipt['phase'], 'publication')
        self.assertEqual(receipt['error'], 'TimeoutError: lost response')
        self.assertTrue(receipt['publication_may_have_happened'])
        self.assertEqual(len(receipt['checks']), 4)


if __name__ == '__main__':
    unittest.main()
