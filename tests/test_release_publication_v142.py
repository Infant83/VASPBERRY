"""Publication must consume exact-source CI and never move or rewrite releases."""
import copy
import importlib.util
import os
from pathlib import Path
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('publish_v142', ROOT/'.github/scripts/publish_v142.py')
publisher = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(publisher)
TARGET = 'a'*40


class ReleaseChecks(unittest.TestCase):
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
            if path == 'releases?per_page=100': return [release] if existing else []
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
        Path('VERSION').write_text('1.4.2\n')
        Path('CHANGELOG.md').write_text('## [1.4.2] - 2026-09-26\n')
        Path('CITATION.cff').write_text('version: 1.4.2\nurl: https://example.test/releases/tag/v1.4.2\n')
        for name in ('vaspberry.f','vaspberry_gfortran_serial.f'):
            Path(name).write_text('! PROGRAM VASPBERRY Version 1.4.2\nVASPBERRY (Ver 1.4.2)\n')
        for name in ('docs/VALIDATION_1.4.2.md','docs/releases/v1.4.2.md','examples/features/procar-character/README.md'):
            path=Path(name); path.parent.mkdir(parents=True,exist_ok=True); path.write_text('release\n')

    def test_consistent_metadata(self):
        publisher.validate_metadata()

    def test_version_runtime_changelog_citation_or_missing_example_rejected(self):
        for path,replacement in [('VERSION','1.3.0'),('CHANGELOG.md','## [Unreleased]'),
                                 ('CITATION.cff','version: 1.4.2rc1'),
                                 ('vaspberry.f','! Version 1.4.2\nVASPBERRY (Ver 1.3.0)'),
                                 ('examples/features/procar-character/README.md',None)]:
            with self.subTest(path=path):
                target=Path(path); previous=target.read_text()
                if replacement is None: target.unlink()
                else: target.write_text(replacement)
                with self.assertRaises(RuntimeError): publisher.validate_metadata()
                target.write_text(previous)

    def test_release_links_pin_tag_and_reject_escaping_repo(self):
        notes=Path('docs/releases/v1.4.2.md')
        notes.write_text('[example](../../examples/features/procar-character/README.md#commands)\n')
        result=publisher.release_body(notes,TARGET,[])
        self.assertIn('/blob/v1.4.2/examples/features/procar-character/README.md#commands',result)
        notes.write_text('[bad](../../../outside.md)\n')
        with self.assertRaisesRegex(RuntimeError,'Invalid release link'):
            publisher.release_body(notes,TARGET,[])


if __name__ == '__main__':
    unittest.main()
