"""Exercise public help requests on complete serial/MPI executables.

Help must work from an empty directory, avoid calculation I/O, fail clearly
for bad requests, and let every MPI rank exit without interleaved output.
"""
import os
import re
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
TASKS = ("chern", "spin-chern", "spin-kubo", "z2", "kubo", "kubo-integral",
         "kubo-pairs", "optical", "spectrum", "wavefunction", "velocity")
OPTIONS = ("bands", "mesh", "input-dir", "wavecar", "spinor", "output", "curvature-csv",
           "pairs-csv", "per-band", "spin-axis", "outcar", "energy-gap-tol",
           "spin-gap-tol", "sum-bands", "wavefunction-band", "kpoint",
           "real-grid", "imaginary", "theta", "phi", "kubo-source", "waveder", "incar")


def compile_executable(work, compiler, name, source="vaspberry.f", mpi=False):
    binary = work / name
    command = [compiler, "-cpp", "-O0", "-fcheck=all", "-ffixed-line-length-none",
               "-fallow-argument-mismatch"]
    if mpi:
        command.append("-DMPI_USE")
    command += [str(ROOT / source), "-llapack", "-lblas", "-o", str(binary)]
    result = subprocess.run(command, capture_output=True, text=True, timeout=90)
    if result.returncode:
        raise AssertionError(result.stdout + result.stderr)
    return binary


class HelpAssertions:
    def request(self, *arguments, success=True, prefix=(), binary=None):
        with tempfile.TemporaryDirectory(prefix="vaspberry-help-empty-") as folder:
            result = subprocess.run([*prefix, str(binary or self.binary), *arguments],
                                    cwd=folder, env=self.environment,
                                    capture_output=True, text=True, timeout=30)
            self.assertEqual(list(Path(folder).iterdir()), [],
                             "Help must not create calculation files")
        if success:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertEqual(result.stdout.count("PROGRAM INSTRUCTION"), 1,
                             result.stdout + result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        return result


@unittest.skipUnless(shutil.which("gfortran"), "gfortran required for native help")
class NativeHelpTests(HelpAssertions, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix="vaspberry-help-build-")
        cls.work = Path(cls.temp.name)
        cls.environment = dict(os.environ)
        cls.binary = compile_executable(cls.work, "gfortran", "vaspberry")
        cls.legacy = compile_executable(cls.work, "gfortran", "legacy",
                                       source="vaspberry_gfortran_serial.f")

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_overview_is_short_and_discoverable(self):
        result = self.request("--help")
        self.assertLessEqual(len(result.stdout.splitlines()), 50)
        for text in ("--help task", "--help kubo", "--help bands", "--help all",
                     "WAVECAR", "docs/NATIVE_COMMANDS.md"):
            self.assertIn(text, result.stdout)
        self.assertIn("Ver " + (ROOT / "VERSION").read_text().strip(), result.stdout)
        self.assertEqual(result.stdout, self.request("-h").stdout)

    def test_task_catalog_lists_supported_calculations(self):
        result = self.request("--help", "task")
        for task in TASKS:
            self.assertIn(task, result.stdout)
        self.assertEqual(result.stdout, self.request("--help", "tasks").stdout)
        self.assertEqual(result.stdout, self.request("--help", "--task").stdout)

    def test_each_task_has_a_specific_command_example(self):
        for task in TASKS:
            with self.subTest(task=task):
                result = self.request("--help", task)
                self.assertIn("--task " + task, result.stdout)
                self.assertIn("WAVECAR", result.stdout)

    def test_wrapped_examples_are_shell_continuations(self):
        for topic in (*TASKS, "examples"):
            with self.subTest(topic=topic):
                lines = self.request("--help", topic).stdout.splitlines()
                commands = 0
                for index, line in enumerate(lines):
                    if not re.match(r'\s*(?:vaspberry|build/vaspberry[\w-]*|mpiexec|"\$VB")\s', line):
                        continue
                    commands += 1
                    while index + 1 < len(lines) and re.match(r"\s+--", lines[index + 1]):
                        self.assertTrue(lines[index].rstrip().endswith("\\"),
                                        "Wrapped command cannot be pasted into a shell:\n" +
                                        lines[index] + "\n" + lines[index + 1])
                        index += 1
                self.assertGreater(commands, 0, "No executable example in task help")

    def test_option_keywords_and_option_spelling_are_equivalent(self):
        for option in OPTIONS:
            with self.subTest(option=option):
                result = self.request("--help", option)
                self.assertEqual(result.stdout, self.request("--help", "--" + option).stdout)
                self.assertIn("--" + option, result.stdout)
        self.assertEqual(self.request("--help", "kubo").stdout,
                         self.request("-h", "kubo").stdout)

    def test_kubo_help_describes_input_output_and_degeneracy(self):
        text = self.request("--help", "kubo").stdout
        for item in ("WAVECAR", "--bands", "--curvature-csv", "--per-band",
                     "KUBO_WAVEDER.csv", "Kubo", "PAW"):
            self.assertIn(item, text)
        self.assertIn("degenera", text.lower())
        self.assertIn("path", text.lower())

    def test_kubo_source_help_requires_explicit_approximation(self):
        source = self.request("--help", "kubo-source").stdout
        for term in ("waveder", "wavecar", "--kubo-source"):
            self.assertIn(term, source.lower())
        self.assertIn("approx", source.lower())
        self.assertIn("default", source.lower())

    def test_input_directory_help_describes_independent_cwd_relative_overrides(self):
        for topic in ('input-dir', 'kubo', 'wavecar', 'waveder', 'incar', 'outcar'):
            text = self.request('--help', topic).stdout
            self.assertIn('--input-dir', text)
            self.assertIn('cwd', text.lower())
            self.assertNotIn('siblings of --wavecar', text)
            self.assertNotIn('beside --wavecar', text)
        text = self.request('--help', 'input-dir', binary=self.legacy).stdout
        self.assertIn('-f PATH', text)
        self.assertIn('outputs stay in cwd', text)

    def test_mesh_and_band_help_distinguish_explicit_control_requests(self):
        mesh = ' '.join(self.request('--help', 'mesh').stdout.split())
        self.assertIn('only explicit --mesh requests integration', mesh)
        self.assertIn('legacy -kx/-ky alone', mesh)
        bands = ' '.join(self.request('--help', 'bands').stdout.split())
        self.assertIn('Prefer one --bands selection', bands)
        self.assertIn('historical endpoint/singleton precedence', bands)
        kubo = ' '.join(self.request('--help', 'kubo').stdout.split())
        self.assertIn('rejects wavefunction/optical/folding controls', kubo)
        self.assertIn('-sigma', kubo)

    def test_pairs_help_points_to_executable_postprocessor(self):
        for topic in ("kubo-pairs", "outputs"):
            text = self.request("--help", topic).stdout
            self.assertIn("tools/vaspberry_kubo.py", text)
            self.assertNotIn("tools/kubo_hall_workflow.py", text)

    def test_spin_tasks_name_their_actual_files(self):
        chern = self.request("--help", "spin-chern").stdout
        for name in ("WAVECAR", "OUTCAR", "SPIN_CHERN.csv", "SPIN_BERRY.csv",
                     "SPIN_SPECTRUM.csv", "--spin-axis"):
            self.assertIn(name, chern)
        kubo = self.request("--help", "spin-kubo").stdout
        for name in ("WAVECAR", "OUTCAR", "SPIN_KUBO.csv", "SPIN_KUBO_SPECTRUM.csv",
                     "SPIN_KUBO_INTEGRAL.csv", "--mesh", "--sum-bands"):
            self.assertIn(name, kubo)

    def test_build_outputs_and_examples_are_discoverable(self):
        self.assertIn("make ifx-mpi", self.request("--help", "build").stdout)
        self.assertIn("--curvature-csv", self.request("--help", "outputs").stdout)
        self.assertIn("--task", self.request("--help", "examples").stdout)
        self.assertIn("--bands", self.request("--help", "options").stdout)

    def test_full_legacy_reference_remains_accessible(self):
        for keyword in ("all", "legacy"):
            with self.subTest(keyword=keyword):
                result = self.request("--help", keyword)
                for text in ("-ii", "--per-band", "-kubo_pairs", "-assume byterecl"):
                    self.assertIn(text, result.stdout)
                self.assertGreater(len(result.stdout.splitlines()), 100)

    def test_removed_option_help_and_preflight(self):
        text = self.request("--help", "bundle").stdout
        self.assertIn("removed", text)
        self.assertIn("--per-band 1", text)
        self.assertIn("For either source", text)
        self.assertIn("--bands 31,33:34", text)
        self.assertNotIn("Individual output requires --kubo-source wavecar", text)
        for binary in (self.binary, self.legacy):
            for flag in ("--bundle", "-kubo_bundle"):
                for suffix in ((), ("0",), ("1",)):
                    result = self.request(flag, *suffix, binary=binary, success=False)
                    self.assertIn("removed", result.stderr)
                    self.assertNotIn("Error opening", result.stdout + result.stderr)

    def test_topics_have_readable_sections_and_bounded_width(self):
        for topic in (*TASKS, *OPTIONS, 'task', 'options', 'outputs', 'examples',
                      'ien', 'fen', 'nediv', 'sigma', 'build', 'bundle', 'all'):
            with self.subTest(topic=topic):
                text = self.request('--help', topic).stdout
                self.assertNotIn('\x1b', text)
                self.assertTrue(all(len(line.expandtabs(8)) <= 80
                                    for line in text.splitlines()), text)
                self.assertGreater(text.count('\n\n'), 1, text)
                if topic in TASKS:
                    for heading in ('Run', 'Inputs', 'Outputs', 'Conditions', 'More'):
                        self.assertIn('\n'+heading+'\n', text)

    def test_corrected_help_contracts_and_required_auxiliary_inputs(self):
        full = ' '.join(self.request('--help', 'all').stdout.split())
        for value in ('-sw_file', 'BERRYCURV.filename.dat', '(2*np-1)',
                      'LEFT/RIGHT intensity', '--kubo-source wavecar'):
            self.assertIn(value, full)
        self.assertNotIn(' -sw ', full)
        self.assertNotIn("unless '-cd' is not 1", full)
        waveform = ' '.join(self.request('--help', 'wavefunction').stdout.split())
        for value in ('POSCAR/EIGENVAL are required', '--input-dir does not relocate',
                      'amplitudes', 'Gamma'):
            self.assertIn(value, waveform)
        bands = ' '.join(self.request('--help', 'bands').stdout.split())
        self.assertIn('N:N is a singleton', bands)
        self.assertIn('With explicit --kubo-source wavecar', bands)
        for topic in ('spin-chern', 'spin-kubo'):
            self.assertIn('ISPIN=1', self.request('--help', topic).stdout)

    def test_bad_keyword_does_not_silently_fall_back(self):
        for keyword in ("typo-topic", "", "k" * 400):
            with self.subTest(keyword=keyword):
                result = self.request("--help", keyword, success=False)
                message = result.stdout + result.stderr
                self.assertIn("help", message.lower())
                self.assertNotIn("Error opening", message)

    def test_extra_help_arguments_fail_before_calculation(self):
        for arguments in (("--help", "kubo", "bands"),
                          ("--help", "kubo", "--bands", "18"),
                          ("-h", "kubo", "extra")):
            with self.subTest(arguments=arguments):
                result = self.request(*arguments, success=False)
                self.assertIn("help", (result.stdout + result.stderr).lower())

    def test_help_mixed_with_calculation_flags_is_rejected(self):
        for arguments in (("--task", "kubo", "--help"),
                          ("--task", "kubo", "--help", "bands"),
                          ("-kubo", "2", "-h", "kubo")):
            with self.subTest(arguments=arguments):
                result = self.request(*arguments, success=False)
                self.assertIn("help", (result.stdout + result.stderr).lower())

    def test_help_text_in_a_value_position_is_not_a_help_request(self):
        result = self.request("--wavecar", "--help", success=False)
        self.assertNotIn("PROGRAM INSTRUCTION", result.stdout)
        self.assertIn("--wavecar", result.stdout + result.stderr)

    def test_historical_serial_source_does_not_advertise_kubo_support(self):
        overview = self.request("--help", binary=self.legacy)
        self.assertLessEqual(len(overview.stdout.splitlines()), 50)
        self.request("--help", "all", binary=self.legacy)
        for topic in ("kubo", "spin-chern", "spin-kubo"):
            with self.subTest(topic=topic):
                result = self.request("--help", topic, success=False, binary=self.legacy)
                self.assertIn("legacy", (result.stdout + result.stderr).lower())


@unittest.skipUnless(shutil.which("mpifort") and shutil.which("mpiexec"),
                     "MPI compiler/runtime required for multi-rank help")
class MpiHelpTests(HelpAssertions, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        wrapper = subprocess.run(["mpifort", "--version"], capture_output=True,
                                 text=True, timeout=15)
        if "GNU Fortran" not in wrapper.stdout + wrapper.stderr:
            raise unittest.SkipTest("GNU MPI wrapper required; Intel help uses make checks")
        cls.temp = tempfile.TemporaryDirectory(prefix="vaspberry-mpi-help-build-")
        cls.work = Path(cls.temp.name)
        cls.environment = dict(os.environ, OMPI_ALLOW_RUN_AS_ROOT="1",
                               OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1",
                               OMPI_MCA_rmaps_base_oversubscribe="1")
        cls.binary = compile_executable(cls.work, "mpifort", "vaspberry-mpi", mpi=True)
        cls.prefix = ("mpiexec", "-n", "2")

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_overview_and_task_help_print_once_and_all_ranks_exit(self):
        for arguments in (("--help",), ("--help", "kubo"), ("--help", "task"),
                          ("--help", "all")):
            with self.subTest(arguments=arguments):
                self.request(*arguments, prefix=self.prefix)

    def test_removed_selectors_fail_once_before_input_on_all_ranks(self):
        for arguments in (("--bundle",), ("--bundle", "0"),
                          ("--task", "kubo", "-kubo_bundle", "1")):
            result = self.request(*arguments, success=False, prefix=self.prefix)
            self.assertEqual(result.stderr.count("was removed"), 1, result.stderr)
            self.assertNotIn("Error opening", result.stdout + result.stderr)

    def test_invalid_help_does_not_hang_or_start_a_calculation(self):
        for arguments in (("--help", "typo-topic"),
                          ("--help", "kubo", "extra"),
                          ("--task", "kubo", "--help")):
            with self.subTest(arguments=arguments):
                result = self.request(*arguments, success=False, prefix=self.prefix)
                self.assertIn("help", (result.stdout + result.stderr).lower())
                self.assertNotIn("MPI_ABORT", result.stderr)


if __name__ == "__main__":
    unittest.main()
