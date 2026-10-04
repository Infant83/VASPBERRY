"""Native aliases and full-executable output contracts.

Parser tests reject malformed commands before file I/O. Independent tiny
WAVECAR fixtures check occupation-free pairs and SI velocity expectations;
actual material equivalence runs are recorded separately in runs/.
"""
import re
import os
import time
import shutil
import struct
import cmath
import math
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def run_legacy_contract(command, *args, **kwargs):
    """These regressions test the retained WAVECAR operator, explicitly selected.

    Standard WAVEDER dispatch and its no-fallback guards have independent
    full-executable coverage in test_fortran_waveder_default.py.
    """
    kubo_task = any(command[i] == '--task' and command[i+1].startswith('kubo')
                    for i in range(len(command)-1))
    if (kubo_task or '-kubo' in command) and '--kubo-source' not in command:
        command = [*command, '--kubo-source', 'wavecar']
    return subprocess.run(command, *args, **kwargs)



@unittest.skipUnless(shutil.which("gfortran"), "gfortran required for native parser")
class NativeCliTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix="vaspberry-native-cli-")
        cls.work = Path(cls.temp.name)
        source = (ROOT / "vaspberry.f").read_text()
        match = re.search(r"(?ims)^ {6}subroutine parse\b.*?^ {6}end subroutine parse\s*$", source)
        if not match:
            raise AssertionError("production parser not found")
        (cls.work / "parse.f").write_text(match.group(0) + "\n")
        driver = (ROOT / "tests/fortran/test_kubo_parser.f90").read_text()
        driver = driver.replace("  logical flag_atom_project", """
  character(256) spinchern_outcar
  common /spinchern_chars/ spinchern_outcar
  logical flag_atom_project
""")
        driver = driver.replace("  write(*,'(A)')trim(foname)", """
  write(*,'(A)')trim(filename),trim(foname),trim(kubo_csv),trim(kubo_pairs),trim(spinchern_outcar)
  write(*,*)nkx,nky,ispinor,icd,ixt,ivel,ikubo,iz,ihf,nini,nmax,nn,kperiod,it,iskp,ine
  write(*,*)iwf,ikwf,ng,imag,nediv,ikubo_bundle
  write(*,*)theta,phi,init_e,fina_e,sigma,rs
""")
        driver = driver.replace("  stop 1", "  print *, 'HELP'\n  stop 0")
        (cls.work / "driver.f90").write_text(driver)
        cls.binary = cls.work / "parser"
        built = run_legacy_contract(["gfortran", "-cpp", "-O0", "-fcheck=all",
                                "-ffixed-line-length-none", "parse.f", "driver.f90",
                                "-o", str(cls.binary)], cwd=cls.work, capture_output=True, text=True)
        if built.returncode:
            raise AssertionError(built.stdout + built.stderr)

        legacy_source = (ROOT / "vaspberry_gfortran_serial.f").read_text()
        legacy_parse = re.search(r"(?ims)^ {6}subroutine parse\b.*?^ {6}end subroutine parse\s*$", legacy_source)
        if not legacy_parse:
            raise AssertionError("historical production parser not found")
        (cls.work / "legacy-parse.f").write_text(legacy_parse.group(0) + "\n")
        legacy_driver = """program test_legacy_input
 implicit none
 character(75) filename,foname,fbz,ver_tag
 integer nx,ny,spin,cd,ext,vel,z2,hf,first,last,kperiod,it,skip,ne
 integer wf,k,ng(3),imag
 real(8) rs(3)
 nx=2;ny=2;spin=2;kperiod=1;ver_tag='legacy parser'
 call parse(filename,foname,nx,ny,spin,cd,ext,fbz, &
  vel,z2,hf,first,last,kperiod,it,skip,ne,ver_tag,wf,k,ng,rs,imag)
 write(*,'(A)')trim(filename),trim(foname)
end program
subroutine help(tag)
 character(*) tag
 stop 1
end subroutine
subroutine vaspberry_fail
 stop 2
end subroutine
"""
        (cls.work / "legacy-driver.f90").write_text(legacy_driver)
        cls.legacy_parser = cls.work / "legacy-parser"
        built = subprocess.run(["gfortran", "-O0", "-fcheck=all", "-ffixed-line-length-none",
                                "legacy-parse.f", "legacy-driver.f90", "-o", str(cls.legacy_parser)],
                               cwd=cls.work, capture_output=True, text=True)
        if built.returncode:
            raise AssertionError(built.stdout + built.stderr)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def run_parser(self, arguments, success=True):
        result = run_legacy_contract([str(self.binary), *arguments], cwd=self.work,
                                capture_output=True, text=True)
        if success:
            self.assertEqual(result.returncode, 0, result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def equivalent(self, new, legacy):
        self.assertEqual(self.run_parser(new).stdout, self.run_parser(legacy).stdout)

    def test_each_task_dispatch_and_argument_order(self):
        for task, legacy, needed in [
            ("chern", [], []), ("z2", ["-z2", "1"], []),
            ("kubo", ["-kubo", "2"], []), ("kubo-line", ["-kubo", "2"], []),
            ("kubo-integral", ["-kubo", "1"], []),
            ("kubo-pairs", ["-kubo", "2"], ["-kubo_pairs", "pairs.csv"]),
            ("optical", ["-cd", "1"], []), ("spectrum", ["-cd", "2"], []),
            ("velocity", ["-vel", "1"], []),
            ("wavefunction", [], ["-wf", "18"]),
        ]:
            with self.subTest(task=task):
                self.equivalent(["--task", task, *needed], [*legacy, *needed])
                self.equivalent([*needed, "--task", task], [*legacy, *needed])
                self.equivalent(["--task", task, *legacy, *needed], [*legacy, *needed])

    def test_mesh_bands_paths_spinor_and_output_aliases(self):
        self.equivalent(["--task", "chern", "--wavecar", "/tmp/a b/WAVECAR",
                         "--mesh", "12,16", "--bands", "1:18", "--spinor", "1",
                         "--output", "sample"],
                        ["-f", "/tmp/a b/WAVECAR", "-kx", "12", "-ky", "16",
                         "-ii", "1", "-if", "18", "-s", "1", "-o", "sample"])
        self.equivalent(["--bands", "18"], ["-is", "18"])

    def test_input_directory_applies_to_non_kubo_defaults_without_changing_other_paths(self):
        directory = self.work/'input directory';directory.mkdir()
        for task in ('chern', 'velocity', 'optical', 'spin-chern'):
            flags = ['--task', task, '--bands', '1:2', '--input-dir', str(directory)]
            result = self.run_parser(flags)
            self.assertEqual(result.stdout.splitlines()[0], str(directory/'WAVECAR'))
            if task == 'spin-chern':
                self.assertIn(str(directory/'OUTCAR'), result.stdout)
                for order in (0, 1):
                    override = ['--outcar', 'explicit.outcar']
                    actual = self.run_parser([*override, *flags] if order else [*flags, *override])
                    self.assertIn('explicit.outcar', actual.stdout)
                    self.assertNotIn(str(directory/'OUTCAR'), actual.stdout)
            for override in (['--wavecar', 'relative.wave'], ['-f', 'relative.wave']):
                for order in (0, 1):
                    actual = self.run_parser([*override, *flags] if order else [*flags, *override])
                    self.assertEqual(actual.stdout.splitlines()[0], 'relative.wave')
        self.equivalent(['--task', 'chern', '--input-dir', str(directory),
                         '--wavecar', 'custom', '--output', 'result'],
                        ['--task', 'chern', '--wavecar', 'custom', '--output', 'result'])

    def test_historical_parser_input_directory_and_explicit_legacy_file(self):
        directory = self.work/'old-input';directory.mkdir()
        commands = [(['--input-dir', 'old-input'], 'old-input/WAVECAR'),
                    (['--input-dir', 'old-input/', '-o', 'result'], 'old-input/WAVECAR'),
                    (['--input-dir', 'old-input', '-f', 'override.wave'], 'override.wave'),
                    (['-f', 'override.wave', '--input-dir', 'old-input'], 'override.wave')]
        for args, filename in commands:
            result = subprocess.run([str(self.legacy_parser), *args], cwd=self.work,
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(result.stdout.splitlines()[0], filename)
            self.assertEqual(result.stdout.splitlines()[1], 'BERRYCURV.result' if '-o' in args else 'BERRYCURV')
        (self.work/'regular-input-file').write_text('not a directory')
        for args in (['--input-dir', 'missing', '-f', 'override.wave'],
                     ['--input-dir', 'regular-input-file'], ['--input-dir', ''],
                     ['-f', ''], ['--input-dir', 'a'*76], ['-f', 'a'*76]):
            result = subprocess.run([str(self.legacy_parser), *args], cwd=self.work,
                                    capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('input', result.stderr)

    def test_directory_probe_accepts_directory_links_and_rejects_file_links(self):
        base = self.work / 'directory probe'
        base.mkdir()
        (base / 'valid input').mkdir()
        (base / 'directory alias').symlink_to('valid input', target_is_directory=True)
        (base / 'regular file').write_text('not a directory')
        (base / 'file alias').symlink_to('regular file')
        (base / 'broken alias').symlink_to('missing')
        valid = ('.', 'directory probe/valid input', 'directory probe/directory alias')
        invalid = ('directory probe/regular file', 'directory probe/file alias',
                   'directory probe/missing', 'directory probe/broken alias')
        for binary in (self.binary, self.legacy_parser):
            for path in (*valid, *invalid):
                with self.subTest(binary=binary.name, path=path):
                    result = subprocess.run([str(binary), '--input-dir', path,
                                             '-f', 'explicit.wave'], cwd=self.work,
                                            capture_output=True, text=True)
                    if path in valid:
                        self.assertEqual(result.returncode, 0, result.stderr)
                        self.assertEqual(result.stdout.splitlines()[0], 'explicit.wave')
                    else:
                        self.assertNotEqual(result.returncode, 0)
                        self.assertIn('must be an existing directory', result.stderr)

    def test_selected_range_defaults_to_trace_and_per_band_is_explicit(self):
        trace = ["--task", "kubo", "--bands", "1:18"]
        for equivalent in (["--bands", "1:18", "-kubo", "2"],
                           ["--task", "kubo", "-ii", "1", "-if", "18"],
                           [*trace, "--per-band", "0"],
                           [*trace, "--curvature-csv", "KUBO.csv"]):
            self.equivalent(trace, equivalent)
        self.assertIn("KUBO.csv", self.run_parser(trace).stdout)
        for task, flag in (("kubo", "2"), ("kubo-line", "2"), ("kubo-integral", "1")):
            self.equivalent(["--task", task, "--bands", "1:18", "--per-band", "1"],
                            ["-kubo", flag, "-ii", "1", "-if", "18"])
            self.equivalent(["--task", task, "--bands", "18:18"],
                            ["-kubo", flag, "-is", "18"])
        self.equivalent(["--pairs-csv", "pairs.csv", "--task", "kubo-pairs"],
                        ["-kubo", "2", "-kubo_pairs", "pairs.csv"])

    def test_removed_selectors_give_migration_errors(self):
        for flag in ("--bundle", "-kubo_bundle"):
            for suffix in ([], ["0"], ["1"], ["garbage"]):
                for prefix in ([], ["--task", "kubo"]):
                    result = self.run_parser([*prefix, flag, *suffix], False)
                    self.assertIn("removed", result.stderr)
                    self.assertIn("--bands FIRST:LAST", result.stderr)
                    self.assertIn("--per-band 1", result.stderr)
        # Removed-option text at a legacy filename value position is data.
        self.run_parser(["-f", "--bundle", "-kubo", "2"])

    def test_per_band_requires_modern_ordinary_kubo_task(self):
        for value in ("0", "1"):
            for task in ("chern", "spin-chern", "spin-kubo", "z2", "optical", "velocity"):
                result = self.run_parser(["--task", task, "--bands", "1:2",
                                          "--per-band", value], False)
                self.assertIn("--per-band requires", result.stderr)
            for args in (["-kubo", "2", "--bands", "1:2"],
                         ["--task", "kubo-pairs", "--pairs-csv", "p.csv"],
                         ["--task", "kubo", "--pairs-csv", "p.csv"]):
                result = self.run_parser([*args, "--per-band", value], False)
                self.assertIn("--per-band requires", result.stderr)

    def test_modern_band_overrides_warn_and_preserve_effective_selection(self):
        cases = [(['--bands', '1:2', '--bands', '1:2'], (1, 2)),
                 (['--bands', '1:2', '--bands', '3'], (3, 3)),
                 (['--bands', '1:2', '-ii', '2'], (2, 2)),
                 (['--bands', '1:2', '-if', '3'], (1, 3)),
                 (['--bands', '1:2', '-is', '3'], (3, 3)),
                 (['-ii', '3', '-if', '4', '--bands', '1:2'], (1, 2)),
                 (['-is', '3', '--bands', '1:2'], (1, 2))]
        for source in ('waveder', 'wavecar'):
            for selectors, expected in cases:
                with self.subTest(source=source, selectors=selectors):
                    result = self.run_parser(['--task', 'kubo', '--kubo-source', source,
                                              *selectors])
                    self.assertIn('WARNING: repeated/mixed band selectors', result.stderr)
                    self.assertIn('Prefer one --bands', result.stderr)
                    values = result.stdout.splitlines()[-3].split()
                    self.assertEqual(tuple(map(int, values[9:11])), expected)
            # The historical endpoint pair is one selection and has no override warning.
            for selectors in (['-ii', '1', '-if', '2'], ['-if', '2', '-ii', '1']):
                result = self.run_parser(['--task', 'kubo', '--kubo-source', source, *selectors])
                self.assertNotIn('repeated/mixed band selectors', result.stderr)
        # Preserve the published comma-list directionality: a later complete modern
        # selection can replace a legacy one; a legacy endpoint cannot edit a list.
        result = self.run_parser(['--task', 'kubo', '--kubo-source', 'waveder',
                                  '-is', '2', '--bands', '1,3'])
        self.assertIn('WARNING: repeated/mixed band selectors', result.stderr)
        for flag in ('-ii', '-if', '-is'):
            result = self.run_parser(['--task', 'kubo', '--kubo-source', 'waveder',
                                      '--bands', '1,3', flag, '2'], False)
            self.assertIn('comma-list', result.stderr)

    def test_modern_waveder_kubo_rejects_ignored_task_controls(self):
        controls = [('--kpoint', '2'), ('-k', '2'), ('--real-grid', '8,8,8'),
                    ('-ng', '8,8,8'), ('--imaginary', '0'), ('-im', '0'),
                    ('-ishift', '0,0,0'), ('--theta', '0'), ('-theta', '0'),
                    ('--phi', '0'), ('-phi', '0'), ('-ien', '0'), ('-fen', '10'),
                    ('-nediv', '1000'), ('-sigma', '0.01'), ('-kp', '1'), ('-skp', '0')]
        for task in ('kubo', 'kubo-line', 'kubo-integral'):
            base = ['--task', task, '--kubo-source', 'waveder', '--mesh', '2,2', '--bands', '1:2']
            for flag, value in controls:
                for arguments in ([flag, value, *base], [*base, flag, value]):
                    with self.subTest(task=task, arguments=arguments):
                        result = self.run_parser(arguments, False)
                        self.assertIn(flag, result.stderr)
                        self.assertIn('not used by standard WAVEDER Kubo', result.stderr)
                        self.assertNotIn('Error opening', result.stdout + result.stderr)
        # Retain the explicit canonical and unrelated historical parser behavior.
        self.run_parser(['--task', 'kubo', '--kubo-source', 'wavecar', '--theta', '0'])
        self.run_parser(['--task', 'wavefunction', '--wavefunction-band', '1', '--kpoint', '2'])
        self.run_parser(['--task', 'optical', '--theta', '0'])

    def test_wavefunction_and_angle_aliases(self):
        self.equivalent(["--task", "wavefunction", "--wavefunction-band", "18",
                         "--kpoint", "2", "--real-grid", "8,10,12", "--imaginary", "1"],
                        ["-wf", "18", "-k", "2", "-ng", "8,10,12", "-im", "1"])
        self.equivalent(["--task", "optical", "--theta", "-45.5", "--phi", "1e2"],
                        ["-cd", "1", "-theta", "-45.5", "-phi", "1e2"])

    def test_task_conflicts_rejected_in_both_orders(self):
        for task, flags in [("chern", ["-kubo", "2"]), ("kubo", ["-kubo", "1"]),
                            ("z2", ["-cd", "1"]), ("optical", ["-vel", "1"]),
                            ("wavefunction", ["-wf", "18", "-z2", "1"]),
                            ("kubo", ["-wf", "18"]), ("chern", ["-ixt", "1"]),
                            ("chern", ["-t", "1"])]:
            for args in [["--task", task, *flags], [*flags, "--task", task]]:
                with self.subTest(args=args):
                    result = self.run_parser(args, False)
                    self.assertIn("conflicts with legacy task", result.stderr)

    def test_missing_requirements_and_existing_csv_guards(self):
        for args in [["--task", "wavefunction"], ["--task", "kubo-pairs"],
                     ["--task", "chern", "--per-band", "1"],
                     ["--task", "chern", "--curvature-csv", "c.csv"],
                     ["--task", "kubo-pairs", "--pairs-csv", "p.csv", "--bands", "1:2"],
                     ["--task", "kubo", "--wavecar", "same", "--curvature-csv", "same"],
                     ["--task", "kubo-pairs", "--pairs-csv", "p.csv", "--curvature-csv", "c.csv"]]:
            with self.subTest(args=args):
                self.run_parser(args, False)

    def test_bad_task_and_malformed_values_rejected(self):
        for args in [["--task", "unknown"], ["--task", "chern", "--task", "chern"],
                     ["--task"], ["--mesh", "--task", "chern"], ["--unknown", "1"],
                     ["--task", "chern", "-typo", "1"],
                     ["-typo", "1", "--task", "chern"]]:
            self.run_parser(args, False)
        for option, bad_values in {
            "--mesh": ["0,2", "-1,2", "2", "2,3,4", "2,,3", "2,", ",2", "2/3", "2,3x", "2, 3"],
            "--bands": ["0", "3:2", "1:", ":2", "1::2", "1:2:3", "1,2", "1/2", "1x"],
            "--real-grid": ["1,2", "1,2,3,4", "1,,2", "1,0,2"],
            "--spinor": ["0", "3", "1,2", "1/", "abc"],
            "--wavefunction-band": ["0", "-1", "1.2"],
            "--kpoint": ["0", "-1"], "--per-band": ["2", "0,1", "-1", "1.0"],
            "--imaginary": ["2"], "--theta": ["abc", "1,2", "1/", "NaN", "Inf", "1e999"],
        }.items():
            for value in bad_values:
                with self.subTest(option=option, value=value):
                    self.run_parser([option, value], False)

    def test_serial_build_replaces_newer_legacy_binary_with_alias(self):
        project = self.work / "build-contract"
        project.mkdir()
        shutil.copy2(ROOT / "Makefile", project / "Makefile")
        (project / "vaspberry.f").write_text("      program fixture\n      end\n")
        (project / "vaspberry_spin_chern.inc").write_text("! build prerequisite fixture\n")
        (project / "vaspberry_spin_kubo.inc").write_text("! build prerequisite fixture\n")
        (project / "vaspberry_help.inc").write_text("! build prerequisite fixture\n")
        (project / "vaspberry_waveder.inc").write_text("! build prerequisite fixture\n")
        command = ["make", "serial", "FC=gfortran", "GNU_FLAGS=-O0", "GNU_LIBS="]
        for attempt in range(2):
            result = run_legacy_contract(command, cwd=project, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            canonical = project / "build/vaspberry"
            compatibility = project / "build/vaspberry-gfortran"
            self.assertTrue(canonical.is_file())
            self.assertTrue(compatibility.is_symlink())
            self.assertEqual(compatibility.resolve(), canonical.resolve())
            if attempt == 0:
                compatibility.unlink()
                compatibility.write_text("older build under the former canonical name")
                future = time.time() + 600
                os.utime(compatibility, (future, future))

    def test_standalone_help_is_success(self):
        self.assertIn("HELP", self.run_parser(["--help"]).stdout)
        self.assertIn("HELP", self.run_parser(["-h"]).stdout)


@unittest.skipUnless(shutil.which("gfortran"), "gfortran required for native execution")
class NativeExecutionTests(unittest.TestCase):
    """Run the complete executable on independent, tiny WAVECAR records.

    Pair vertices do not use occupations, unlike automatic occupied-subspace
    selection. These are format/dispatch fixtures, not material calculations.
    """

    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix="vaspberry-native-pair-occupations-")
        cls.work = Path(cls.temp.name)
        cls.binary = cls.work / "vaspberry"
        built = run_legacy_contract([
            "gfortran", "-cpp", "-O0", "-fcheck=all", "-ffixed-line-length-none",
            "-fallow-argument-mismatch", str(ROOT / "vaspberry.f"),
            "-llapack", "-lblas", "-o", str(cls.binary),
        ], capture_output=True, text=True)
        if built.returncode:
            raise AssertionError(built.stdout + built.stderr)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def write_fixture(self, path, variable_occupations, spinor_components=2):
        stride, nk, nb = 512, 2, 2
        data = bytearray(stride * (2 + nk * (nb + 1)))
        struct.pack_into("<3d", data, 0, stride, 1, 45200)
        struct.pack_into("<12d", data, stride, nk, nb, 120.,
                         1., 0., 0., 0., 1., 0., 0., 0., 1.)
        for ik, point in enumerate([(0., 0., 0.), (.25, .25, 0.)]):
            # Independently enumerate the scalar plane waves; double for SOC.
            count = sum(sum((2 * math.pi * (point[j] + g[j])) ** 2 for j in range(3))
                        / .262465831 < 120.
                        for g in [(x, y, z) for x in range(-2, 3)
                                  for y in range(-2, 3) for z in range(-2, 3)])
            npw = spinor_components * count
            record = 2 + ik * (nb + 1)
            occ = .1 if ik == 1 and variable_occupations else 1.
            struct.pack_into("<10d", data, record * stride,
                             npw, *point, -1., 0., occ, 1., 0., 0.)
            for band in range(nb):
                # Two orthonormal Fourier states over the retained coefficients.
                for ipw in range(npw):
                    coeff = cmath.exp(2j * math.pi * band * ipw / npw) / math.sqrt(npw)
                    struct.pack_into("<2f", data, (record + 1 + band) * stride + ipw * 8,
                                     coeff.real, coeff.imag)
        path.write_bytes(data)

    def write_gap_fixture(self, path, energies):
        """ISPIN channels/k points/bands, independent orthonormal plane waves."""
        stride, nk, nb = 1024, 2, 3
        data = bytearray(stride * (2 + len(energies) * nk * (nb + 1)))
        struct.pack_into("<3d", data, 0, stride, len(energies), 45200)
        struct.pack_into("<12d", data, stride, nk, nb, 120.,
                         1., 0., 0., 0., 1., 0., 0., 0., 1.)
        for isp, channel in enumerate(energies):
            for ik, point in enumerate(((.25, .25, 0.), (-.25, -.25, 0.))):
                record = 2 + (isp * nk + ik) * (nb + 1)
                values = [3., *point]
                for ib, energy in enumerate(channel[ik]):
                    values += [energy, 0., 1. if ib < 2 else 0.]
                struct.pack_into("<13d", data, record * stride, *values)
                for band in range(nb):
                    for ipw in range(3):
                        coeff = cmath.exp(2j * math.pi * band * ipw / 3) / math.sqrt(3)
                        struct.pack_into("<2f", data, (record + 1 + band) * stride + ipw * 8,
                                         coeff.real, coeff.imag)
        path.write_bytes(data)

    def execute_gap_case(self, name, energies, flags, success=True, prefix=(), binary=None):
        case = self.work / name
        case.mkdir()
        self.write_gap_fixture(case / "WAVECAR", energies)
        result = run_legacy_contract([*prefix, str(binary or self.binary), "--spinor", "1", *flags],
                                cwd=case, capture_output=True, text=True, timeout=30,
                                env=dict(os.environ, OMPI_ALLOW_RUN_AS_ROOT="1",
                                         OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1",
                                         OMPI_MCA_rmaps_base_oversubscribe="1"))
        self.assertEqual(result.returncode == 0, success, result.stdout + result.stderr)
        return case, result

    def test_internal_degeneracy_uses_trace_without_auto_enlargement(self):
        energies = [[[-1., -1., 1.], [-1., -1., 1.]]]
        case, result = self.execute_gap_case("trace-degenerate", energies,
                                             ["--task", "kubo", "--bands", "1:2"])
        text = (case / "KUBO.csv").read_text()
        self.assertIn("VASPBERRY_BARE_MOMENTUM_KUBO_BUNDLE_V2", text)
        self.assertEqual(len([line for line in text.splitlines() if line and line[0].isdigit()]), 2)
        self.assertEqual({p.name for p in case.iterdir()}, {"WAVECAR", "KUBO.csv"})
        for i, flags in enumerate((["--bands", "1"], ["--bands", "1:1"],
                                   ["--bands", "1:2", "--per-band", "1"])):
            case, result = self.execute_gap_case(f"degenerate-single-{i}", energies,
                                                 ["--task", "kubo", *flags,
                                                  "--curvature-csv", "BANDS.csv"], False)
            self.assertEqual({p.name for p in case.iterdir()}, {"WAVECAR"})
            self.assertIn("other_band", result.stderr)
            self.assertIn("1e-5 eV", result.stderr)
            self.assertIn("selection was not enlarged", result.stderr)

    def test_later_spin_k_and_nonfinite_gaps_fail_before_any_output(self):
        clear = [[-2., -1., 1.], [-2., -1., 1.]]
        for i, bad in enumerate((-1., -1. + 1.e-5, -1. + 5.e-6, float("nan"))):
            # Third band is outside the selected range; failure is in the last spin/k.
            energies = [clear, [clear[0], [-2., -1., bad]]]
            for j, selection in enumerate((["--bands", "2"],
                                            ["--bands", "1:2", "--per-band", "1"],
                                            ["--bands", "1:2"])):
                case, result = self.execute_gap_case(f"late-gap-{i}-{j}", energies,
                    ["--task", "kubo", *selection, "--curvature-csv", "FAIL.csv"], False)
                self.assertEqual({p.name for p in case.iterdir()}, {"WAVECAR"})
                self.assertIn("isolated", result.stderr)
        case, result = self.execute_gap_case("legacy-degenerate", [[[-1., -1., 1.]] * 2],
                                             ["-kubo", "2", "-ii", "1", "-if", "2"], False)
        self.assertEqual({p.name for p in case.iterdir()}, {"WAVECAR"})

    def test_individual_outputs_and_trace_agree_when_bands_are_isolated(self):
        energies = [[[-2., -1., 1.], [-2., -1., 1.]]]
        trace, _ = self.execute_gap_case("isolated-trace", energies,
                                         ["--task", "kubo", "--bands", "1:2"])
        per_band, _ = self.execute_gap_case("isolated-per-band", energies,
            ["--task", "kubo", "--bands", "1:2", "--per-band", "1", "--curvature-csv", "BANDS.csv"])
        legacy, _ = self.execute_gap_case("isolated-legacy", energies,
            ["-kubo", "2", "-ii", "1", "-if", "2", "-kubo_csv", "BANDS.csv"])
        self.assertEqual((per_band / "BANDS.csv").read_bytes(), (legacy / "BANDS.csv").read_bytes())
        import csv
        def rows(path):
            return list(csv.DictReader(line for line in path.read_text().splitlines() if not line.startswith("#")))
        bands = rows(per_band / "BANDS.csv")
        for row in rows(trace / "KUBO.csv"):
            total = sum(float(b["omega_z_A2"]) for b in bands if b["k_index"] == row["k_index"])
            self.assertAlmostEqual(float(row["omega_z_A2"]), total, delta=1e-12)
        before = (trace / "KUBO.csv").read_bytes()
        result = run_legacy_contract([str(self.binary), "--spinor", "1", "--task", "kubo", "--bands", "1:2"],
                                cwd=trace, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual((trace / "KUBO.csv").read_bytes(), before)

    @unittest.skipUnless(shutil.which("mpifort") and shutil.which("mpiexec"), "MPI tools required")
    def test_mpi_trace_and_gap_rejection_exit_without_partial_files(self):
        wrapper = run_legacy_contract(["mpifort", "--version"], capture_output=True, text=True)
        if "GNU Fortran" not in wrapper.stdout + wrapper.stderr:
            self.skipTest("GNU MPI wrapper required")
        binary = self.work / "vaspberry-mpi"
        built = run_legacy_contract(["mpifort", "-cpp", "-DMPI_USE", "-O0", "-fcheck=all",
                                "-ffixed-line-length-none", "-fallow-argument-mismatch",
                                str(ROOT / "vaspberry.f"), "-llapack", "-lblas", "-o", str(binary)],
                               capture_output=True, text=True)
        self.assertEqual(built.returncode, 0, built.stderr)
        prefix = ["mpiexec", "-n", "2"]
        energies = [[[-1., -1., 1.], [-1., -1., 1.]]]
        serial, _ = self.execute_gap_case("mpi-reference", energies,
                                          ["--task", "kubo", "--bands", "1:2"])
        parallel, _ = self.execute_gap_case("mpi-trace", energies,
                                            ["--task", "kubo", "--bands", "1:2"],
                                            prefix=prefix, binary=binary)
        self.assertEqual((serial / "KUBO.csv").read_bytes(), (parallel / "KUBO.csv").read_bytes())
        energies = [[[-2., -1., 1.], [-2., -1., 1.]],
                    [[-2., -1., 1.], [-2., -1., -1.]]]
        for i, flags in enumerate((["--bands", "2"],
                                   ["--bands", "1:2", "--per-band", "1"],
                                   ["--bands", "1:2"])):
            case, result = self.execute_gap_case(f"mpi-late-gap-{i}", energies,
                ["--task", "kubo", *flags], False, prefix=prefix, binary=binary)
            self.assertEqual({p.name for p in case.iterdir()}, {"WAVECAR"})
            self.assertIn("isolated", result.stderr)

        # The same native MPI executable must resolve word4 and coordinate corrupt-input failure.
        case = self.work / "mpi-word-recl";case.mkdir()
        self.write_gap_fixture(case / "WAVECAR", [[[-1., -1., 1.], [-1., -1., 1.]]])
        self.word_recl_copy(case / "WAVECAR")
        env = dict(os.environ, OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1",
                   OMPI_MCA_rmaps_base_oversubscribe="1")
        command = [*prefix, str(binary), "--spinor", "1", "--task", "kubo", "--bands", "1:2"]
        result = run_legacy_contract(command, cwd=case, capture_output=True, text=True, timeout=30, env=env)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((serial / "KUBO.csv").read_bytes(), (case / "KUBO.csv").read_bytes())
        for corrupt in ("nan", "inf", "truncate"):
            case = self.work / ("mpi-corrupt-" + corrupt);case.mkdir()
            self.write_gap_fixture(case / "WAVECAR", [[[-2., -1., 1.], [-2., -1., 1.]]] * 2)
            data = bytearray((case / "WAVECAR").read_bytes())
            if corrupt == "truncate":
                del data[-1024:]
            else:
                struct.pack_into("<f", data, 17 * 1024, float(corrupt))
            (case / "WAVECAR").write_bytes(data)
            result = run_legacy_contract(command, cwd=case, capture_output=True, text=True, timeout=30, env=env)
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse((case / "KUBO.csv").exists())
            for partial in case.glob("*.csv.partial"):
                self.assertNotIn("result_status=PASS", partial.read_text())

    @staticmethod
    def word_recl_copy(path):
        data = path.read_bytes()
        recl = int(struct.unpack_from("<d", data)[0])
        out = bytearray(4 * len(data))
        for start in range(0, len(data), recl):
            out[4 * start:4 * start + recl] = data[start:start + recl]
        path.write_bytes(out)

    def test_byte_and_word_recl_have_identical_complete_outputs(self):
        energies = [[[-2., -1., 1.], [-2., -1., 1.]]] * 2
        modes = [("trace", ["--task", "kubo", "--bands", "1:2"], "KUBO.csv", 4),
                 ("bands", ["--task", "kubo", "--bands", "1:2", "--per-band", "1",
                             "--curvature-csv", "BANDS.csv"], "BANDS.csv", 8),
                 ("pairs", ["--task", "kubo-pairs", "--pairs-csv", "PAIRS.csv"], "PAIRS.csv", 12)]
        for mode, flags, output, rows in modes:
            results = []
            for word in (False, True):
                case = self.work / f"recl-{mode}-{word}"
                case.mkdir()
                self.write_gap_fixture(case / "WAVECAR", energies)
                if word:
                    self.word_recl_copy(case / "WAVECAR")
                run = run_legacy_contract([str(self.binary), "--spinor", "1", *flags],
                                     cwd=case, capture_output=True, text=True, timeout=30)
                self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
                self.assertIn("RECL_layout=" + ("word4" if word else "byte"), run.stdout)
                text = (case / output).read_text()
                self.assertIn("# source_nspin=2", text)
                self.assertIn("# source_nkpoints=2", text)
                self.assertIn("# expected_rows=" + str(rows), text)
                self.assertEqual(text.rstrip().splitlines()[-1], "# result_status=PASS")
                self.assertEqual(text.count("result_status=PASS"), 1)
                self.assertFalse((case / (output + ".partial")).exists())
                self.assertEqual(len([line for line in text.splitlines() if line and line[0].isdigit()]), rows)
                results.append(text)
            self.assertEqual(*results)

    def test_corrupt_coefficients_never_publish_complete_csv(self):
        clear = [[[-2., -1., 1.], [-2., -1., 1.]]] * 2
        for component, value in ((0, float("nan")), (1, float("inf"))):
            for mode, flags, output in (
                ("trace", ["--task", "kubo", "--bands", "1:2"], "KUBO.csv"),
                ("full", ["--task", "kubo", "--bands", "1:3"], "KUBO.csv"),
                ("bands", ["--task", "kubo", "--bands", "1:2", "--per-band", "1",
                           "--curvature-csv", "BANDS.csv"], "BANDS.csv"),
                ("pairs", ["--task", "kubo-pairs", "--pairs-csv", "PAIRS.csv"], "PAIRS.csv")):
                case = self.work / f"bad-coeff-{component}-{mode}"
                case.mkdir()
                self.write_gap_fixture(case / "WAVECAR", clear)
                data = bytearray((case / "WAVECAR").read_bytes())
                # Last source spin, last k, last band's real/imaginary coefficient.
                struct.pack_into("<f", data, 17 * 1024 + component * 4, value)
                (case / "WAVECAR").write_bytes(data)
                run = run_legacy_contract([str(self.binary), "--spinor", "1", *flags],
                                     cwd=case, capture_output=True, text=True, timeout=30)
                self.assertNotEqual(run.returncode, 0)
                self.assertIn("nonfinite Kubo coefficients", run.stderr)
                self.assertFalse((case / output).exists())
                # Spin 1 already wrote its rows, but these are staged, not complete output.
                self.assertTrue((case / (output + ".partial")).exists())
                for partial in case.glob("*.csv.partial"):
                    self.assertNotIn("result_status=PASS", partial.read_text())

    def test_truncation_and_invalid_or_ambiguous_headers_fail_before_allocation(self):
        clear = [[[-2., -1., 1.], [-2., -1., 1.]]] * 2
        mutations = {
            "truncated": lambda data: data.__delitem__(slice(-1024, None)),
            "fractional-nk": lambda data: struct.pack_into("<d", data, 1024, 2.5),
            "nan-nk": lambda data: struct.pack_into("<d", data, 1024, float("nan")),
            "nan-cutoff": lambda data: struct.pack_into("<d", data, 1040, float("nan")),
            "negative-cutoff": lambda data: struct.pack_into("<d", data, 1040, -1.),
            "singular": lambda data: struct.pack_into("<9d", data, 1048, *([0.] * 9)),
            "huge-record-count": lambda data: struct.pack_into("<d", data, 1024, 2147483647.),
            "fractional-recl": lambda data: struct.pack_into("<d", data, 0, 1024.5),
        }
        for name, mutate in mutations.items():
            case = self.work / ("header-" + name)
            case.mkdir()
            self.write_gap_fixture(case / "WAVECAR", clear)
            data = bytearray((case / "WAVECAR").read_bytes());mutate(data)
            (case / "WAVECAR").write_bytes(data)
            run = run_legacy_contract([str(self.binary), "--task", "kubo", "--bands", "1:2"],
                                 cwd=case, capture_output=True, text=True, timeout=30)
            self.assertNotEqual(run.returncode, 0, run.stdout + run.stderr)
            self.assertNotIn("TOTAL RECORD LENGTH", run.stdout)
            self.assertEqual({p.name for p in case.iterdir()}, {"WAVECAR"})
        case = self.work / "ambiguous-recl";case.mkdir()
        data = bytearray(16384)
        struct.pack_into("<3d", data, 0, 1024., 1., 45200.)
        for stride in (1024, 4096):
            struct.pack_into("<12d", data, stride, 1., 1., 120.,
                             1., 0., 0., 0., 1., 0., 0., 0., 1.)
        (case / "WAVECAR").write_bytes(data)
        run = run_legacy_contract([str(self.binary), "--task", "kubo", "--bands", "1"],
                             cwd=case, capture_output=True, text=True, timeout=30)
        self.assertNotEqual(run.returncode, 0)
        self.assertIn("valid candidates =", run.stderr)
        self.assertNotIn("TOTAL RECORD LENGTH", run.stdout)
        self.assertEqual({p.name for p in case.iterdir()}, {"WAVECAR"})

    def test_one_source_band_still_validates_its_coefficients(self):
        for corrupt in (False, True):
            case = self.work / f"one-band-{corrupt}";case.mkdir()
            data = bytearray(4 * 512)
            struct.pack_into("<3d", data, 0, 512., 1., 45200.)
            struct.pack_into("<12d", data, 512, 1., 1., 120.,
                             1., 0., 0., 0., 1., 0., 0., 0., 1.)
            struct.pack_into("<7d", data, 1024, 1., 0., 0., 0., -1., 0., 2.)
            struct.pack_into("<2f", data, 1536, float("nan") if corrupt else 1., 0.)
            (case / "WAVECAR").write_bytes(data)
            run = run_legacy_contract([str(self.binary), "--spinor", "1", "--task", "kubo", "--bands", "1",
                                  "--curvature-csv", "BAND.csv"],
                                 cwd=case, capture_output=True, text=True, timeout=30)
            self.assertEqual(run.returncode == 0, not corrupt, run.stdout + run.stderr)
            if corrupt:
                self.assertIn("nonfinite Kubo coefficients", run.stderr)
                self.assertFalse((case / "BAND.csv").exists())

    def test_existing_partial_csv_is_preserved(self):
        case = self.work / "existing-partial";case.mkdir()
        self.write_gap_fixture(case / "WAVECAR", [[[-2., -1., 1.], [-2., -1., 1.]]])
        partial = case / "KUBO.csv.partial";partial.write_text("previous interrupted run")
        run = run_legacy_contract([str(self.binary), "--spinor", "1", "--task", "kubo", "--bands", "1:2"],
                             cwd=case, capture_output=True, text=True, timeout=30)
        self.assertNotEqual(run.returncode, 0)
        self.assertFalse((case / "KUBO.csv").exists())
        self.assertEqual(partial.read_text(), "previous interrupted run")

    def test_pairs_ignore_source_occupations_but_subspace_guard_remains(self):
        outputs = []
        for variable in (False, True):
            case = self.work / ("variable" if variable else "uniform")
            case.mkdir()
            self.write_fixture(case / "WAVECAR", variable)
            result = run_legacy_contract([
                str(self.binary), "--task", "kubo-pairs", "--wavecar", "WAVECAR",
                "--spinor", "2", "--pairs-csv", "PAIRS.csv",
            ], cwd=case, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn("Source occupations: not used for pair export", result.stdout)
            outputs.append((case / "PAIRS.csv").read_bytes())
            self.assertTrue(outputs[-1].endswith(b"# result_status=PASS\n"))
        self.assertEqual(outputs[0], outputs[1])
        # Pair independence must not loosen the occupied-subspace validation.
        result = run_legacy_contract([
            str(self.binary), "--task", "chern", "--mesh", "2,1", "--bands", "1",
        ], cwd=self.work / "variable", capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("inconsistent occupations", result.stderr)

    def test_velocity_matches_independent_si_plane_wave_expectation(self):
        for spinor in (1, 2):
            with self.subTest(spinor=spinor):
                self.check_velocity(spinor)

    def check_velocity(self, spinor):
        case = self.work / f"velocity-{spinor}"
        case.mkdir()
        self.write_fixture(case / "WAVECAR", False, spinor_components=spinor)
        result = run_legacy_contract([
            str(self.binary), "--task", "velocity", "--bands", "1",
            "--spinor", str(spinor), "--mesh", "2,1", "-kp", "1",
        ], cwd=case, capture_output=True, text=True)
        # Full executable is compiled with -fcheck=all: this also guards the
        # historical extrema calculation reading ik=nk+1 after the k loop.
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        output = (case / "VEL_EXPT.dat").read_text()
        rows = [list(map(float, line.split())) for line in output.splitlines()
                if line.strip() and not line.startswith("#")]
        self.assertEqual(len(rows), 2)
        # At k=(1/4,1/4,0), the equally weighted G are 0,-x,-y. Thus
        # <kx>=<ky>=-pi/6 A^-1. Use SI constants, independent of the
        # production conversion through electron rest energy in eV.
        hbar_si, mass_si = 1.054571817e-34, 9.1093837015e-31
        amplitude = struct.unpack("<f", struct.pack("<f", 1 / math.sqrt(3 * spinor)))[0]
        expected = hbar_si / mass_si * 1.e10 * (-math.pi / 6) * (3 * spinor * amplitude**2)
        for row in rows:
            target = expected if abs(row[5] - .25) < 1.e-6 else 0.
            self.assertAlmostEqual(row[3], target, delta=5.e-2)
            self.assertAlmostEqual(row[4], target, delta=5.e-2)
        self.assertGreater(abs(expected), 1.e4)
        for axis in ("x", "y"):
            min_line = next(line for line in output.splitlines()
                            if f"MINVAL of VEL_EXPT <v_{axis}>" in line)
            max_line = next(line for line in output.splitlines()
                            if f"MAXVAL of VEL_EXPT <v_{axis}>" in line)
            self.assertAlmostEqual(float(min_line.split()[-1]), expected, delta=5.e-2)
            self.assertAlmostEqual(float(max_line.split()[-1]), 0., delta=1.e-6)


if __name__ == "__main__":
    unittest.main()
