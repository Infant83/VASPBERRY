"""Exercise the build contract without requiring a particular Fortran compiler."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(shutil.which("make"), "GNU make required")
class MakefileTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="vaspberry-make-")
        self.addCleanup(self.temporary.cleanup)
        self.project = Path(self.temporary.name)
        shutil.copy2(ROOT / "Makefile", self.project / "Makefile")
        (self.project / "vaspberry.f").write_text("      end\n")
        (self.project / "legacy.f").write_text("      end\n")
        os.utime(self.project / "legacy.f", (1, 1))
        (self.project / "tests/fortran").mkdir(parents=True)
        (self.project / "tests/fortran/test_mpi_runtime.f90").write_text("end\n")
        self.calls = self.project / "calls.jsonl"
        self.compilers = []
        for name in ("compiler-a", "compiler-b"):
            path = self.project / name
            path.write_text(
                "#!" + sys.executable + "\n"
                "import json, os, pathlib, sys\n"
                "with open(os.environ['COMPILER_CALLS'], 'a') as stream:\n"
                "    stream.write(json.dumps(sys.argv) + '\\n')\n"
                "output = pathlib.Path(sys.argv[sys.argv.index('-o') + 1])\n"
                "output.write_text('#!/bin/sh\\nexit 0\\n')\n"
                "output.chmod(0o755)\n"
                "sys.exit(int(os.environ.get('COMPILER_EXIT', '0')))\n"
            )
            path.chmod(0o755)
            self.compilers.append(str(path))
        self.environment = dict(os.environ, COMPILER_CALLS=str(self.calls))
        for variable in ("FC", "MPIFC", "MPIEXEC", "IFX", "IFORT", "MPIIFX", "MPIIFORT",
                         "INTEL_MPIEXEC", "FFLAGS", "LDFLAGS", "MAKEFLAGS", "MFLAGS",
                         "OMPI_FC", "OMPI_FCFLAGS", "OMPI_LDFLAGS", "OMPI_LIBS", "MPICH_FC",
                         "GNU_FLAGS", "GNU_LIBS", "INTEL_FLAGS", "IFX_MKL_FLAGS", "IFORT_MKL_FLAGS"):
            self.environment.pop(variable, None)

    def make(self, *arguments, environment=None, success=True):
        result = subprocess.run(["make", *arguments], cwd=self.project,
                                env=environment or self.environment, text=True, capture_output=True)
        if success:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def invocations(self):
        return [json.loads(line) for line in self.calls.read_text().splitlines()] if self.calls.exists() else []

    def test_default_gnu_compiler_replaces_make_builtin_f77(self):
        result = self.make("-n", "serial")
        self.assertIn("gfortran -cpp -O2 -ffixed-line-length-none -fallow-argument-mismatch", result.stdout)

    def test_environment_tool_overrides_and_command_line_precedence(self):
        replacements = {name: "site-" + name.lower() for name in
                        ("FC", "MPIFC", "MPIEXEC", "IFX", "IFORT", "MPIIFX", "MPIIFORT", "INTEL_MPIEXEC")}
        env = dict(self.environment, **replacements)
        result = self.make("-n", "check-gnu", "ifx", "ifort", "check-ifx-mpi", "check-ifort-mpi",
                           environment=env)
        for command in replacements.values():
            self.assertIn(command + " ", result.stdout)
        result = self.make("-n", "serial", "FC=explicit-fc", environment=env)
        self.assertIn("explicit-fc -cpp", result.stdout)
        self.assertNotIn("site-fc -cpp", result.stdout)

    def test_optional_flags_are_additive_and_linker_flags_precede_libraries(self):
        env = dict(self.environment, FFLAGS="-g", LDFLAGS="-L/site/lib")
        self.make("serial", "FC=" + self.compilers[0], environment=env)
        args = self.invocations()[0]
        for flag in ("-cpp", "-ffixed-line-length-none", "-fallow-argument-mismatch", "-g", "-L/site/lib"):
            self.assertIn(flag, args)
        self.assertLess(args.index("-L/site/lib"), args.index("-llapack"))

    def test_serial_rebuilds_for_configuration_and_older_source_changes(self):
        common = ["serial", "FC=" + self.compilers[0]]
        self.make(*common)
        self.make(*common)
        self.assertEqual(len(self.invocations()), 1)
        variants = [
            [*common, "GNU_FLAGS=-cpp -O0 -ffixed-line-length-none"],
            [*common, "FFLAGS=-g"],
            [*common, "LDFLAGS=-L/changed"],
            [*common, "GNU_LIBS=-lopenblas"],
            ["serial", "FC=" + self.compilers[1]],
            [*common, "SERIAL_SOURCE=legacy.f"],
        ]
        for count, arguments in enumerate(variants, start=2):
            with self.subTest(arguments=arguments):
                self.make(*arguments)
                self.assertEqual(len(self.invocations()), count)
                self.make(*arguments)
                self.assertEqual(len(self.invocations()), count)
        self.assertIn("legacy.f", self.invocations()[-1])

    def test_makefile_change_rebuilds_binary(self):
        args = ("serial", "FC=" + self.compilers[0])
        self.make(*args)
        makefile = self.project / "Makefile"
        makefile.write_text(makefile.read_text() + "\n# Updated build instructions\n")
        os.utime(makefile, (1, 1))
        self.make(*args)
        self.assertEqual(len(self.invocations()), 2)

    def test_parallel_gnu_build_and_runtime_probe_track_wrapper_changes(self):
        args = ["-j4", "gnu", "build/test-mpi-runtime", "FC=" + self.compilers[0],
                "MPIFC=" + self.compilers[0]]
        self.make(*args)
        self.assertEqual(len(self.invocations()), 3)
        alias = self.project / "build/vaspberry-gfortran"
        self.assertEqual(alias.resolve(), (self.project / "build/vaspberry").resolve())
        alias.unlink()
        alias.write_text("old independent executable")
        self.make(*args)
        self.assertEqual(len(self.invocations()), 3)
        self.assertTrue(alias.is_symlink())
        self.make(*args[:-1], "MPIFC=" + self.compilers[1])
        self.assertEqual(len(self.invocations()), 5)
        self.assertTrue(all(row[0] == self.compilers[1] for row in self.invocations()[-2:]))

    def test_failed_compile_is_not_reused_as_success(self):
        args = ("serial", "FC=" + self.compilers[0])
        self.make(*args, environment=dict(self.environment, COMPILER_EXIT="1"), success=False)
        self.make(*args)
        self.assertEqual(len(self.invocations()), 2)

    def test_mpi_wrapper_environment_changes_rebuild_binary_and_probe(self):
        args = ("-j2", "mpi", "build/test-mpi-runtime", "MPIFC=" + self.compilers[0])
        self.make(*args)
        self.assertEqual(len(self.invocations()), 2)
        env = self.environment.copy()
        for count, variable in enumerate(("OMPI_FC", "MPICH_FC", "OMPI_FCFLAGS", "OMPI_LDFLAGS", "OMPI_LIBS"), 2):
            env[variable] = "changed-wrapper-setting"
            self.make(*args, environment=env)
            self.assertEqual(len(self.invocations()), count * 2)
            self.make(*args, environment=env)
            self.assertEqual(len(self.invocations()), count * 2)

    def test_changed_source_with_old_timestamp_and_failed_rebuild_are_detected(self):
        args = ("serial", "FC=" + self.compilers[0])
        self.make(*args)
        source = self.project / "vaspberry.f"
        source.write_text("      program changed\n      end\n")
        os.utime(source, (1, 1))
        self.make(*args, environment=dict(self.environment, COMPILER_EXIT="1"), success=False)
        self.assertFalse((self.project / "build/vaspberry").exists())
        self.make(*args)
        self.assertEqual(len(self.invocations()), 3)
        self.make(*args)
        self.assertEqual(len(self.invocations()), 3)

    def test_intel_serial_checks_and_help_require_no_gnu_toolchain(self):
        for compiler in ("ifx", "ifort"):
            result = self.make("-n", "check-" + compiler)
            self.assertIn("build/vaspberry-" + compiler + " --help", result.stdout)
            self.assertNotIn("gfortran -", result.stdout)
        result = self.make("help")
        self.assertIn("ifx-mpi", result.stdout)
        self.assertIn("FFLAGS", result.stdout)

    def test_symlink_build_directory_is_rejected_without_touching_destination(self):
        outside = self.project / "external"
        outside.mkdir()
        sentinel = outside / "keep.txt"
        sentinel.write_text("keep\n")
        (self.project / "build-external").symlink_to(outside, target_is_directory=True)
        for target in ("serial", "clean"):
            result = self.make(target, "BUILD_DIR=build-external", "FC=" + self.compilers[0], success=False)
            self.assertIn("symlink build directories are not supported", result.stderr)
            self.assertEqual(sentinel.read_text(), "keep\n")
        self.assertEqual(self.invocations(), [])


if __name__ == "__main__":
    unittest.main()
