SHELL := /bin/sh

BUILD_DIR ?= build
REPO_ROOT := $(abspath .)
BUILD_DIR_ABS := $(abspath $(BUILD_DIR))
BUILD_DIR_PARENT := $(patsubst %/,%,$(dir $(BUILD_DIR_ABS)))
BUILD_DIR_NAME := $(notdir $(BUILD_DIR_ABS))
ifeq ($(words $(BUILD_DIR_ABS)),1)
ifeq ($(BUILD_DIR_PARENT),$(REPO_ROOT))
ifneq ($(filter build build-%,$(BUILD_DIR_NAME)),)
SAFE_BUILD_DIR := $(BUILD_DIR_ABS)
endif
endif
endif

# GNU make predefines FC=f77. Replace only that built-in (or an unset FC),
# while preserving an exported site/compiler choice and command-line overrides.
ifneq ($(filter default undefined,$(origin FC)),)
FC := gfortran
endif
MPIFC ?= mpifort
MPIEXEC ?= mpiexec
IFX ?= ifx
MPIIFX ?= mpiifx
IFORT ?= ifort
MPIIFORT ?= mpiifort
INTEL_MPIEXEC ?= mpiexec.hydra

GNU_FLAGS ?= -cpp -O2 -ffixed-line-length-none -fallow-argument-mismatch
GNU_LIBS ?= -llapack -lblas
INTEL_FLAGS ?= -fpp -O2 -extend-source -assume byterecl
IFX_MKL_FLAGS ?= -qmkl=sequential
IFORT_MKL_FLAGS ?= -mkl=sequential
MPIEXEC_FLAGS ?=
INTEL_MPIEXEC_FLAGS ?=
# Optional additions; the compiler-specific defaults above remain in effect.
FFLAGS ?=
LDFLAGS ?=

SERIAL_SOURCE := vaspberry.f
# Historical no-Kubo source is available explicitly via
# make serial SERIAL_SOURCE=vaspberry_gfortran_serial.f BUILD_DIR=build-legacy
MPI_SOURCE := vaspberry.f
GNU_SERIAL_BIN := $(BUILD_DIR)/vaspberry
GNU_SERIAL_COMPAT_BIN := $(BUILD_DIR)/vaspberry-gfortran
GNU_MPI_BIN := $(BUILD_DIR)/vaspberry-mpi
MPI_RUNTIME_TEST := $(BUILD_DIR)/test-mpi-runtime
SERIAL_CONFIG := $(BUILD_DIR)/.serial-config
MPI_CONFIG := $(BUILD_DIR)/.mpi-config
MPI_RUNTIME_CONFIG := $(BUILD_DIR)/.mpi-runtime-config

# Compare full commands on every invocation and rebuild only if configuration or
# prerequisites changed. Avoid stamp-mtime comparisons: make 3.81 may lose
# subsecond precision. A successful build records its command, source/Makefile
# checksums and PATH. Use clean or a new BUILD_DIR after replacing a compiler or
# library in place while retaining the same command, flags and search paths.
SERIAL_COMMAND = $(FC) $(GNU_FLAGS) $(FFLAGS) -o $(GNU_SERIAL_BIN) $(SERIAL_SOURCE) $(LDFLAGS) $(GNU_LIBS)
MPI_COMMAND = $(MPIFC) $(GNU_FLAGS) $(FFLAGS) -DMPI_USE -o $(GNU_MPI_BIN) $(MPI_SOURCE) $(LDFLAGS) $(GNU_LIBS)
MPI_RUNTIME_COMMAND = $(MPIFC) -O2 $(FFLAGS) -o $(MPI_RUNTIME_TEST) tests/fortran/test_mpi_runtime.f90 $(LDFLAGS)
shell_quote = '$(subst ','"'"',$(1))'
define build_if_changed
@set -e; tmp="$(2).tmp.$$$$"; trap 'rm -f "$$tmp"' 0 1 2 3 15; \
  printf '%s\n' $(call shell_quote,$(1)) $(call shell_quote,PATH=$(PATH)) \
    $(call shell_quote,LIBRARY_PATH=$(LIBRARY_PATH)) $(call shell_quote,CPATH=$(CPATH)) \
    $(call shell_quote,OMPI_FC=$(OMPI_FC)) $(call shell_quote,OMPI_FCFLAGS=$(OMPI_FCFLAGS)) \
    $(call shell_quote,OMPI_LDFLAGS=$(OMPI_LDFLAGS)) $(call shell_quote,OMPI_LIBS=$(OMPI_LIBS)) \
    $(call shell_quote,MPICH_FC=$(MPICH_FC)) > "$$tmp"; \
  cksum $(call shell_quote,$<) $(foreach file,$(MAKEFILE_LIST),$(call shell_quote,$(file))) >> "$$tmp"; \
  if [ ! -f "$@" ] || [ -n $(call shell_quote,$(filter-out force-build-config,$?)) ] || \
      ! cmp -s "$$tmp" "$(2)"; then \
    printf '%s\n' $(call shell_quote,$(1)); \
    $(1) || { status=$$?; rm -f "$@"; exit "$$status"; }; \
    mv -f "$$tmp" "$(2)"; \
  fi
endef

.DELETE_ON_ERROR:
.PHONY: all gnu serial mpi help check check-gnu check-serial-help \
	check-mpi-help check-mpi-runtime check-build-dir ifx ifx-mpi \
	ifort ifort-mpi check-ifx check-ifort check-ifx-mpi check-ifort-mpi \
	clean force-serial-compat force-build-config

all: gnu

gnu: serial mpi

serial: $(GNU_SERIAL_BIN) $(GNU_SERIAL_COMPAT_BIN)

mpi: $(GNU_MPI_BIN)

help:
	@printf '%s\n' \
	  'VASPBERRY native Fortran build (no Python required)' \
	  '  make serial / mpi / gnu       GNU serial / MPI / both (default: gnu)' \
	  '  make ifx / ifx-mpi            Intel oneAPI serial / Intel MPI' \
	  '  make ifort / ifort-mpi        Retained Intel classic compiler' \
	  '  make check-serial-help       GNU serial executable/help check' \
	  '  make check-gnu               GNU serial + two-rank MPI checks' \
	  '  make check-ifx[-mpi]         Intel ifx serial or two-rank MPI check' \
	  '  make check-ifort[-mpi]       Intel ifort serial or two-rank MPI check' \
	  '  make clean                   Remove only the selected build directory' \
	  'Tools: FC MPIFC MPIEXEC IFX MPIIFX IFORT MPIIFORT INTEL_MPIEXEC' \
	  'Flags: FFLAGS LDFLAGS add to compiler-specific defaults; GNU_LIBS sets BLAS/LAPACK' \
	  'Overrides: export variables or pass NAME=value to make; BUILD_DIR=build or build-*' \
	  'Install prerequisites and compiler/runtime matching: docs/BUILD.md'

check-build-dir:
	@if [ "$(words $(BUILD_DIR_ABS))" -ne 1 ] || \
	    [ -z "$(SAFE_BUILD_DIR)" ] || [ -L "$(BUILD_DIR_ABS)" ]; then \
	  echo "error: BUILD_DIR must be build or a build-* directory at repository root" >&2; \
	  echo "       symlink build directories are not supported" >&2; \
	  exit 2; \
	fi

$(BUILD_DIR): | check-build-dir
	mkdir -p -- "$@"

force-build-config:

$(GNU_SERIAL_BIN): $(SERIAL_SOURCE) $(MAKEFILE_LIST) force-build-config | $(BUILD_DIR)
	$(call build_if_changed,$(SERIAL_COMMAND),$(SERIAL_CONFIG))

# Keep the old compiler-specific name without a second compiled binary.
force-serial-compat:

$(GNU_SERIAL_COMPAT_BIN): $(GNU_SERIAL_BIN) force-serial-compat
	ln -sf vaspberry "$@"

$(GNU_MPI_BIN): $(MPI_SOURCE) $(MAKEFILE_LIST) force-build-config | $(BUILD_DIR)
	$(call build_if_changed,$(MPI_COMMAND),$(MPI_CONFIG))

$(MPI_RUNTIME_TEST): tests/fortran/test_mpi_runtime.f90 $(MAKEFILE_LIST) force-build-config | $(BUILD_DIR)
	$(call build_if_changed,$(MPI_RUNTIME_COMMAND),$(MPI_RUNTIME_CONFIG))

check: check-gnu

check-gnu: check-serial-help check-mpi-runtime check-mpi-help

check-serial-help: $(GNU_SERIAL_BIN)
	$< --help > $(BUILD_DIR)/help-serial.txt
	sh tests/check_fortran_help.sh $(BUILD_DIR)/help-serial.txt

check-mpi-runtime: $(MPI_RUNTIME_TEST)
	$(MPIEXEC) $(MPIEXEC_FLAGS) -n 2 $<

check-mpi-help: $(GNU_MPI_BIN)
	$(MPIEXEC) $(MPIEXEC_FLAGS) -n 2 $< -h > $(BUILD_DIR)/help-mpi.txt
	sh tests/check_fortran_help.sh $(BUILD_DIR)/help-mpi.txt

# Load the Intel compiler, MPI SDK and oneMKL environment before these targets.
# The dedicated Intel MPI workflow exercises the default oneMKL link flags.
ifx: | $(BUILD_DIR)
	$(IFX) $(INTEL_FLAGS) $(FFLAGS) -o $(BUILD_DIR)/vaspberry-ifx \
	  $(SERIAL_SOURCE) $(LDFLAGS) $(IFX_MKL_FLAGS)

ifx-mpi: | $(BUILD_DIR)
	$(MPIIFX) $(INTEL_FLAGS) -DMPI_USE $(FFLAGS) -o $(BUILD_DIR)/vaspberry-ifx-mpi \
	  $(MPI_SOURCE) $(LDFLAGS) $(IFX_MKL_FLAGS)

ifort: | $(BUILD_DIR)
	$(IFORT) $(INTEL_FLAGS) $(FFLAGS) -o $(BUILD_DIR)/vaspberry-ifort \
	  $(SERIAL_SOURCE) $(LDFLAGS) $(IFORT_MKL_FLAGS)

ifort-mpi: | $(BUILD_DIR)
	$(MPIIFORT) $(INTEL_FLAGS) -DMPI_USE $(FFLAGS) \
	  -o $(BUILD_DIR)/vaspberry-ifort-mpi \
	  $(MPI_SOURCE) $(LDFLAGS) $(IFORT_MKL_FLAGS)

check-ifx: ifx
	$(BUILD_DIR)/vaspberry-ifx --help > $(BUILD_DIR)/help-ifx.txt
	sh tests/check_fortran_help.sh $(BUILD_DIR)/help-ifx.txt

check-ifort: ifort
	$(BUILD_DIR)/vaspberry-ifort --help > $(BUILD_DIR)/help-ifort.txt
	sh tests/check_fortran_help.sh $(BUILD_DIR)/help-ifort.txt

check-ifx-mpi: ifx-mpi
	$(MPIIFX) -O2 -o $(BUILD_DIR)/test-ifx-mpi-runtime tests/fortran/test_mpi_runtime.f90 $(FFLAGS) $(LDFLAGS)
	$(INTEL_MPIEXEC) $(INTEL_MPIEXEC_FLAGS) -n 2 $(BUILD_DIR)/test-ifx-mpi-runtime
	$(INTEL_MPIEXEC) $(INTEL_MPIEXEC_FLAGS) -n 2 $(BUILD_DIR)/vaspberry-ifx-mpi --help > $(BUILD_DIR)/help-ifx-mpi.txt
	sh tests/check_fortran_help.sh $(BUILD_DIR)/help-ifx-mpi.txt

check-ifort-mpi: ifort-mpi
	$(MPIIFORT) -O2 -o $(BUILD_DIR)/test-ifort-mpi-runtime tests/fortran/test_mpi_runtime.f90 $(FFLAGS) $(LDFLAGS)
	$(INTEL_MPIEXEC) $(INTEL_MPIEXEC_FLAGS) -n 2 $(BUILD_DIR)/test-ifort-mpi-runtime
	$(INTEL_MPIEXEC) $(INTEL_MPIEXEC_FLAGS) -n 2 $(BUILD_DIR)/vaspberry-ifort-mpi --help > $(BUILD_DIR)/help-ifort-mpi.txt
	sh tests/check_fortran_help.sh $(BUILD_DIR)/help-ifort-mpi.txt

clean: check-build-dir
	rm -rf -- "$(SAFE_BUILD_DIR)"
