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

FC = gfortran
MPIFC = mpifort
MPIEXEC = mpiexec
IFX = ifx
MPIIFX = mpiifx
IFORT = ifort
MPIIFORT = mpiifort
INTEL_MPIEXEC = mpiexec.hydra

GNU_FLAGS ?= -cpp -O2 -ffixed-line-length-none -fallow-argument-mismatch
GNU_LIBS ?= -llapack -lblas
INTEL_FLAGS ?= -fpp -O2 -extend-source -assume byterecl
IFX_MKL_FLAGS ?= -qmkl=sequential
IFORT_MKL_FLAGS ?= -mkl=sequential
MPIEXEC_FLAGS ?=
INTEL_MPIEXEC_FLAGS ?=

SERIAL_SOURCE := vaspberry.f
# Historical no-Kubo source is available explicitly via
# make serial SERIAL_SOURCE=vaspberry_gfortran_serial.f BUILD_DIR=build-legacy
MPI_SOURCE := vaspberry.f
GNU_SERIAL_BIN := $(BUILD_DIR)/vaspberry
GNU_SERIAL_COMPAT_BIN := $(BUILD_DIR)/vaspberry-gfortran
GNU_MPI_BIN := $(BUILD_DIR)/vaspberry-mpi
MPI_RUNTIME_TEST := $(BUILD_DIR)/test-mpi-runtime

.PHONY: all gnu serial mpi check check-gnu check-serial-help \
	check-mpi-help check-mpi-runtime check-build-dir ifx ifx-mpi \
	ifort ifort-mpi check-ifx-mpi check-ifort-mpi clean force-serial-compat

all: gnu

gnu: serial mpi

serial: $(GNU_SERIAL_BIN) $(GNU_SERIAL_COMPAT_BIN)

mpi: $(GNU_MPI_BIN)

check-build-dir:
	@if [ "$(words $(BUILD_DIR_ABS))" -ne 1 ] || \
	    [ -z "$(SAFE_BUILD_DIR)" ]; then \
	  echo "error: BUILD_DIR must be build or a build-* directory at repository root" >&2; \
	  exit 2; \
	fi

$(BUILD_DIR): | check-build-dir
	mkdir -p -- "$@"

$(GNU_SERIAL_BIN): $(SERIAL_SOURCE) | $(BUILD_DIR)
	$(FC) $(GNU_FLAGS) -o $@ $< $(GNU_LIBS)

# Keep the old compiler-specific name without a second compiled binary.
force-serial-compat:

$(GNU_SERIAL_COMPAT_BIN): $(GNU_SERIAL_BIN) force-serial-compat
	ln -sf vaspberry "$@"

$(GNU_MPI_BIN): $(MPI_SOURCE) | $(BUILD_DIR)
	$(MPIFC) $(GNU_FLAGS) -DMPI_USE -o $@ $< $(GNU_LIBS)

$(MPI_RUNTIME_TEST): tests/fortran/test_mpi_runtime.f90 | $(BUILD_DIR)
	$(MPIFC) -O2 -o $@ $<

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
	$(IFX) $(INTEL_FLAGS) -o $(BUILD_DIR)/vaspberry-ifx \
	  $(SERIAL_SOURCE) $(IFX_MKL_FLAGS)

ifx-mpi: | $(BUILD_DIR)
	$(MPIIFX) $(INTEL_FLAGS) -DMPI_USE -o $(BUILD_DIR)/vaspberry-ifx-mpi \
	  $(MPI_SOURCE) $(IFX_MKL_FLAGS)

ifort: | $(BUILD_DIR)
	$(IFORT) $(INTEL_FLAGS) -o $(BUILD_DIR)/vaspberry-ifort \
	  $(SERIAL_SOURCE) $(IFORT_MKL_FLAGS)

ifort-mpi: | $(BUILD_DIR)
	$(MPIIFORT) $(INTEL_FLAGS) -DMPI_USE \
	  -o $(BUILD_DIR)/vaspberry-ifort-mpi \
	  $(MPI_SOURCE) $(IFORT_MKL_FLAGS)

check-ifx-mpi: ifx-mpi
	$(MPIIFX) -O2 -o $(BUILD_DIR)/test-ifx-mpi-runtime tests/fortran/test_mpi_runtime.f90
	$(INTEL_MPIEXEC) $(INTEL_MPIEXEC_FLAGS) -n 2 $(BUILD_DIR)/test-ifx-mpi-runtime
	$(INTEL_MPIEXEC) $(INTEL_MPIEXEC_FLAGS) -n 2 $(BUILD_DIR)/vaspberry-ifx-mpi --help > $(BUILD_DIR)/help-ifx-mpi.txt
	sh tests/check_fortran_help.sh $(BUILD_DIR)/help-ifx-mpi.txt

check-ifort-mpi: ifort-mpi
	$(MPIIFORT) -O2 -o $(BUILD_DIR)/test-ifort-mpi-runtime tests/fortran/test_mpi_runtime.f90
	$(INTEL_MPIEXEC) $(INTEL_MPIEXEC_FLAGS) -n 2 $(BUILD_DIR)/test-ifort-mpi-runtime
	$(INTEL_MPIEXEC) $(INTEL_MPIEXEC_FLAGS) -n 2 $(BUILD_DIR)/vaspberry-ifort-mpi --help > $(BUILD_DIR)/help-ifort-mpi.txt
	sh tests/check_fortran_help.sh $(BUILD_DIR)/help-ifort-mpi.txt

clean: check-build-dir
	rm -rf -- "$(SAFE_BUILD_DIR)"
