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

GNU_FLAGS ?= -cpp -O2 -ffixed-line-length-none -fallow-argument-mismatch
GNU_LIBS ?= -llapack -lblas
INTEL_FLAGS ?= -fpp -O2 -extend-source -assume byterecl
IFX_MKL_FLAGS ?= -qmkl=sequential
IFORT_MKL_FLAGS ?= -mkl=sequential
MPIEXEC_FLAGS ?=

SERIAL_SOURCE := vaspberry_gfortran_serial.f
MPI_SOURCE := vaspberry.f
GNU_SERIAL_BIN := $(BUILD_DIR)/vaspberry-gfortran
GNU_MPI_BIN := $(BUILD_DIR)/vaspberry-mpi
MPI_RUNTIME_TEST := $(BUILD_DIR)/test-mpi-runtime

.PHONY: all gnu serial mpi check check-gnu check-serial-help \
	check-mpi-help check-mpi-runtime check-build-dir ifx ifx-mpi \
	ifort ifort-mpi clean

all: gnu

gnu: serial mpi

serial: $(GNU_SERIAL_BIN)

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

$(GNU_MPI_BIN): $(MPI_SOURCE) | $(BUILD_DIR)
	$(MPIFC) $(GNU_FLAGS) -DMPI_USE -o $@ $< $(GNU_LIBS)

$(MPI_RUNTIME_TEST): tests/fortran/test_mpi_runtime.f90 | $(BUILD_DIR)
	$(MPIFC) -O2 -o $@ $<

check: check-gnu

check-gnu: check-serial-help check-mpi-runtime check-mpi-help

check-serial-help: $(GNU_SERIAL_BIN)
	$< -h > $(BUILD_DIR)/help-serial.txt
	sh tests/check_fortran_help.sh $(BUILD_DIR)/help-serial.txt

check-mpi-runtime: $(MPI_RUNTIME_TEST)
	$(MPIEXEC) $(MPIEXEC_FLAGS) -n 2 $<

check-mpi-help: $(GNU_MPI_BIN)
	$(MPIEXEC) $(MPIEXEC_FLAGS) -n 2 $< -h > $(BUILD_DIR)/help-mpi.txt
	sh tests/check_fortran_help.sh $(BUILD_DIR)/help-mpi.txt

# Intel MPI and oneMKL recipes remain opt-in site targets. CI provisions the
# serial compilers separately and overrides the link flags with system LP64
# BLAS/LAPACK.
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

clean: check-build-dir
	rm -rf -- "$(SAFE_BUILD_DIR)"
