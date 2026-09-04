#!/usr/bin/env bash
#SBATCH --job-name=mos2-scf
#SBATCH --output=%x.o%j
#SBATCH --error=%x.e%j
#SBATCH --nodes=1
#SBATCH --time=24:00:00

set -euo pipefail

: "${VASP_BIN:?Set VASP_BIN to your licensed VASP executable}"

# Add site-specific partition, account, task count, and module commands above.
# This template deliberately contains no private paths, email addresses, or
# scratch-cleanup commands.
srun "$VASP_BIN" > log.out
