#!/bin/bash

#SBATCH --job-name A_{{ momentum|join('') }}_{{ irrep }}
#SBATCH --cpus-per-task=1
#SBATCH --mem=500MB
#SBATCH --mail-user=ueding@hiskp.uni-bonn.de
#SBATCH --mail-type=FAIL

#SBATCH --output=batch_output/analytic_{{ momentum|join('') }}_{{ irrep }}_slurm_%j.txt

set -e
set -u
set -x

wolframscript=/usr/remote/Wolfram/Mathematica/11.3/SystemFiles/Kernel/Binaries/Linux-x86-64/wolframscript
if ! [[ -f "$wolframscript" ]]; then
    wolframscript=/usr/local/Wolfram/Mathematica/12.0/SystemFiles/Kernel/Binaries/Linux-x86-64/wolframscript
fi

hostname
date -Iseconds

mkdir -p prescriptions
time $wolframscript -script {{ srcdir }}/analytic_projection/driver_rho.wl {{ momentum|join(' ') }} {{ irrep }}

date -Iseconds
