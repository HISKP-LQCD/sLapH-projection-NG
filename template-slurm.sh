#!/bin/bash

#SBATCH --job-name {{ momentum|join('') }}-{{ irrep }}
#SBATCH --cpus-per-task=1
#SBATCH --mem=500MB
#SBATCH --mail-user=ueding@hiskp.uni-bonn.de
#SBATCH --mail-type=FAIL

set -e
set -u
set -x

wolframscript=/usr/remote/Wolfram/Mathematica/11.3/SystemFiles/Kernel/Binaries/Linux-x86-64/wolframscript
if ! [[ -f "$wolframscript" ]]; then
    wolframscript=/usr/local/Wolfram/Mathematica/11.3/SystemFiles/Kernel/Binaries/Linux-x86-64/wolframscript
fi

date -Iseconds

$wolframscript -script driver.wl {{ momentum|join(' ') }} {{ irrep }}

date -Iseconds
