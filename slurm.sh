#!/bin/bash

#SBATCH --job-name 001_A1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8GB
#SBATCH --mail-user=ueding@hiskp.uni-bonn.de
#SBATCH --mail-type=ALL

set -e
set -u
set -x

wolframscript=/usr/remote/Wolfram/Mathematica/11.3/SystemFiles/Kernel/Binaries/Linux-x86-64/wolframscript
if ! [[ -f "$wolframscript" ]]; then
    wolframscript=/usr/local/Wolfram/Mathematica/11.3/SystemFiles/Kernel/Binaries/Linux-x86-64/wolframscript
fi

date -Iseconds

$wolframscript -script driver.wl 1 0 0 A1

date -Iseconds
