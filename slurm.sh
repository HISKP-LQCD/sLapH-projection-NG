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

date -I seconds

$wolframscript -script driver.wls "{0, 0, 1}" "A1"

date -I seconds
