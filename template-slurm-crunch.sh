#!/bin/bash

#SBATCH --job-name {{ momentum|join('') }}-{{ irrep }}
#SBATCH --cpus-per-task=1
#SBATCH --mem=500MB
#SBATCH --mail-user=ueding@hiskp.uni-bonn.de
#SBATCH --mail-type=FAIL

set -e
set -u
set -x

hostname
date -Iseconds

./number_crunching.R {{ momentum|join(' ') }} {{ irrep }}

date -Iseconds
