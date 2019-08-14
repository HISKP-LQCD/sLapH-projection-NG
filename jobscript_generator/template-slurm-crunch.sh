#!/bin/bash

#SBATCH --job-name Cr{{ momentum|join('') }}-{{ irrep }}
#SBATCH --cpus-per-task=1
#SBATCH --mem=500MB
#SBATCH --mail-user=ueding@hiskp.uni-bonn.de
#SBATCH --mail-type=FAIL

#SBATCH --output=batch_output/crunch-{{ momentum|join('') }}-{{ irrep }}-slurm_%j.txt

set -e
set -u
set -x

hostname
date -Iseconds

/usr/bin/time {{ srcdir }}/numeric_projection/number_crunching.R {{ momentum|join(' ') }} {{ irrep }}

date -Iseconds
