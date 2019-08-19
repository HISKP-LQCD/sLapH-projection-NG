#!/bin/bash

#SBATCH --job-name N_{{ momentum|join('') }}-{{ irrep }}-{{ '%04d'|format(config_number) }}
#SBATCH --cpus-per-task=1
#SBATCH --mem=500MB
#SBATCH --mail-user=ueding@hiskp.uni-bonn.de
#SBATCH --mail-type=FAIL

#SBATCH --output=batch_output/numeric-{{ momentum|join('') }}-{{ irrep }}-slurm_%j.txt

set -e
set -u
set -x

hostname
date -Iseconds

/usr/bin/time {{ srcdir }}/numeric_projection/number_crunching.R {{ momentum|join(' ') }} {{ irrep }} {{ config_number }}

date -Iseconds
