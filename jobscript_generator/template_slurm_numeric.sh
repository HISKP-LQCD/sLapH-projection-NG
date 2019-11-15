#!/bin/bash

#SBATCH --job-name N_{{ '%04d'|format(config_number) }}
#SBATCH --cpus-per-task=1
#SBATCH --mem=500MB
#SBATCH --mail-user=ueding@hiskp.uni-bonn.de
#SBATCH --mail-type=FAIL

#SBATCH --output=batch_output/numeric_{{ config_number }}_slurm_%j.txt

set -e
set -u
set -x

hostname
date -Iseconds

/usr/bin/time {{ srcdir }}/numeric_projection/driver.R {{ config_number }}

date -Iseconds
