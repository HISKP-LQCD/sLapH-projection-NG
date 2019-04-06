#!/bin/bash

#SBATCH --job-name Proj1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8GB
#SBATCH --mail-user=ueding@hiskp.uni-bonn.de
#SBATCH --mail-type=ALL

set -e
set -u
set -x

date -I seconds

math -script Projection_for_rho.wls

date -I seconds
