#!/bin/bash

#SBATCH --job-name projected_merge
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=50GB
#SBATCH --mail-user=fischer@hiskp.uni-bonn.de
#SBATCH --mail-type=FAIL

#SBATCH --output=batch_output/projected_merge_slurm_%j.txt

/usr/bin/time Rscript -e 'numericprojection::projected_merge()'
