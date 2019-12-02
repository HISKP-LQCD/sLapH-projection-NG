#!/bin/bash

#SBATCH --job-name projected_merge
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=950MB
#SBATCH --mail-user=ueding@hiskp.uni-bonn.de
#SBATCH --mail-type=FAIL

#SBATCH --output=batch_output/projected_merge_slurm_%j.txt

r_libs_user="$(Rscript -e "cat(Sys.getenv('R_LIBS_USER'))")"
script="$r_libs_user/numericprojection/exec/projected_merge.R"
script="${script/#\~/$HOME}"  # https://stackoverflow.com/a/27485157/653152

"$script"
