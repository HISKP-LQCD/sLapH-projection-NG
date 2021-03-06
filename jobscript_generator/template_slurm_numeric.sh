#!/bin/bash

#SBATCH --job-name N_{{ '%04d'|format(config_number) }}
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=5GB
#SBATCH --mail-user=fischer@hiskp.uni-bonn.de
#SBATCH --mail-type=FAIL

#SBATCH --output=batch_output/numeric_{{ config_number }}_slurm_%j.txt

#set -e
set -u
set -x

hostname
date -Iseconds

{% for momentum_sq, grouped in grouped2 %}
###############################################################################
#                                    P² = {{ momentum_sq }}                                   #
###############################################################################
{% for irrep, values in grouped %}
# {{ irrep }}
{% for _, momentum in values -%}
if ! [[ -f "projected/resolved_{{ momentum|join('') }}_{{ irrep }}_{{ '%04d'|format(config_number) }}.js" ]]; then
    /usr/bin/time Rscript -e 'numericprojection::numeric_projection(c({{ momentum|join(', ') }}), "{{ irrep }}", {{ config_number }})'
    echo "Status: ${?}"
fi
{% endfor -%}
{% endfor -%}
{% endfor -%}

date -Iseconds
