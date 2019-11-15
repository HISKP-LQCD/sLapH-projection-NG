#!/bin/bash

#SBATCH --job-name N_{{ '%04d'|format(config_number) }}
#SBATCH --cpus-per-task=1
#SBATCH --mem=2000MB
#SBATCH --mail-user=ueding@hiskp.uni-bonn.de
#SBATCH --mail-type=FAIL

#SBATCH --output=batch_output/numeric_{{ config_number }}_slurm_%j.txt

set -e
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
/usr/bin/time {{ srcdir }}/numeric_projection/driver.R {{ momentum|join(' ') }} {{ irrep }} {{ config_number }}
{% endfor -%}
{% endfor -%}
{% endfor -%}

date -Iseconds
