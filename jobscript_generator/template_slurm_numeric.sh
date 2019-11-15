#!/bin/bash

#SBATCH --job-name N_{{ '%04d'|format(config_number) }}
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=ueding@hiskp.uni-bonn.de
#SBATCH --mail-type=FAIL

#SBATCH --output=batch_output/numeric_{{ config_number }}_slurm_%j.txt

set -e
set -u
set -x

hostname
date -Iseconds

config="{{ '%04d'|format(config_number) }}"

tempdir="/storage/ueding/correlators/${config}"
mkdir -p "$tempdir"
cp correlators/*_cnfg${config}.h5 "$tempdir"

cleanup() {
    rm -rf "$tempdir"
}

trap cleanup EXIT


{% for momentum_sq, grouped in grouped2 %}
###############################################################################
#                                    P² = {{ momentum_sq }}                                   #
###############################################################################
{% for irrep, values in grouped %}
# {{ irrep }}
{% for _, momentum in values -%}
/usr/bin/time {{ srcdir }}/numeric_projection/driver.R {{ momentum|join(' ') }} {{ irrep }} {{ config_number }} "$tempdir"
{% endfor -%}
{% endfor -%}
{% endfor -%}

date -Iseconds
