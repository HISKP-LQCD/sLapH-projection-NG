#!/bin/bash

#SBATCH --job-name N_{{ '%04d'|format(config_number) }}
#SBATCH --cpus-per-task=1
#SBATCH --mem=950MB
#SBATCH --gres=disk:1
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

cleanup() {
    rm -rf "$tempdir"
}

has_copied=false

copy-configs() {
    if [[ "$has_copied" = "false" ]]; then
        cp correlators/*_cnfg${config}.h5 "$tempdir"
    fi
    has_copied=true
}

trap cleanup EXIT


{% for momentum_sq, grouped in grouped2 %}
###############################################################################
#                                    P² = {{ momentum_sq }}                                   #
###############################################################################
{% for irrep, values in grouped %}
# {{ irrep }}
{% for _, momentum in values -%}
if ! [[ -f "projected/resolved_{{ momentum|join('') }}_{{ irrep }}_{{ '%04d'|format(config_number) }}.js" ]]; then
    copy-configs
    /usr/bin/time {{ srcdir }}/numeric_projection/driver.R {{ momentum|join(' ') }} {{ irrep }} {{ config_number }} "$tempdir"
fi
{% endfor -%}
{% endfor -%}
{% endfor -%}

date -Iseconds
