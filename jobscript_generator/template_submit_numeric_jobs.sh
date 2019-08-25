#!/bin/bash

set -e
set -u

{% for momentum_sq, grouped in grouped2 %}
###############################################################################
#                                    P² = {{ momentum_sq }}                                   #
###############################################################################
{% for irrep, values in grouped %}
# {{ irrep }}
{% for _, momentum in values -%}
{% for config_number in config_numbers -%}
sbatch jobscripts/slurm_numeric_{{ momentum }}_{{ irrep }}_{{ '%04d'|format(config_number) }}.sh
{% endfor -%}
{% endfor -%}
{% endfor -%}
{% endfor -%}
