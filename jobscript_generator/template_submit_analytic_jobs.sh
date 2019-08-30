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
sbatch jobscripts/slurm_analytic_{{ momentum|join('') }}_{{ irrep }}.sh
{% endfor -%}
{% endfor -%}
{% endfor -%}
