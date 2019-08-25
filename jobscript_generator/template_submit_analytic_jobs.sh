#!/bin/bash

set -e
set -u

{% for momentum_sq, grouped in grouped2 %}
# P² = {{ momentum_sq }}
{%- for irrep, values in grouped %}
{% for _, momentum in values -%}
sbatch jobscripts/slurm_analytic_{{ momentum }}_{{ irrep }}.sh
{% endfor -%}
{% endfor -%}
{% endfor -%}
