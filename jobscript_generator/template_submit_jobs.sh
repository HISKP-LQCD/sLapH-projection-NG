#!/bin/bash

set -e
set -u

{% for momentum_sq, grouped in grouped2 %}
# P² = {{ momentum_sq }}
{%- for irrep, values in grouped %}
{% for _, momentum in values -%}
sbatch jobscripts/slurm_{{ action }}_{{ momentum }}_{{ irrep }}.sh
{% endfor -%}
{% endfor -%}
{% endfor -%}
