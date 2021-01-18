#!/bin/bash

set -e
set -u

{% for config_number in config_numbers -%}
sbatch -x lcpunode02 jobscripts/slurm_numeric_{{ '%04d'|format(config_number) }}.sh
{% endfor -%}
