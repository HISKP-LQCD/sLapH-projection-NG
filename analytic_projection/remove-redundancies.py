#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright © 2019 Martin Ueding <martin-ueding.de>

import argparse
import json
import os
import random
import re
import sys


traces = {
    'C6cD': [[0, 1], [2, 3], [4, 5]],
    'C4cD': [[0, 1], [2, 3]],
    'C6cC': [[0, 1, 2, 3, 4, 5]],
    'C4cC': [[0, 1, 2, 3]],
    'C6cCD': [[0, 1, 2, 3], [4, 5]],
}


def get_ordering(sequence):
    return [index for element, index in sorted(zip(sequence, range(len(sequence))))]


def rotate_trace(operators):
    ordering = get_ordering(operators[::2])
    rotate = 2 * ordering[0]
    rotated = operators[rotate:] + operators[:rotate]
    return rotated


def process_name(name, spec):
    parts = name.split('_')
    operators = parts[2:]
    new_operators = operators[:]

    # Rotate the elements within each trace.
    for trace in spec:
        ops = [operators[i] for i in trace]
        new_ops = rotate_trace(ops)

        for i, op in zip(trace, new_ops):
            new_operators[i] = op

    # For the direct (D) diagrams we can also exchange the traces because there
    # are multiple of them. We first need to identify the traces of length 2.
    tr2 = [i for i, tr in enumerate(spec) if len(tr) == 2]
    if len(tr2) > 1:
        chunks = [new_operators[2*i:2*(i+1)] for i in range(len(tr2))]
        old_chunks = chunks[:]
        chunks.sort()

        new_operators = [op for chunk in chunks for op in chunk]

    new_name = '_'.join(parts[:2] + new_operators)

    return new_name


def print_difference_table(names, new_names):
    differences = [(old, new) for old, new in zip(names, new_names) if old != new]

    print()
    random.seed(0)
    for old, new in random.sample(differences, 30):
        print('| `{}` | `{}` |'.format(old, new))
    print()


def main():
    options = _parse_args()

    if not os.path.isdir(options.out):
        os.makedirs(options.out)

    for path in options.prescription:
        print(path)
        with open(path) as f:
            prescription = json.load(f)

        for total_momentum, value in prescription.items():
            print(' ', total_momentum)
            for irrep, value in value.items():
                print('   ', irrep)
                for irrep_col, value in value.items():
                    for irrep_row, value in value.items():
                        for gevp_row, value in value.items():
                            for gevp_col, value in value.items():
                                for i, value in enumerate(value):
                                    kind = value['datasetname'].split('_')[0]
                                    prescription[total_momentum][irrep][irrep_col][irrep_row][gevp_row][gevp_col][i]['datasetname'] = value['datasetname'] = process_name(value['datasetname'], traces[kind])

        basename = os.path.basename(path)
        path_out = os.path.join(options.out, basename)

        with open(path_out, 'w') as f:
            json.dump(prescription, f, indent=2)


def _parse_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('prescription', nargs='+')
    parser.add_argument('--out', required=True)
    options = parser.parse_args()

    return options


if __name__ == '__main__':
    main()
