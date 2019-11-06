#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright © 2019 Martin Ueding <dev@martin-ueding.de>

import argparse
import json
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

    with open(options.correlator_names) as f:
        all_names = json.load(f)

    total_before = 0
    total_after = 0

    for trace, spec in traces.items():
        names = all_names[trace]
        new_names = [process_name(name, spec) for name in names]

        #print_difference_table(names, new_names)

        print('| {:6} | {:6} | {:6} |'.format(trace, len(names), len(set(new_names))))

        total_before += len(names)
        total_after += len(set(new_names))

        all_names[trace] = list(set(new_names))

    print()
    print('Reduction: {} to {}, saving of {:.1f} %.'.format(total_before, total_after, 100 * (1 - total_after/total_before)))

    with open(options.out, 'w') as f:
        json.dump(all_names, f, indent=2)


def _parse_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('correlator_names')
    parser.add_argument('out')
    options = parser.parse_args()

    return options


if __name__ == '__main__':
    main()
