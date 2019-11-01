#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright © 2019 Martin Ueding <dev@martin-ueding.de>

import argparse
import json
import re
import sys


traces = {
    'C6cC': [[0, 1, 2, 3, 4, 5]],
    'C4cC': [[0, 1, 2, 3]],
}


def get_ordering(sequence):
    return [index for element, index in sorted(zip(sequence, range(len(sequence))))]


def process_name(name):
    parts = name.split('_')
    operators = parts[2:]

    ordering = get_ordering(operators[::2])
    rotate = 2 * ordering[0]

    rotated = operators[rotate:] + operators[:rotate]
    new_name = '_'.join(parts[:2] + rotated)

    return new_name


def main():
    options = _parse_args()

    with open(options.correlator_names) as f:
        all_names = json.load(f)

    for trace in traces:
        names = all_names[trace]
        new_names = [process_name(name) for name in names]

        print(len(names))
        print(len(set(new_names)))




def _parse_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('correlator_names')
    options = parser.parse_args()

    return options


if __name__ == '__main__':
    main()
