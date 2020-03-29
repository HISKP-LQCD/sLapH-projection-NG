#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright © 2019 Martin Ueding <mu@martin-ueding.de>

import argparse
import itertools
import json


def make_momenta():
    base_momenta = [
        (0, 0, 0),

        (0, 0, 1),
        (0, 1, 0),
        (1, 0, 0),

        (0, 1, 1),
        (1, 0, 1),
        (1, 1, 0),

        (1, 1, 1),

        (0, 0, 2),
        (0, 2, 0),
        (2, 0, 0),
    ]

    momenta = []

    for px, py, pz in base_momenta:
        for a, b, c in itertools.product(*([[+1, -1]] * 3)):
            p = (a * px, b * py, c * pz)
            momenta.append(p)

    return momenta


def main():
    options = _parse_args()

    correlators = []

    for px, py, pz in make_momenta():
        correlators.append('C2c_uu_p{}{}{}.d000.g5_p{}{}{}.d000.g5'.format(
            px, py, pz,
            -px, -py, -pz,
        ))

    print(json.dumps(correlators, indent=2))

def _parse_args():
    parser = argparse.ArgumentParser(description='')
    options = parser.parse_args()

    return options


if __name__ == '__main__':
    main()
