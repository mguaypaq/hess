#!/usr/bin/env python
# -*- coding: utf-8 -*-

r"""
Code for computing chromatic symmetric functions.
"""

import itertools as it
import time
from collections import defaultdict, Counter, namedtuple

from path import *
from perm import *
from util import *

# ---------------------------------------------------------

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# ---------------------------------------------------------

def compute_csfs(n):
    r"""
    Compute the coefficients of the q-csf for everything of size n,
    in the monomial basis.
    """
    logger.info('starting size %d', n)
    paths = list(iter_path(n))
    perms = list(iter_blist(n))
    parts = list(partitions(n))
    for k, path in enumerate(paths):
        logger.info('doing size %d path %d', n, k)
        csf = {
            part: [0]*(n*(n-1)//2+1)
            for part in parts
            }
        for perm in perms:
            degree = inversions(path, perm)
            for part in parts:
                if contractible(path, perm, part):
                    csf[part][degree] += 1
        yield path, csf
    logger.info('done with size %d', n)

# ---------------------------------------------------------

def contractible(path, perm, composition):
    r"""
    Check whether contracting a permutation colouring according to
    a composition gives a valid chain colouring for a Dyck path.

    `path` describes a unit interval order

    `perm` gives a colouring of the poset

    `composition` gives chain lengths for the contracted colouring

    >>> any(contractible((0, 0, 0), perm, (2, 1))
    ...     for perm in it.permutations(range(3)))
    False
    >>> contractible((0, 0, 1), (0, 1, 2), (2, 1))
    False
    >>> contractible((0, 0, 1), (0, 2, 1), (2, 1))
    True
    >>> contractible((0, 0, 1), (1, 0, 2), (2, 1))
    False
    >>> contractible((0, 0, 1), (1, 2, 0), (2, 1))
    False
    >>> contractible((0, 0, 1), (2, 0, 1), (2, 1))
    False
    >>> contractible((0, 0, 1), (2, 1, 0), (2, 1))
    False
    >>> contractible((0, 1, 2), (0, 1, 2), (3,))
    True
    >>> contractible((0, 1, 2), (0, 2, 1), (3,))
    False
    >>> contractible((0, 1, 2), (1, 0, 2), (3,))
    False
    >>> contractible((0, 1, 2), (1, 2, 0), (3,))
    False
    >>> contractible((0, 1, 2), (2, 0, 1), (3,))
    False
    >>> contractible((0, 1, 2), (2, 1, 0), (3,))
    False
    """
    n = len(path)
#    assert is_path(path)
#    assert is_blist(perm) and len(perm) == n
    assert isinstance(composition, tuple) and sum(composition) == n
    comparable = boxes_under_path(path)
    total = 0
    for part in composition:
        for col1 in range(total, total+part-1):
            col2 = col1 + 1
            pos1 = perm.index(col1)
            pos2 = perm.index(col2)
            if pos1 > pos2 or (pos1, pos2) in comparable:
                return False
        total += part
    return True

# ---------------------------------------------------------

def doctest():
    import doctest
    doctest.testmod(verbose=False)

# ---------------------------------------------------------

def setup_logging():
    handler = logging.StreamHandler()
    handler.setFormatter(
        logging.Formatter(
            '%(module)s (elapsed time %(relativeCreated)d): %(message)s'))
    logger.addHandler(handler)

# ---------------------------------------------------------

def argparse():
    import argparse
    parser = argparse.ArgumentParser(
        description='Compute graded chromatic symmetric functions for all unit interval orders of size n.',
        )
    parser.add_argument(
        'n',
        type=int,
        help='The size of unit interval orders to consider.',
        )
    args = parser.parse_args()
    return args.n

# ---------------------------------------------------------

output_header = r"""csf[{path}] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
"""

output_footer = r"""    ])

"""

def save(path, csf):
    filename = 'output/csf-' + ''.join(map(str, path)) + '.py'
    with open(filename, 'w') as f:
        f.write(output_header.format(path=path))
        for index, coeffs in sorted(csf.iteritems()):
            while coeffs and coeffs[-1] == 0:
                coeffs.pop()
            if coeffs:
                f.write("    ({}, {}),\n".format(list(index), coeffs))
        f.write(output_footer)

# ---------------------------------------------------------

if __name__ == '__main__':
    doctest()
    setup_logging()
    n = argparse()
    for path, csf in compute_csfs(n):
        save(path, csf)

# ---------------------------------------------------------

