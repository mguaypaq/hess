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
log = logging.getLogger(__name__)

log.setLevel(logging.INFO)
h = logging.StreamHandler()
h.setFormatter(logging.Formatter('%(relativeCreated)d: %(message)s'))
log.addHandler(h)
del h

log.info('starting script')

# ---------------------------------------------------------

def doit(n):
    r"""
    Compute the coefficients of the q-csf for everything of size n,
    in the monomial basis.
    """
    log.info('starting size %d', n)
    result = {}
    paths = list(iter_path(n))
    perms = list(iter_blist(n))
    parts = list(partitions(n))
    for k, path in enumerate(paths):
        log.info('doing size %d path %d', n, k)
        current = {
            part: [0]*(n*(n-1)//2+1)
            for part in parts
            }
        for perm in perms:
            degree = inversions(path, perm)
            for part in parts:
                if contractible(path, perm, part):
                    current[part][degree] += 1
        result[path] = current
    log.info('done with size %d', n)
    return result

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

def output_all(sizes):
    start = time.time()
    result = ""
    for n in sizes:
        data = doit(n)
        for path in iter_path(n):
            result += '''csf[{}] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
'''.format(path)
            for index, coeffs in sorted(data[path].iteritems()):
                while coeffs and coeffs[-1] == 0:
                    coeffs.pop()
                if coeffs:
                    result += '    ({}, {}),\n'.format(list(index), coeffs)
            result += '    ])\n\n'
    return result

# ---------------------------------------------------------

def short_computation():
    with open('csf-123456.py', 'w') as f:
        f.write(output_all([1, 2, 3, 4, 5, 6]))

def long_computation():
    for n in [1, 2, 3, 4, 5, 6, 7, 8]:
        with open('csf-{}.py'.format(n), 'w') as f:
            f.write(output_all([n]))

# ---------------------------------------------------------
import doctest
doctest.testmod()
# ---------------------------------------------------------

if __name__ == '__main__':
    long_computation()

