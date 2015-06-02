#!/usr/bin/env python
# -*- coding: utf-8 -*-

r"""
Code for computing chromatic symmetric functions.
"""

import itertools as it
import math
import random
import numpy as np
import operator
import time
from collections import defaultdict, Counter, namedtuple

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

def inversions(path, perm):
    r"""
    Compute the number of path-inversions for a given permutation.

    >>> inversions((0, 0, 0), (0, 1, 2))
    0
    >>> inversions((0, 0, 0), (0, 2, 1))
    1
    >>> inversions((0, 0, 0), (1, 0, 2))
    1
    >>> inversions((0, 0, 0), (1, 2, 0))
    2
    >>> inversions((0, 0, 0), (2, 0, 1))
    2
    >>> inversions((0, 0, 0), (2, 1, 0))
    3
    >>> inversions((0, 0, 1), (0, 1, 2))
    0
    >>> inversions((0, 0, 1), (0, 2, 1))
    1
    >>> inversions((0, 0, 1), (1, 0, 2))
    1
    >>> inversions((0, 0, 1), (1, 2, 0))
    1
    >>> inversions((0, 0, 1), (2, 0, 1))
    1
    >>> inversions((0, 0, 1), (2, 1, 0))
    2
    """
    return sum((perm[i] > perm[j]) for i, j in boxes_under_path(path))

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

def is_blist(bl):
    r"""
    Check whether `bl` is a valid bijection, in the index list format.

    In the index list format, `bl[i] == j` if ``L_i \mapsto R_j``.

    >>> is_blist((0, 0, 0, 0, 0))
    False
    >>> is_blist((0, 1, 2, 3, 4))
    True
    >>> is_blist((4, 3, 2, 1, 0))
    True
    """
    return (isinstance(bl, tuple) and
            sorted(bl) == range(len(bl)))

def iter_blist(n):
    r"""
    Return an iterator over all valid blists of length `n`.
    """
    return it.permutations(range(n))

def test_is_blist(below=7):
    r"""
    Test that `is_blist` corresponds to `iter_blist`.

    >>> test_is_blist()
    """
    for n in range(below):
        expected = list(iter_blist(n))
        actual = filter(is_blist, it.product(range(n), repeat=n))
        assert expected == actual

# ---------------------------------------------------------

def is_path(p):
    r"""
    Check whether `p` is a valid Dyck path.

    The entries of `p` give the number of boxes above the path in each row.

    >>> is_path((0, 0, 0))
    True
    >>> is_path((0, 0, 1))
    True
    >>> is_path((0, 0, 2))
    True
    >>> is_path((0, 1, 1))
    True
    >>> is_path((0, 1, 2))
    True
    >>> is_path((0, 2, 1))
    False
    """
    return (isinstance(p, tuple) and
            sorted(p) == list(p) and
            all(p[i] in range(i+1)
                for i in range(len(p))))

def iter_path(n):
    r"""
    Return an iterator over all valid Dyck paths of length `n`.
    """
    if n == 0:
        yield ()
        return
    elif n == 1:
        yield (0,)
        return
    elif n >= 2:
        for head in iter_path(n-1):
            for tail in range(head[-1], n):
                yield head + (tail,)

def test_is_path(below=7):
    r"""
    Test that `is_path` corresponds to `iter_path`.

    >>> test_is_path()
    """
    for n in range(below):
        expected = list(iter_path(n))
        actual = filter(is_path, it.product(range(n), repeat=n))
        assert expected == actual

_boxes_cache = {}
def boxes_under_path(p):
    try:
        return _boxes_cache[p]
    except KeyError:
        result = _boxes_compute(p)
        _boxes_cache[p] = result
        return result
def _boxes_compute(p):
    r"""
    Return the set of boxes `(i, j)` below the path `p`.

    >>> sorted(boxes_under_path((0, 0, 0)))
    [(0, 1), (0, 2), (1, 2)]
    >>> sorted(boxes_under_path((0, 0, 1)))
    [(0, 1), (1, 2)]
    >>> sorted(boxes_under_path((0, 0, 2)))
    [(0, 1)]
    >>> sorted(boxes_under_path((0, 1, 1)))
    [(1, 2)]
    >>> sorted(boxes_under_path((0, 1, 2)))
    []
    """
    return {(i, j) for j, k in enumerate(p) for i in range(k, j)}

# ---------------------------------------------------------

def compositions(n):
    r"""
    Iterator for the compositions of the integer `n`.

    >>> sorted(list(compositions(3)))
    [(1, 1, 1), (1, 2), (2, 1), (3,)]
    >>> len(list(compositions(8)))
    128
    """
    if n == 0:
        yield ()
        return
    for head in range(1, n+1):
        for tail in compositions(n-head):
            yield (head,) + tail

# ---------------------------------------------------------

def partitions(n):
    r"""
    Iterator for the partitions of the integer `n`.

    >>> sorted(list(partitions(3)))
    [(1, 1, 1), (2, 1), (3,)]
    >>> len(list(partitions(8)))
    22
    """
    for c in compositions(n):
        if c == tuple(sorted(c, reverse=True)):
            yield c

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

