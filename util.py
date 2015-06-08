#!/usr/bin/env python
# -*- coding: utf-8 -*-

r"""
Various combinatorial utilities that don't fit elsewhere.
"""

__all__ = [
    'compositions',
    'inversions',
    'partitions',
    ]

# ---------------------------------------------------------

import itertools as it
import numpy as np
import operator
from collections import defaultdict, Counter, namedtuple

from path import *
from perm import *

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

def partitions(n, bound=None):
    r"""
    Iterator for the partitions of the integer `n`.

    Optional argument `bound` is an upper bound on the part sizes.

    >>> sorted(list(partitions(3)))
    [(1, 1, 1), (2, 1), (3,)]
    >>> len(list(partitions(8)))
    22
    >>> sorted(list(partitions(4, 2)))
    [(1, 1, 1, 1), (2, 1, 1), (2, 2)]
    """
    if n == 0:
        yield ()
        return
    if bound is None or bound > n:
        bound = n
    for head in range(1, bound+1):
        for tail in partitions(n-head, head):
            yield (head,) + tail

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
if __name__ == '__main__':
    import doctest
    doctest.testmod()
# ---------------------------------------------------------

