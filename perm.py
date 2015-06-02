#!/usr/bin/env python
# -*- coding: utf-8 -*-

r"""
Code for working with permutations.

There are two useful representations:

 - blist: tuple of images
 - bfact: the code of a permutation
"""

__all__ = [
    'bfact_from_blist',
    'blist_from_bfact',
    'is_bfact',
    'is_blist',
    'iter_bfact',
    'iter_blist',
    ]

# ---------------------------------------------------------
import itertools as it
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

def is_bfact(bf):
    r"""
    Check whether `bf` is a valid bijection, in the factorial base format.

    In the factorial base format, `bf[i]` is the number of variables
    ``R_j`` with `j < i` that appear after ``R_i``.

    >>> is_bfact((0, 0, 0, 0, 0))
    True
    >>> is_bfact((0, 1, 2, 3, 4))
    True
    >>> is_bfact((4, 3, 2, 1, 0))
    False
    """
    return (isinstance(bf, tuple) and
            all(bf[i] in range(i+1)
                for i in range(len(bf))))

def iter_bfact(n):
    r"""
    Return an iterator over all valid bfacts of length `n`.
    """
    return it.product(*[range(i+1) for i in range(n)])

def bfact_from_blist(bl):
    r"""
    Convert a blist into a bfact.

    >>> bfact_from_blist((0, 1, 2, 3, 4))
    (0, 0, 0, 0, 0)
    >>> bfact_from_blist((4, 3, 2, 1, 0))
    (0, 1, 2, 3, 4)
    >>> bfact_from_blist((1, 0, 2, 3, 4))
    (0, 1, 0, 0, 0)
    """
    result = [0] * len(bl)
    for i, j in it.combinations(bl, 2):
        if j < i:
            result[i] += 1
    return tuple(result)

def blist_from_bfact(bf):
    r"""
    Convert a bfact into a blist.

    >>> blist_from_bfact((0, 0, 0 ,0, 0))
    (0, 1, 2, 3, 4)
    >>> blist_from_bfact((0, 1, 2, 3, 4))
    (4, 3, 2, 1, 0)
    >>> blist_from_bfact((0, 1, 0, 0, 0))
    (1, 0, 2, 3, 4)
    """
    result = []
    for i, inversions in enumerate(bf):
        result.insert(len(result)-inversions, i)
    return tuple(result)

# ---------------------------------------------------------

def test_is_blist(below=7):
    r"""
    Test that `is_blist` corresponds to `iter_blist`.

    >>> test_is_blist()
    """
    for n in range(below):
        expected = list(iter_blist(n))
        actual = filter(is_blist, it.product(range(n), repeat=n))
        assert expected == actual

def test_is_bfact(below=7):
    r"""
    Test that `is_bfact` corresponds to `iter_bfact`.

    >>> test_is_bfact()
    """
    for n in range(below):
        expected = list(iter_bfact(n))
        actual = filter(is_bfact, it.product(range(n), repeat=n))
        assert expected == actual

def test_convert(below=7):
    r"""
    Test that converting blist -> bfact -> blist is the
    identity function on valid blists, and similarly for
    bfact -> blist -> bfact on valid bfacts.

    >>> test_convert()
    """
    for n in range(below):
        for bl in iter_blist(n):
            assert bl == blist_from_bfact(bfact_from_blist(bl))
        for bf in iter_bfact(n):
            assert bf == bfact_from_blist(blist_from_bfact(bf))

def test_truncate(below=7):
    r"""
    Test that truncating a bfact corresponds to removing
    the largest elements in the equivalent blist.

    >>> test_truncate()
    """
    for n in range(1, below):
        for bf in iter_bfact(n):
            bl1 = blist_from_bfact(bf[:-1])
            bl2 = blist_from_bfact(bf)
            bl2 = list(bl2)
            bl2.remove(n-1)
            bl2 = tuple(bl2)
            assert bl1 == bl2

# ---------------------------------------------------------
if __name__ == '__main__':
    import doctest
    doctest.testmod()
# ---------------------------------------------------------

