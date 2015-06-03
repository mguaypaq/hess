#!/usr/bin/env python
# -*- coding: utf-8 -*-

r"""
Code for working with Dyck paths.

There are two useful representations:

Dyck paths are represented as numbers of missing boxes in
each downward row. For example, the following Dyck paths:

  /\
 /\/\
/\/\/\
(0, 0, 0)

 /\/\
/\/\/\
(0, 0, 1)

 /\
/\/\/\
(0, 0, 2)

   /\
/\/\/\
(0, 1, 1)

/\/\/\
(0, 1, 2)


To run some tests for this module, use the command:
$ python path.py
For more verbose output, use:
$ python path.py -v

"""

__all__ = [
    'boxes_under_path',
    'is_path',
    'iter_path',
    ]

# ---------------------------------------------------------

import itertools as it

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

def test_is_path(below=7):
    r"""
    Test that `is_path` corresponds to `iter_path`.

    >>> test_is_path()
    """
    for n in range(below):
        expected = list(iter_path(n))
        actual = filter(is_path, it.product(range(n), repeat=n))
        assert expected == actual

# ---------------------------------------------------------
if __name__ == '__main__':
    import doctest
    doctest.testmod()
# ---------------------------------------------------------

