#!/usr/bin/env python
# -*- coding: utf-8 -*-

r"""
Code for working with Dyck paths.
"""

__all__ = [
    'is_fragment',
    'lvaluated_fragment',
    'rvaluated_fragment',
    'translated_fragment',
    ]

# ---------------------------------------------------------

import itertools as it
import operator

from path import *
from perm import *

# ---------------------------------------------------------

def normalize_root_product(roots):
    r"""
    Return a `sign, roots` pair where `sign` is ``\pm 1`` or 0
    and `roots` is a sorted tuple of roots `(i, j)` with `i < j`.

    >>> normalize_root_product([])
    (1, ())
    >>> normalize_root_product([(0, 1), (0, 2), (1, 2)])
    (1, ((0, 1), (0, 2), (1, 2)))
    >>> normalize_root_product([(0, 1), (0, 2), (2, 1)])
    (-1, ((0, 1), (0, 2), (1, 2)))
    >>> normalize_root_product([(0, 1), (2, 0), (2, 1)])
    (1, ((0, 1), (0, 2), (1, 2)))
    >>> normalize_root_product([(0, 1), (1, 2), (0, 2)])
    (1, ((0, 1), (0, 2), (1, 2)))
    >>> normalize_root_product([(0, 1), (0, 0), (1, 2)])
    (0, ())
    """
    sign = 1
    result = []
    for (i, j) in roots:
        if i < j:
            result.append((i, j))
        elif i > j:
            result.append((j, i))
            sign = -sign
        else:
            return 0, ()
    result.sort()
    return sign, tuple(result)

def project_root_product(roots, fro, to):
    r"""
    Return a `sign, roots` pair `sign` is ``\pm 1`` or 0
    and `roots` is obtained by replacing every instance
    of `fro` by `to` in the given roots (and normalizing).

    >>> project_root_product([(0, 1)], fro=1, to=3)
    (1, ((0, 3),))
    >>> project_root_product([(0, 1)], fro=2, to=3)
    (1, ((0, 1),))
    >>> project_root_product([(0, 1)], fro=3, to=1)
    (1, ((0, 1),))
    >>> project_root_product([(0, 1)], fro=0, to=3)
    (-1, ((1, 3),))
    >>> project_root_product([(0, 1)], fro=0, to=1)
    (0, ())
    """
    return normalize_root_product(
        (i if i != fro else to,
         j if j != fro else to)
        for (i, j) in roots)

# ---------------------------------------------------------

def is_fragment(frag, path):
    r"""
    Check whether `frag` is a valid vector fragment.

    A vector fragment is given as a dict of (blist->root_product).
    The coefficient for a missing blist is assumed to be zero or don't-care.
    The fragment is valid if the coefficients for the given
    blists satisfy the divisibility conditions associated with `path`
    with all blists that are *below* them.

    >>> test_is_fragment()
    """
    boxes = boxes_under_path(path)
    for bl in frag:
        if not is_blist(bl): return False
        for i, j in it.combinations(range(len(bl)), 2):
            if bl[i] > bl[j] and (bl[j], bl[i]) in boxes:
                tmp = list(bl)
                tmp[i], tmp[j] = tmp[j], tmp[i]
                tmp = tuple(tmp)
                root_prod_above = frag[bl]
                root_prod_below = frag.get(tmp, ((0, 0),))
                if (project_root_product(root_prod_above, fro=bl[i], to=bl[j]) !=
                    project_root_product(root_prod_below, fro=bl[i], to=bl[j])):
                        return False
    return True

def prod_lperm(lp1, lp2):
    r"""
    Return the composition of the permutations of L variables `lp1` and `lp2`.

    >>> lp1 = (1, 0, 2, 3, 4)
    >>> lp2 = (1, 2, 3, 4, 0)
    >>> prod_lperm(lp1, lp2)
    (2, 1, 3, 4, 0)
    >>> prod_lperm(lp2, lp1)
    (0, 2, 3, 4, 1)
    """
    n = len(lp1)
    assert n == len(lp2)
    assert is_blist(lp1)
    assert is_blist(lp2)
    return tuple(lp2[lp1[i]] for i in range(n))

def translated_fragment(lperm, frag):
    r"""
    Return the result of acting by the permutation of L variables `lperm`
    on the vector fragment `frag`.

    Note that the result may not be valid according to `is_fragment`,
    but this may be fine anyway.

    >>> lp1 = (1, 0, 2)
    >>> lp2 = (1, 2, 0)
    >>> frag = {
    ...     (0, 1, 2): (),
    ...     (1, 0, 2): ((0, 1),),
    ...     (0, 2, 1): ((0, 1), (0, 1)),
    ...     (1, 2, 0): ((0, 1), (0, 1), (0, 1)),
    ...     }
    >>> translated_fragment(lp1, frag) == {
    ...     (1, 0, 2): (),
    ...     (0, 1, 2): ((0, 1),),
    ...     (2, 0, 1): ((0, 1), (0, 1)),
    ...     (2, 1, 0): ((0, 1), (0, 1), (0, 1)),
    ...     }
    True
    >>> translated_fragment(lp2, frag) == {
    ...     (1, 2, 0): (),
    ...     (0, 2, 1): ((0, 1),),
    ...     (2, 1, 0): ((0, 1), (0, 1)),
    ...     (2, 0, 1): ((0, 1), (0, 1), (0, 1)),
    ...     }
    True
    """
    result = {}
    for bl in frag:
        tmp = tuple(bl[lperm[i]] for i in range(len(lperm)))
        result[tmp] = frag[bl]
    return result

def lvaluated_fragment(frag):
    r"""
    Return the result of evaluating all coordinates of `frag` at ``L_i = i``.

    >>> frag = {
    ...     (0, 1, 2): ((0, 1), (0, 2)),
    ...     (1, 0, 2): ((0, 1), (0, 2)),
    ...     (0, 2, 1): ((0, 1), (0, 2)),
    ...     (1, 2, 0): ((0, 1), (0, 2)),
    ...     (2, 0, 1): ((0, 1), (0, 2)),
    ...     (2, 1, 0): ((0, 1), (0, 2)),
    ...     }
    >>> lvaluated_fragment(frag) == {
    ...     (0, 1, 2): 2,
    ...     (1, 0, 2): -1,
    ...     (0, 2, 1): 2,
    ...     (1, 2, 0): 2,
    ...     (2, 0, 1): -1,
    ...     (2, 1, 0): 2,
    ...     }
    True
    """
    result = {}
    for bl, coeff in frag.iteritems():
        result[bl] = reduce(
            operator.mul,
            [bl.index(i)-bl.index(j) for (i, j) in coeff],
            1,
            )
    return result

def rvaluated_fragment(frag):
    r"""
    Return the result of evaluating all coordinates of `frag` at ``R_i = i``.

    >>> frag = {
    ...     (0, 1, 2): (),
    ...     (1, 0, 2): ((0, 1), (0, 2)),
    ...     (0, 2, 1): ((0, 2), (0, 2), (0, 2)),
    ...     (1, 2, 0): ((0, 1), (1, 2)),
    ...     (2, 0, 1): ((0, 2),),
    ...     (2, 1, 0): ((0, 1), (0, 2), (0, 2), (0, 2), (0, 2)),
    ...     }
    >>> rvaluated_fragment(frag) == {
    ...     (0, 1, 2): 1,
    ...     (1, 0, 2): 2,
    ...     (0, 2, 1): -8,
    ...     (1, 2, 0): 1,
    ...     (2, 0, 1): -2,
    ...     (2, 1, 0): -16,
    ...     }
    True
    """
    result = {}
    for bl, coeff in frag.iteritems():
        result[bl] = reduce(
            operator.mul,
            [i-j for (i, j) in coeff],
            1,
            )
    return result

# ---------------------------------------------------------

def test_is_fragment():
    r"""
    Test the function `is_fragment` on various cases.

    >>> test_is_fragment()
    """
    paths = [
        (0, 0, 0),
        (0, 0, 1),
        (0, 0, 2),
        (0, 1, 1),
        (0, 1, 2),
        ]
    frags = [
        {#0
            (0, 1, 2): (),
            (1, 0, 2): (),
            (0, 2, 1): (),
            (1, 2, 0): (),
            (2, 0, 1): (),
            (2, 1, 0): (),
        },
        {#1
            (1, 0, 2): ((0, 1),),
            (1, 2, 0): ((0, 1),),
            (2, 0, 1): ((0, 2),),
            (2, 1, 0): ((0, 2),),
        },
        {#2
            (0, 2, 1): ((1, 2),),
            (1, 2, 0): ((0, 2),),
            (2, 0, 1): ((1, 2),),
            (2, 1, 0): ((0, 2),),
        },
        {#3
            (1, 2, 0): ((0, 1), (0, 2)),
            (2, 1, 0): ((0, 1), (0, 2)),
        },
        {#4
            (2, 0, 1): ((0, 2), (1, 2)),
            (2, 1, 0): ((0, 2), (1, 2)),
        },
        {#5
            (2, 1, 0): ((0, 1), (0, 2), (1, 2)),
        },
        {#6
            (1, 2, 0): ((0, 1),),
            (2, 1, 0): ((0, 1),),
        },
        {#7
            (2, 0, 1): ((1, 2),),
            (2, 1, 0): ((1, 2),),
        },
        {#8
            (2, 1, 0): ((0, 1), (1, 2)),
        },
        {#9
            (1, 0, 2): (),
        },
        {#10
            (1, 2, 0): (),
        },
        {#11
            (2, 1, 0): (),
        },
        {#12
            (1, 2, 0): ((0, 1),),
            (2, 1, 0): ((0, 0),),
        },
        ]
    assert ([is_fragment(frag, path) for frag in frags for path in paths] ==
        [
        1, 1, 1, 1, 1, #0
        1, 1, 1, 1, 1, #1
        1, 1, 1, 1, 1, #2
        1, 1, 1, 1, 1, #3
        1, 1, 1, 1, 1, #4
        1, 1, 1, 1, 1, #5
        0, 1, 1, 1, 1, #6
        0, 1, 1, 1, 1, #7
        0, 1, 1, 1, 1, #8
        0, 0, 0, 1, 1, #9
        0, 0, 0, 1, 1, #10
        0, 0, 0, 0, 1, #11
        0, 0, 1, 0, 1, #12
        ])

def test_associativity():
    r"""
    Test that `prod_lperm` and `translated_fragment` are associative
    (as a group and a group action).

    >>> test_associativity()
    """
    perms = [
        (0, 1, 2),
        (1, 0, 2),
        (0, 2, 1),
        (1, 2, 0),
        (2, 0, 1),
        (2, 1, 0),
        ]
    frag = {
        (0, 1, 2): (),
        (1, 0, 2): ((0, 1),),
        (0, 2, 1): ((0, 1), (0, 1)),
        (1, 2, 0): ((0, 1), (0, 1), (0, 1)),
        }
    assert all(
        prod_lperm(x, prod_lperm(y, z)) == prod_lperm(prod_lperm(x, y), z)
        for x in perms
        for y in perms
        for z in perms
        )
    assert all(
        translated_fragment(x, translated_fragment(y, frag)) ==
        translated_fragment(prod_lperm(x, y), frag)
        for x in perms
        for y in perms
        )

# ---------------------------------------------------------
if __name__ == '__main__':
    import doctest
    doctest.testmod()
# ---------------------------------------------------------

