#!/usr/bin/env python
# -*- coding: utf-8 -*-

r"""
Code for computing characters of Hessenberg varieties.
"""

import itertools as it
import numpy as np
import operator
import time
from collections import defaultdict, Counter, namedtuple

from path import *
from perm import *
from util import *

# ---------------------------------------------------------

def translators(n):
    r"""
    Return the lperms for which we can compute character values.

    >>> translators(3) == {
    ...     (0, 1, 2): (1, 1, 1),
    ...     (1, 0, 2): (2, 1),
    ...     (2, 0, 1): (3,),
    ...     }
    True
    """
    result = {}
    for cc in compositions(n):
        if cc != tuple(sorted(cc, reverse=True)):
            continue
        cycle_type = tuple(sorted(cc, reverse=True))
        bfact = []
        for part in cc:
            bfact.extend([0]*(part-1))
            bfact.append(part-1)
        lperm = blist_from_bfact(tuple(bfact))
        result[lperm] = cycle_type
    return result

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

def rpspec(bfact, path):
    r"""
    Do some voodoo on a bijection and a path to get a root product spec.
    """
    mobile, fixed = [], []
    blist = []
    n = len(bfact)
    assert is_bfact(bfact)
    assert is_path(path)
    assert n == len(path)
    for i in range(n):
        j = len(blist) - bfact[i]
        blist.insert(j, i)
        crossed = False
        for jj in range(j+1, len(blist)):
            if blist[jj] < path[i]:
                crossed = True
            elif blist[jj] < i:
                if crossed:
                    fixed.append((blist[jj], i, 'f'))
                else:
                    mobile.append((jj, i, 'm'))
    return mobile, fixed

def rp(bfact, spec):
    r"""
    Do some more voodoo to transform a spec into an actual root product.
    """
    mobile, fixed = spec
    result = []
    for (j, i, f) in fixed:
        assert f == 'f'
        result.append((j, i))
    blist = []
    for i in range(len(bfact)):
        j = len(blist) - bfact[i]
        blist.insert(j, i)
        for (jj, ii, m) in mobile:
            assert m == 'm'
            if ii == i:
                result.append((blist[jj], i))
    return result

_flowup_cache = {}
def flowup(bfact, path):
    try:
        return _flowup_cache[bfact, path]
    except KeyError:
        result = _flowup_compute(bfact, path)
        _flowup_cache[bfact, path] = result
        return result
def _flowup_compute(bfact, path):
    r"""
    Return a fragment of a flowup basis vector.
    """
    spec = rpspec(bfact, path)
    result = {}
    projections = [
        ((i, i+1) if i < k else (i,))
        for k, i in enumerate(bfact)
        ]
    coords = list(it.product(*projections))
    assert all(is_bfact(c) for c in coords)
    for c in coords:
        result[blist_from_bfact(c)] = rp(c, spec)
    return result

# ---------------------------------------------------------

_indices_above_cache = {}
def indices_above(bfact):
    try:
        return _indices_above_cache[bfact]
    except KeyError:
        result = _indices_above_compute(bfact)
        _indices_above_cache[bfact] = result
        return result
def _indices_above_compute(bfact):
    result = {}
    assert is_bfact(bfact)
    n = len(bfact)
    for offset in it.product(*([(0,)] + [(0, 1)]*(n-1))):
        ofact = tuple(b+o for b, o in zip(bfact, offset))
        if is_bfact(ofact):
            olist = blist_from_bfact(ofact)
            result[offset] = olist
    return result

_indices_below_cache = {}
def indices_below(bfact):
    try:
        return _indices_below_cache[bfact]
    except KeyError:
        result = _indices_below_compute(bfact)
        _indices_below_cache[bfact] = result
        return result
def _indices_below_compute(bfact):
    result = {}
    assert is_bfact(bfact)
    n = len(bfact)
    maxoff = (0,) + (1,)*(n-1)
    bfact = tuple(b-m for b, m in zip(bfact, maxoff))
    for offset in it.product(*([(0,)] + [(0, 1)]*(n-1))):
        ofact = tuple(b+o for b, o in zip(bfact, offset))
        if is_bfact(ofact):
            olist = blist_from_bfact(ofact)
            result[offset] = olist
    return result

def frag_at(n, frag, indices):
    result = np.zeros(
        (1,) + (2,)*(n-1),
        dtype=object,
        )
    for offset, olist in indices.items():
        if olist in frag:
            result[offset] = frag[olist]
    return result

# ---------------------------------------------------------

def doleft(path):
    assert is_path(path)
    n = len(path)
    maxoff = (0,) + (1,)*(n-1)
    basis = {}
    for bfact in iter_bfact(n):
        basis[bfact] = frag_at(
            n,
            lvaluated_fragment(flowup(bfact, path)),
            indices_above(bfact),
            )
    csf = defaultdict(int)
    for t in translators(n):
        for bfact in iter_bfact(n):
            f = flowup(bfact, path)
            deg = len(f[blist_from_bfact(bfact)])
            work_array = frag_at(
                n,
                lvaluated_fragment(translated_fragment(t, f)),
                indices_below(bfact),
                )
            for offset in it.product(*([(0,)] + [(0, 1)]*(n-1))):
                coeff = work_array[offset]
                if coeff == 0:
                    quo = 0
                else:
                    ofact = tuple(b+o-m for b, o, m in zip(bfact, offset, maxoff))
                    ovect = basis[ofact]
                    olead = ovect[(0,)*n]
                    quo, rem = divmod(coeff, olead)
                    assert rem == 0
                    wa_indices = [
                        slice(None, None, None) if i == 0 else 1
                        for i in offset
                        ]
                    ov_indices = [
                        slice(None, None, None) if i == 0 else 0
                        for i in offset
                        ]
                    work_array[wa_indices] -= quo * ovect[ov_indices]
            csf[t,deg] += quo
    return csf

def doright(path):
    assert is_path(path)
    n = len(path)
    maxoff = (0,) + (1,)*(n-1)
    basis = {}
    for bfact in iter_bfact(n):
        basis[bfact] = frag_at(
            n,
            rvaluated_fragment(flowup(bfact, path)),
            indices_above(bfact),
            )
    csf = defaultdict(int)
    for t in translators(n):
        for bfact in iter_bfact(n):
            f = flowup(bfact, path)
            deg = len(f[blist_from_bfact(bfact)])
            work_array = frag_at(
                n,
                rvaluated_fragment(translated_fragment(t, f)),
                indices_below(bfact),
                )
            for offset in it.product(*([(0,)] + [(0, 1)]*(n-1))):
                coeff = work_array[offset]
                if coeff == 0:
                    quo = 0
                else:
                    ofact = tuple(b+o-m for b, o, m in zip(bfact, offset, maxoff))
                    ovect = basis[ofact]
                    olead = ovect[(0,)*n]
                    quo, rem = divmod(coeff, olead)
                    assert rem == 0
                    wa_indices = [
                        slice(None, None, None) if i == 0 else 1
                        for i in offset
                        ]
                    ov_indices = [
                        slice(None, None, None) if i == 0 else 0
                        for i in offset
                        ]
                    work_array[wa_indices] -= quo * ovect[ov_indices]
            csf[t,deg] += quo
    return csf

def check_rreg(path):
    r"""
    >>> all(check_rreg(path) == 24 for path in iter_path(4))
    True
    """
    assert is_path(path)
    sums = defaultdict(int)
    for ((t, _), coeff) in doright(path).iteritems():
        sums[t] += coeff
    result = sums.pop(tuple(range(len(path))))
    if all(s == 0 for s in sums.itervalues()):
        return result
    else:
        return False

# ---------------------------------------------------------

def output_all(sizes):
    start = time.time()
    result = ""
    for n in sizes:
        cycle_type = translators(n)
        k = 0
        for path in iter_path(n):
            k += 1
            print 'size', n, 'path', k, 'time', time.time()-start
            #right
            tmp = defaultdict(lambda: [0]*(1+n*(n-1)//2))
            for ((lperm, deg), coeff) in doright(path).iteritems():
                tmp[cycle_type[lperm]][deg] = coeff
            result += '''hess_right[{}] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
'''.format(path)
            for index, coeffs in sorted(tmp.iteritems()):
                while coeffs and coeffs[-1] == 0:
                    coeffs.pop()
                if coeffs:
                    result += '    ({}, {}),\n'.format(list(index), coeffs)
            result += '    ])\n\n'
            #left
            tmp = defaultdict(lambda: [0]*(1+n*(n-1)//2))
            for ((lperm, deg), coeff) in doleft(path).iteritems():
                tmp[cycle_type[lperm]][deg] = coeff
            result += '''hess_left[{}] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
'''.format(path)
            for index, coeffs in sorted(tmp.iteritems()):
                while coeffs and coeffs[-1] == 0:
                    coeffs.pop()
                if coeffs:
                    result += '    ({}, {}),\n'.format(list(index), coeffs)
            result += '    ])\n\n'
    return result

# ---------------------------------------------------------

def short_computation():
    with open('hess-12345.py', 'w') as f:
        f.write(output_all([1, 2, 3, 4, 5]))

def long_computation():
    with open('hess-123456.py', 'w') as f:
        f.write(output_all([1, 2, 3, 4, 5, 6]))
    with open('hess-7.py', 'w') as f:
        f.write(output_all([7]))
    with open('hess-8.py', 'w') as f:
        f.write(output_all([8]))
    print 'all done! :D'

# ---------------------------------------------------------

#check translation classes?
#test vector supports?

# ---------------------------------------------------------
import doctest
doctest.testmod()
# ---------------------------------------------------------

if __name__ == '__main__':
    short_computation()
#    long_computation()

