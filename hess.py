#!/usr/bin/env python
# -*- coding: utf-8 -*-

r"""
Code for computing characters of Hessenberg varieties.
"""

import itertools as it
import math
import random
import numpy as np
import operator
import time
from collections import defaultdict, Counter, namedtuple

# ---------------------------------------------------------

def is_blist(bl):
    """
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
    """
    Return an iterator over all valid blists of length `n`.
    """
    return it.permutations(range(n))

def is_bfact(bf):
    """
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
    """
    Return an iterator over all valid bfacts of length `n`.
    """
    return it.product(*[range(i+1) for i in range(n)])

def bfact_from_blist(bl):
    """
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
    """
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

def test_is_blist(below=7):
    """
    Test that `is_blist` corresponds to `iter_blist`.

    >>> test_is_blist()
    """
    for n in range(below):
        expected = list(iter_blist(n))
        actual = filter(is_blist, it.product(range(n), repeat=n))
        assert expected == actual

def test_is_bfact(below=7):
    """
    Test that `is_bfact` corresponds to `iter_bfact`.

    >>> test_is_bfact()
    """
    for n in range(below):
        expected = list(iter_bfact(n))
        actual = filter(is_bfact, it.product(range(n), repeat=n))
        assert expected == actual

def test_convert(below=7):
    """
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
    """
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

def normalize_root_product(roots):
    """
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
    """
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

def is_path(p):
    """
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
    """
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
    """
    Test that `is_path` corresponds to `iter_path`.

    >>> test_is_path()
    """
    for n in range(below):
        expected = list(iter_path(n))
        actual = filter(is_path, it.product(range(n), repeat=n))
        assert expected == actual

def boxes_under_path(p):
    """
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

def is_fragment(frag, path):
    """
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
    """
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
    """
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
    """
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
    """
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
    """
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
    """
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
    """
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
    """
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
    """
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

def compositions(n):
    """
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

def translators(n):
    """
    Return the lperms for which we can compute character values.

    >>> translators(3) == {
    ...     (0, 1, 2): (1, 1, 1),
    ...     (1, 0, 2): (2, 1),
    ...     (0, 2, 1): (2, 1),
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
    """
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

output_preamble = """
from sage.combinat.sf.sfa import zee
R, q = QQ['q'].objgen()
sym = SymmetricFunctions(R)
p = sym.p()
s = sym.s()
h = sym.h()

# Frobenius characters for left and right Hessenberg representations.
# There is one of each for each Dyck path.
# Dyck paths are represented as numbers of missing boxes in each row.
# For example, the following Dyck paths:

#   /\
#  /  \
# /    \
# >>> hess_left[0, 0, 0]
# >>> hess_right[0, 0, 0]

#  /\/\
# /    \
# >>> hess_left[0, 0, 1]
# >>> hess_right[0, 0, 1]

#  /\
# /  \/\
# >>> hess_left[0, 0, 2]
# >>> hess_right[0, 0, 2]

#    /\
# /\/  \
# >>> hess_left[0, 1, 1]
# >>> hess_right[0, 1, 1]

# /\/\/\
# >>> hess_left[0, 1, 2]
# >>> hess_right[0, 1, 2]

#--------------------------------
# conjecture 1
#--------------------------------

# If we remove the grading, then the right Hessenberg representation
# is simply a regular representation.

def eval_q_1(symfunc):
    return symfunc.parent().sum_of_terms(
        (partition, sum(polynomial))
        for partition, polynomial in symfunc
        )

def test_conjecture_1():
    return all(
        p(eval_q_1(hess_right[path])) == p[[1]*len(path)]
        for path in hess_right
        )

#--------------------------------
# conjecture 2
#--------------------------------

# The left Hessenberg representation has an h-positive character
# (so that it should actually be a graded permutation representation).

def is_h_positive(symfunc):
    return all(
        coefficient >= 0
        for partition, polynomial in h(symfunc)
        for coefficient in polynomial
        )

def test_conjecture_2():
    return all(
        is_h_positive(hess_left[path])
        for path in hess_left
        )

#--------------------------------
# conjecture 3
#--------------------------------

# The left and right Hessenberg representations are related by some
# specified Kronecker product, so that the character for (say) the
# left rep can be obtained from the character for the right rep
# by multiplying each coefficient in the power sum basis by
# a rational function of `q` (which only depends on the index of
# the power sum in question).

def f(partition):
    return (1-q)**sum(partition)
F = p.module_morphism(
    codomain=p,
    diagonal=f,
    )

def g(partition):
    return R.prod([(1-q**k) for k in partition])
G = p.module_morphism(
    codomain=p,
    diagonal=g,
    )

def test_conjecture_3():
    return all(
        F(hess_left[path]) == G(hess_right[path])
        for path in hess_left
        )

#--------------------------------
# data
#--------------------------------

hess_right = {}
hess_left = {}

"""

def output_all(sizes):
    start = time.time()
    result = output_preamble
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
    with open('output-12345.py', 'w') as f:
        f.write(output_all([1, 2, 3, 4, 5]))

def long_computation():
    with open('output-123456.py', 'w') as f:
        f.write(output_all([1, 2, 3, 4, 5, 6]))
    with open('output-7.py', 'w') as f:
        f.write(output_all([7]))
    with open('output-8.py', 'w') as f:
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
#    short_computation()
    long_computation()

