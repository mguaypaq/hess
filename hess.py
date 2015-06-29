#!/usr/bin/env python
# -*- coding: utf-8 -*-

r"""
Code for computing characters of Hessenberg varieties.
"""

import itertools as it
import numpy as np
from collections import defaultdict

from fragment import *
from path import *
from perm import *
from util import *

# ---------------------------------------------------------

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

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

def compute_left(path):
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

def compute_right(path):
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
    for ((t, _), coeff) in compute_right(path).iteritems():
        sums[t] += coeff
    result = sums.pop(tuple(range(len(path))))
    if all(s == 0 for s in sums.itervalues()):
        return result
    else:
        return False

# ---------------------------------------------------------

#check translation classes?
#test vector supports?

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
        description='Compute left and right Hessenberg characters for a given Dyck path.',
        )
    parser.add_argument(
        'path',
        help='The Dyck path (e.g. triangle is 000, fully disconnected is 012).',
        )
    args = parser.parse_args()
    path = tuple(map(int, args.path))
    assert is_path(path)
    return path

# ---------------------------------------------------------

left_output_header = r"""hess_left[{path}] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
"""

right_output_header = r"""hess_right[{path}] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
"""

left_output_footer = right_output_footer = r"""    ])

"""

def save(path, left, right):
    n = len(path)
    filename = 'output/hess-' + ''.join(map(str, path)) + '.py'
    with open(filename, 'w') as f:
        f.write(left_output_header.format(path=path))
        tmp = defaultdict(lambda: [0]*(1+n*(n-1)//2))
        for ((lperm, deg), coeff) in left.iteritems():
            tmp[cycle_type[lperm]][deg] = coeff
        for index, coeffs in sorted(tmp.iteritems()):
            while coeffs and coeffs[-1] == 0:
                coeffs.pop()
            if coeffs:
                f.write("    ({}, {}),\n".format(list(index), coeffs))
        f.write(left_output_footer)
        f.write(right_output_header.format(path=path))
        tmp = defaultdict(lambda: [0]*(1+n*(n-1)//2))
        for ((lperm, deg), coeff) in right.iteritems():
            tmp[cycle_type[lperm]][deg] = coeff
        for index, coeffs in sorted(tmp.iteritems()):
            while coeffs and coeffs[-1] == 0:
                coeffs.pop()
            if coeffs:
                f.write("    ({}, {}),\n".format(list(index), coeffs))
        f.write(right_output_footer)

# ---------------------------------------------------------

if __name__ == '__main__':
    doctest()
    setup_logging()
    path = argparse()
    logger.info('starting left computation for path %s', path)
    left = compute_left(path)
    logger.info('starting right computation for path %s', path)
    right = compute_right(path)
    save(path, left, right)
    logger.info('done with path %s', path)

# ---------------------------------------------------------

