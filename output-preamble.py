r"""
Sage script which contains precomputed data on q-chromatic symmetric
functions and q-Frobenius characters of Tymoczko's action on
Hessenberg varieties, along with some code to check certain
conjectures about this data.
"""

from sage.combinat.sf.sfa import zee
R, q = QQ['q'].objgen()
sym = SymmetricFunctions(R)
e = sym.e()
h = sym.h()
m = sym.m()
p = sym.p()
s = sym.s()

# q-chromatic symmetric functions for unit interval orders.
# There is one of each for each Dyck path.

# Frobenius characters for left and right Hessenberg representations.
# There is one of each for each Dyck path.

# Dyck paths are represented as numbers of missing boxes in each downward row.
# For example, the following Dyck paths:

#   /\
#  /\/\
# /\/\/\
# >>> csf[0, 0, 0]
# >>> hess_left[0, 0, 0]
# >>> hess_right[0, 0, 0]

#  /\/\
# /\/\/\
# >>> csf[0, 0, 1]
# >>> hess_left[0, 0, 1]
# >>> hess_right[0, 0, 1]

#  /\
# /\/\/\
# >>> csf[0, 0, 2]
# >>> hess_left[0, 0, 2]
# >>> hess_right[0, 0, 2]

#    /\
# /\/\/\
# >>> csf[0, 1, 1]
# >>> hess_left[0, 1, 1]
# >>> hess_right[0, 1, 1]

# /\/\/\
# >>> csf[0, 1, 2]
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
    r"""
    Test whether ungraded right Hessenberg is the regular representation.
    """
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
    r"""
    Test whether left Hessenberg is h-positive.
    """
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
    r"""
    Test whether left and right Hessenberg are related by Kronecker.
    """
    return all(
        F(hess_left[path]) == G(hess_right[path])
        for path in hess_left
        )

#--------------------------------
# conjecture 4
#--------------------------------

# The q-csf is e-positive
# (so that it should actually be a graded permutation representation).

def is_e_positive(symfunc):
    return all(
        coefficient >= 0
        for partition, polynomial in e(symfunc)
        for coefficient in polynomial
        )

def test_conjecture_4():
    r"""
    Test whether q-csf is e-positive.
    """
    return all(
        is_e_positive(csf[path])
        for path in csf
        )

#--------------------------------
# conjecture 5
#--------------------------------

# The q-csf and the Frobenius character of the left Hessenberg
# should be the same up to fundamental involution.

def test_conjecture_5():
    r"""
    Test whether q-csf and left Hessenberg are related by omega.
    """
    return all(
        s(csf[path]) == s(hess_left[path]).omega()
        for path in csf
        if path in hess_left
        )

#--------------------------------
# data
#--------------------------------

csf = {}
hess_right = {}
hess_left = {}


