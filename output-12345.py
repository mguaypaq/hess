
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

hess_right[(0,)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1], [1]),
    ])

hess_left[(0,)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1], [1]),
    ])

hess_right[(0, 0)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1], [1, 1]),
    ([2], [1, -1]),
    ])

hess_left[(0, 0)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1], [1, 1]),
    ([2], [1, 1]),
    ])

hess_right[(0, 1)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1], [2]),
    ])

hess_left[(0, 1)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1], [2]),
    ])

hess_right[(0, 0, 0)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1], [1, 2, 2, 1]),
    ([2, 1], [1, 0, 0, -1]),
    ([3], [1, -1, -1, 1]),
    ])

hess_left[(0, 0, 0)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1], [1, 2, 2, 1]),
    ([2, 1], [1, 2, 2, 1]),
    ([3], [1, 2, 2, 1]),
    ])

hess_right[(0, 0, 1)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1], [1, 4, 1]),
    ([2, 1], [1, 0, -1]),
    ([3], [1, -2, 1]),
    ])

hess_left[(0, 0, 1)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1], [1, 4, 1]),
    ([2, 1], [1, 2, 1]),
    ([3], [1, 1, 1]),
    ])

hess_right[(0, 0, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1], [3, 3]),
    ([2, 1], [1, -1]),
    ])

hess_left[(0, 0, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1], [3, 3]),
    ([2, 1], [1, 1]),
    ])

hess_right[(0, 1, 1)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1], [3, 3]),
    ([2, 1], [1, -1]),
    ])

hess_left[(0, 1, 1)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1], [3, 3]),
    ([2, 1], [1, 1]),
    ])

hess_right[(0, 1, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1], [6]),
    ])

hess_left[(0, 1, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1], [6]),
    ])

hess_right[(0, 0, 0, 0)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1], [1, 3, 5, 6, 5, 3, 1]),
    ([2, 1, 1], [1, 1, 1, 0, -1, -1, -1]),
    ([2, 2], [1, -1, 1, -2, 1, -1, 1]),
    ([3, 1], [1, 0, -1, 0, -1, 0, 1]),
    ([4], [1, -1, -1, 0, 1, 1, -1]),
    ])

hess_left[(0, 0, 0, 0)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1], [1, 3, 5, 6, 5, 3, 1]),
    ([2, 1, 1], [1, 3, 5, 6, 5, 3, 1]),
    ([2, 2], [1, 3, 5, 6, 5, 3, 1]),
    ([3, 1], [1, 3, 5, 6, 5, 3, 1]),
    ([4], [1, 3, 5, 6, 5, 3, 1]),
    ])

hess_right[(0, 0, 0, 1)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1], [1, 3, 8, 8, 3, 1]),
    ([2, 1, 1], [1, 1, 2, -2, -1, -1]),
    ([2, 2], [1, -1, 0, 0, -1, 1]),
    ([3, 1], [1, 0, -1, -1, 0, 1]),
    ([4], [1, -1, -2, 2, 1, -1]),
    ])

hess_left[(0, 0, 0, 1)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1], [1, 3, 8, 8, 3, 1]),
    ([2, 1, 1], [1, 3, 6, 6, 3, 1]),
    ([2, 2], [1, 3, 4, 4, 3, 1]),
    ([3, 1], [1, 3, 5, 5, 3, 1]),
    ([4], [1, 3, 4, 4, 3, 1]),
    ])

hess_right[(0, 0, 0, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1], [1, 6, 10, 6, 1]),
    ([2, 1, 1], [1, 2, 0, -2, -1]),
    ([2, 2], [1, -2, 2, -2, 1]),
    ([3, 1], [1, 0, -2, 0, 1]),
    ([4], [1, -2, 0, 2, -1]),
    ])

hess_left[(0, 0, 0, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1], [1, 6, 10, 6, 1]),
    ([2, 1, 1], [1, 4, 6, 4, 1]),
    ([2, 2], [1, 2, 2, 2, 1]),
    ([3, 1], [1, 3, 4, 3, 1]),
    ([4], [1, 2, 2, 2, 1]),
    ])

hess_right[(0, 0, 0, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1], [4, 8, 8, 4]),
    ([2, 1, 1], [2, 0, 0, -2]),
    ([3, 1], [1, -1, -1, 1]),
    ])

hess_left[(0, 0, 0, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1], [4, 8, 8, 4]),
    ([2, 1, 1], [2, 4, 4, 2]),
    ([3, 1], [1, 2, 2, 1]),
    ])

hess_right[(0, 0, 1, 1)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1], [1, 6, 10, 6, 1]),
    ([2, 1, 1], [1, 2, 0, -2, -1]),
    ([2, 2], [1, -2, 2, -2, 1]),
    ([3, 1], [1, 0, -2, 0, 1]),
    ([4], [1, -2, 0, 2, -1]),
    ])

hess_left[(0, 0, 1, 1)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1], [1, 6, 10, 6, 1]),
    ([2, 1, 1], [1, 4, 6, 4, 1]),
    ([2, 2], [1, 2, 2, 2, 1]),
    ([3, 1], [1, 3, 4, 3, 1]),
    ([4], [1, 2, 2, 2, 1]),
    ])

hess_right[(0, 0, 1, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1], [1, 11, 11, 1]),
    ([2, 1, 1], [1, 3, -3, -1]),
    ([2, 2], [1, -1, -1, 1]),
    ([3, 1], [1, -1, -1, 1]),
    ([4], [1, -3, 3, -1]),
    ])

hess_left[(0, 0, 1, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1], [1, 11, 11, 1]),
    ([2, 1, 1], [1, 5, 5, 1]),
    ([2, 2], [1, 3, 3, 1]),
    ([3, 1], [1, 2, 2, 1]),
    ([4], [1, 1, 1, 1]),
    ])

hess_right[(0, 0, 1, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1], [4, 16, 4]),
    ([2, 1, 1], [2, 0, -2]),
    ([3, 1], [1, -2, 1]),
    ])

hess_left[(0, 0, 1, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1], [4, 16, 4]),
    ([2, 1, 1], [2, 4, 2]),
    ([3, 1], [1, 1, 1]),
    ])

hess_right[(0, 0, 2, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1], [6, 12, 6]),
    ([2, 1, 1], [2, 0, -2]),
    ([2, 2], [2, -4, 2]),
    ])

hess_left[(0, 0, 2, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1], [6, 12, 6]),
    ([2, 1, 1], [2, 4, 2]),
    ([2, 2], [2, 4, 2]),
    ])

hess_right[(0, 0, 2, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1], [12, 12]),
    ([2, 1, 1], [2, -2]),
    ])

hess_left[(0, 0, 2, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1], [12, 12]),
    ([2, 1, 1], [2, 2]),
    ])

hess_right[(0, 1, 1, 1)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1], [4, 8, 8, 4]),
    ([2, 1, 1], [2, 0, 0, -2]),
    ([3, 1], [1, -1, -1, 1]),
    ])

hess_left[(0, 1, 1, 1)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1], [4, 8, 8, 4]),
    ([2, 1, 1], [2, 4, 4, 2]),
    ([3, 1], [1, 2, 2, 1]),
    ])

hess_right[(0, 1, 1, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1], [4, 16, 4]),
    ([2, 1, 1], [2, 0, -2]),
    ([3, 1], [1, -2, 1]),
    ])

hess_left[(0, 1, 1, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1], [4, 16, 4]),
    ([2, 1, 1], [2, 4, 2]),
    ([3, 1], [1, 1, 1]),
    ])

hess_right[(0, 1, 1, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1], [12, 12]),
    ([2, 1, 1], [2, -2]),
    ])

hess_left[(0, 1, 1, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1], [12, 12]),
    ([2, 1, 1], [2, 2]),
    ])

hess_right[(0, 1, 2, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1], [12, 12]),
    ([2, 1, 1], [2, -2]),
    ])

hess_left[(0, 1, 2, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1], [12, 12]),
    ([2, 1, 1], [2, 2]),
    ])

hess_right[(0, 1, 2, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1], [24]),
    ])

hess_left[(0, 1, 2, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1], [24]),
    ])

hess_right[(0, 0, 0, 0, 0)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 4, 9, 15, 20, 22, 20, 15, 9, 4, 1]),
    ([2, 1, 1, 1], [1, 2, 3, 3, 2, 0, -2, -3, -3, -2, -1]),
    ([2, 2, 1], [1, 0, 1, -1, 0, -2, 0, -1, 1, 0, 1]),
    ([3, 1, 1], [1, 1, 0, 0, -1, -2, -1, 0, 0, 1, 1]),
    ([3, 2], [1, -1, 0, 0, -1, 0, 1, 0, 0, 1, -1]),
    ([4, 1], [1, 0, -1, -1, 0, 0, 0, 1, 1, 0, -1]),
    ([5], [1, -1, -1, 0, 0, 2, 0, 0, -1, -1, 1]),
    ])

hess_left[(0, 0, 0, 0, 0)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 4, 9, 15, 20, 22, 20, 15, 9, 4, 1]),
    ([2, 1, 1, 1], [1, 4, 9, 15, 20, 22, 20, 15, 9, 4, 1]),
    ([2, 2, 1], [1, 4, 9, 15, 20, 22, 20, 15, 9, 4, 1]),
    ([3, 1, 1], [1, 4, 9, 15, 20, 22, 20, 15, 9, 4, 1]),
    ([3, 2], [1, 4, 9, 15, 20, 22, 20, 15, 9, 4, 1]),
    ([4, 1], [1, 4, 9, 15, 20, 22, 20, 15, 9, 4, 1]),
    ([5], [1, 4, 9, 15, 20, 22, 20, 15, 9, 4, 1]),
    ])

hess_right[(0, 0, 0, 0, 1)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 4, 9, 19, 27, 27, 19, 9, 4, 1]),
    ([2, 1, 1, 1], [1, 2, 3, 5, 1, -1, -5, -3, -2, -1]),
    ([2, 2, 1], [1, 0, 1, -1, -1, -1, -1, 1, 0, 1]),
    ([3, 1, 1], [1, 1, 0, 1, -3, -3, 1, 0, 1, 1]),
    ([3, 2], [1, -1, 0, -1, 1, -1, 1, 0, 1, -1]),
    ([4, 1], [1, 0, -1, -1, -1, 1, 1, 1, 0, -1]),
    ([5], [1, -1, -1, -1, 2, 2, -1, -1, -1, 1]),
    ])

hess_left[(0, 0, 0, 0, 1)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 4, 9, 19, 27, 27, 19, 9, 4, 1]),
    ([2, 1, 1, 1], [1, 4, 9, 17, 23, 23, 17, 9, 4, 1]),
    ([2, 2, 1], [1, 4, 9, 15, 19, 19, 15, 9, 4, 1]),
    ([3, 1, 1], [1, 4, 9, 16, 21, 21, 16, 9, 4, 1]),
    ([3, 2], [1, 4, 9, 14, 17, 17, 14, 9, 4, 1]),
    ([4, 1], [1, 4, 9, 15, 19, 19, 15, 9, 4, 1]),
    ([5], [1, 4, 9, 14, 17, 17, 14, 9, 4, 1]),
    ])

hess_right[(0, 0, 0, 0, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 4, 13, 26, 32, 26, 13, 4, 1]),
    ([2, 1, 1, 1], [1, 2, 5, 4, 0, -4, -5, -2, -1]),
    ([2, 2, 1], [1, 0, 1, -2, 0, -2, 1, 0, 1]),
    ([3, 1, 1], [1, 1, 1, -1, -4, -1, 1, 1, 1]),
    ([3, 2], [1, -1, -1, 1, 0, -1, 1, 1, -1]),
    ([4, 1], [1, 0, -1, -2, 0, 2, 1, 0, -1]),
    ([5], [1, -1, -2, 1, 2, 1, -2, -1, 1]),
    ])

hess_left[(0, 0, 0, 0, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 4, 13, 26, 32, 26, 13, 4, 1]),
    ([2, 1, 1, 1], [1, 4, 11, 20, 24, 20, 11, 4, 1]),
    ([2, 2, 1], [1, 4, 9, 14, 16, 14, 9, 4, 1]),
    ([3, 1, 1], [1, 4, 10, 17, 20, 17, 10, 4, 1]),
    ([3, 2], [1, 4, 8, 11, 12, 11, 8, 4, 1]),
    ([4, 1], [1, 4, 9, 14, 16, 14, 9, 4, 1]),
    ([5], [1, 4, 8, 11, 12, 11, 8, 4, 1]),
    ])

hess_right[(0, 0, 0, 0, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 8, 20, 31, 31, 20, 8, 1]),
    ([2, 1, 1, 1], [1, 4, 4, 3, -3, -4, -4, -1]),
    ([2, 2, 1], [1, 0, 0, -1, -1, 0, 0, 1]),
    ([3, 1, 1], [1, 2, -1, -2, -2, -1, 2, 1]),
    ([3, 2], [1, -2, 1, 0, 0, -1, 2, -1]),
    ([4, 1], [1, 0, -2, -1, 1, 2, 0, -1]),
    ([5], [1, -2, 0, 1, 1, 0, -2, 1]),
    ])

hess_left[(0, 0, 0, 0, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 8, 20, 31, 31, 20, 8, 1]),
    ([2, 1, 1, 1], [1, 6, 14, 21, 21, 14, 6, 1]),
    ([2, 2, 1], [1, 4, 8, 11, 11, 8, 4, 1]),
    ([3, 1, 1], [1, 5, 11, 16, 16, 11, 5, 1]),
    ([3, 2], [1, 3, 5, 6, 6, 5, 3, 1]),
    ([4, 1], [1, 4, 8, 11, 11, 8, 4, 1]),
    ([5], [1, 3, 5, 6, 6, 5, 3, 1]),
    ])

hess_right[(0, 0, 0, 0, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [5, 15, 25, 30, 25, 15, 5]),
    ([2, 1, 1, 1], [3, 3, 3, 0, -3, -3, -3]),
    ([2, 2, 1], [1, -1, 1, -2, 1, -1, 1]),
    ([3, 1, 1], [2, 0, -2, 0, -2, 0, 2]),
    ([4, 1], [1, -1, -1, 0, 1, 1, -1]),
    ])

hess_left[(0, 0, 0, 0, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [5, 15, 25, 30, 25, 15, 5]),
    ([2, 1, 1, 1], [3, 9, 15, 18, 15, 9, 3]),
    ([2, 2, 1], [1, 3, 5, 6, 5, 3, 1]),
    ([3, 1, 1], [2, 6, 10, 12, 10, 6, 2]),
    ([4, 1], [1, 3, 5, 6, 5, 3, 1]),
    ])

hess_right[(0, 0, 0, 1, 1)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 4, 13, 26, 32, 26, 13, 4, 1]),
    ([2, 1, 1, 1], [1, 2, 5, 4, 0, -4, -5, -2, -1]),
    ([2, 2, 1], [1, 0, 1, -2, 0, -2, 1, 0, 1]),
    ([3, 1, 1], [1, 1, 1, -1, -4, -1, 1, 1, 1]),
    ([3, 2], [1, -1, -1, 1, 0, -1, 1, 1, -1]),
    ([4, 1], [1, 0, -1, -2, 0, 2, 1, 0, -1]),
    ([5], [1, -1, -2, 1, 2, 1, -2, -1, 1]),
    ])

hess_left[(0, 0, 0, 1, 1)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 4, 13, 26, 32, 26, 13, 4, 1]),
    ([2, 1, 1, 1], [1, 4, 11, 20, 24, 20, 11, 4, 1]),
    ([2, 2, 1], [1, 4, 9, 14, 16, 14, 9, 4, 1]),
    ([3, 1, 1], [1, 4, 10, 17, 20, 17, 10, 4, 1]),
    ([3, 2], [1, 4, 8, 11, 12, 11, 8, 4, 1]),
    ([4, 1], [1, 4, 9, 14, 16, 14, 9, 4, 1]),
    ([5], [1, 4, 8, 11, 12, 11, 8, 4, 1]),
    ])

hess_right[(0, 0, 0, 1, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 4, 17, 38, 38, 17, 4, 1]),
    ([2, 1, 1, 1], [1, 2, 7, 4, -4, -7, -2, -1]),
    ([2, 2, 1], [1, 0, 1, -2, -2, 1, 0, 1]),
    ([3, 1, 1], [1, 1, 2, -4, -4, 2, 1, 1]),
    ([3, 2], [1, -1, -2, 4, -4, 2, 1, -1]),
    ([4, 1], [1, 0, -1, -4, 4, 1, 0, -1]),
    ([5], [1, -1, -3, 3, 3, -3, -1, 1]),
    ])

hess_left[(0, 0, 0, 1, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 4, 17, 38, 38, 17, 4, 1]),
    ([2, 1, 1, 1], [1, 4, 13, 24, 24, 13, 4, 1]),
    ([2, 2, 1], [1, 4, 9, 14, 14, 9, 4, 1]),
    ([3, 1, 1], [1, 4, 11, 17, 17, 11, 4, 1]),
    ([3, 2], [1, 4, 7, 9, 9, 7, 4, 1]),
    ([4, 1], [1, 4, 9, 12, 12, 9, 4, 1]),
    ([5], [1, 4, 7, 8, 8, 7, 4, 1]),
    ])

hess_right[(0, 0, 0, 1, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 8, 29, 44, 29, 8, 1]),
    ([2, 1, 1, 1], [1, 4, 7, 0, -7, -4, -1]),
    ([2, 2, 1], [1, 0, 1, -4, 1, 0, 1]),
    ([3, 1, 1], [1, 2, -1, -4, -1, 2, 1]),
    ([3, 2], [1, -2, 1, 0, -1, 2, -1]),
    ([4, 1], [1, 0, -3, 0, 3, 0, -1]),
    ([5], [1, -2, -1, 4, -1, -2, 1]),
    ])

hess_left[(0, 0, 0, 1, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 8, 29, 44, 29, 8, 1]),
    ([2, 1, 1, 1], [1, 6, 17, 24, 17, 6, 1]),
    ([2, 2, 1], [1, 4, 9, 12, 9, 4, 1]),
    ([3, 1, 1], [1, 5, 11, 14, 11, 5, 1]),
    ([3, 2], [1, 3, 5, 6, 5, 3, 1]),
    ([4, 1], [1, 4, 7, 8, 7, 4, 1]),
    ([5], [1, 3, 4, 4, 4, 3, 1]),
    ])

hess_right[(0, 0, 0, 1, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [5, 15, 40, 40, 15, 5]),
    ([2, 1, 1, 1], [3, 3, 6, -6, -3, -3]),
    ([2, 2, 1], [1, -1, 0, 0, -1, 1]),
    ([3, 1, 1], [2, 0, -2, -2, 0, 2]),
    ([4, 1], [1, -1, -2, 2, 1, -1]),
    ])

hess_left[(0, 0, 0, 1, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [5, 15, 40, 40, 15, 5]),
    ([2, 1, 1, 1], [3, 9, 18, 18, 9, 3]),
    ([2, 2, 1], [1, 3, 4, 4, 3, 1]),
    ([3, 1, 1], [2, 6, 10, 10, 6, 2]),
    ([4, 1], [1, 3, 4, 4, 3, 1]),
    ])

hess_right[(0, 0, 0, 2, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 8, 29, 44, 29, 8, 1]),
    ([2, 1, 1, 1], [1, 4, 7, 0, -7, -4, -1]),
    ([2, 2, 1], [1, 0, 1, -4, 1, 0, 1]),
    ([3, 1, 1], [1, 2, -1, -4, -1, 2, 1]),
    ([3, 2], [1, -2, 1, 0, -1, 2, -1]),
    ([4, 1], [1, 0, -3, 0, 3, 0, -1]),
    ([5], [1, -2, -1, 4, -1, -2, 1]),
    ])

hess_left[(0, 0, 0, 2, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 8, 29, 44, 29, 8, 1]),
    ([2, 1, 1, 1], [1, 6, 17, 24, 17, 6, 1]),
    ([2, 2, 1], [1, 4, 9, 12, 9, 4, 1]),
    ([3, 1, 1], [1, 5, 11, 14, 11, 5, 1]),
    ([3, 2], [1, 3, 5, 6, 5, 3, 1]),
    ([4, 1], [1, 4, 7, 8, 7, 4, 1]),
    ([5], [1, 3, 4, 4, 4, 3, 1]),
    ])

hess_right[(0, 0, 0, 2, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 17, 42, 42, 17, 1]),
    ([2, 1, 1, 1], [1, 7, 4, -4, -7, -1]),
    ([2, 2, 1], [1, 1, -2, -2, 1, 1]),
    ([3, 1, 1], [1, 2, -3, -3, 2, 1]),
    ([3, 2], [1, -2, 1, -1, 2, -1]),
    ([4, 1], [1, -1, -2, 2, 1, -1]),
    ([5], [1, -3, 2, 2, -3, 1]),
    ])

hess_left[(0, 0, 0, 2, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 17, 42, 42, 17, 1]),
    ([2, 1, 1, 1], [1, 9, 20, 20, 9, 1]),
    ([2, 2, 1], [1, 5, 10, 10, 5, 1]),
    ([3, 1, 1], [1, 5, 9, 9, 5, 1]),
    ([3, 2], [1, 3, 5, 5, 3, 1]),
    ([4, 1], [1, 3, 4, 4, 3, 1]),
    ([5], [1, 2, 2, 2, 2, 1]),
    ])

hess_right[(0, 0, 0, 2, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [5, 30, 50, 30, 5]),
    ([2, 1, 1, 1], [3, 6, 0, -6, -3]),
    ([2, 2, 1], [1, -2, 2, -2, 1]),
    ([3, 1, 1], [2, 0, -4, 0, 2]),
    ([4, 1], [1, -2, 0, 2, -1]),
    ])

hess_left[(0, 0, 0, 2, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [5, 30, 50, 30, 5]),
    ([2, 1, 1, 1], [3, 12, 18, 12, 3]),
    ([2, 2, 1], [1, 2, 2, 2, 1]),
    ([3, 1, 1], [2, 6, 8, 6, 2]),
    ([4, 1], [1, 2, 2, 2, 1]),
    ])

hess_right[(0, 0, 0, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [10, 30, 40, 30, 10]),
    ([2, 1, 1, 1], [4, 4, 0, -4, -4]),
    ([2, 2, 1], [2, -2, 0, -2, 2]),
    ([3, 1, 1], [1, 0, -2, 0, 1]),
    ([3, 2], [1, -2, 0, 2, -1]),
    ])

hess_left[(0, 0, 0, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [10, 30, 40, 30, 10]),
    ([2, 1, 1, 1], [4, 12, 16, 12, 4]),
    ([2, 2, 1], [2, 6, 8, 6, 2]),
    ([3, 1, 1], [1, 3, 4, 3, 1]),
    ([3, 2], [1, 3, 4, 3, 1]),
    ])

hess_right[(0, 0, 0, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [20, 40, 40, 20]),
    ([2, 1, 1, 1], [6, 0, 0, -6]),
    ([3, 1, 1], [2, -2, -2, 2]),
    ])

hess_left[(0, 0, 0, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [20, 40, 40, 20]),
    ([2, 1, 1, 1], [6, 12, 12, 6]),
    ([3, 1, 1], [2, 4, 4, 2]),
    ])

hess_right[(0, 0, 1, 1, 1)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 8, 20, 31, 31, 20, 8, 1]),
    ([2, 1, 1, 1], [1, 4, 4, 3, -3, -4, -4, -1]),
    ([2, 2, 1], [1, 0, 0, -1, -1, 0, 0, 1]),
    ([3, 1, 1], [1, 2, -1, -2, -2, -1, 2, 1]),
    ([3, 2], [1, -2, 1, 0, 0, -1, 2, -1]),
    ([4, 1], [1, 0, -2, -1, 1, 2, 0, -1]),
    ([5], [1, -2, 0, 1, 1, 0, -2, 1]),
    ])

hess_left[(0, 0, 1, 1, 1)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 8, 20, 31, 31, 20, 8, 1]),
    ([2, 1, 1, 1], [1, 6, 14, 21, 21, 14, 6, 1]),
    ([2, 2, 1], [1, 4, 8, 11, 11, 8, 4, 1]),
    ([3, 1, 1], [1, 5, 11, 16, 16, 11, 5, 1]),
    ([3, 2], [1, 3, 5, 6, 6, 5, 3, 1]),
    ([4, 1], [1, 4, 8, 11, 11, 8, 4, 1]),
    ([5], [1, 3, 5, 6, 6, 5, 3, 1]),
    ])

hess_right[(0, 0, 1, 1, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 8, 29, 44, 29, 8, 1]),
    ([2, 1, 1, 1], [1, 4, 7, 0, -7, -4, -1]),
    ([2, 2, 1], [1, 0, 1, -4, 1, 0, 1]),
    ([3, 1, 1], [1, 2, -1, -4, -1, 2, 1]),
    ([3, 2], [1, -2, 1, 0, -1, 2, -1]),
    ([4, 1], [1, 0, -3, 0, 3, 0, -1]),
    ([5], [1, -2, -1, 4, -1, -2, 1]),
    ])

hess_left[(0, 0, 1, 1, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 8, 29, 44, 29, 8, 1]),
    ([2, 1, 1, 1], [1, 6, 17, 24, 17, 6, 1]),
    ([2, 2, 1], [1, 4, 9, 12, 9, 4, 1]),
    ([3, 1, 1], [1, 5, 11, 14, 11, 5, 1]),
    ([3, 2], [1, 3, 5, 6, 5, 3, 1]),
    ([4, 1], [1, 4, 7, 8, 7, 4, 1]),
    ([5], [1, 3, 4, 4, 4, 3, 1]),
    ])

hess_right[(0, 0, 1, 1, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 12, 47, 47, 12, 1]),
    ([2, 1, 1, 1], [1, 6, 7, -7, -6, -1]),
    ([2, 2, 1], [1, 0, -1, -1, 0, 1]),
    ([3, 1, 1], [1, 3, -4, -4, 3, 1]),
    ([3, 2], [1, -3, 4, -4, 3, -1]),
    ([4, 1], [1, 0, -5, 5, 0, -1]),
    ([5], [1, -3, 2, 2, -3, 1]),
    ])

hess_left[(0, 0, 1, 1, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 12, 47, 47, 12, 1]),
    ([2, 1, 1, 1], [1, 8, 21, 21, 8, 1]),
    ([2, 2, 1], [1, 4, 7, 7, 4, 1]),
    ([3, 1, 1], [1, 6, 11, 11, 6, 1]),
    ([3, 2], [1, 2, 3, 3, 2, 1]),
    ([4, 1], [1, 4, 5, 5, 4, 1]),
    ([5], [1, 2, 2, 2, 2, 1]),
    ])

hess_right[(0, 0, 1, 1, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [5, 30, 50, 30, 5]),
    ([2, 1, 1, 1], [3, 6, 0, -6, -3]),
    ([2, 2, 1], [1, -2, 2, -2, 1]),
    ([3, 1, 1], [2, 0, -4, 0, 2]),
    ([4, 1], [1, -2, 0, 2, -1]),
    ])

hess_left[(0, 0, 1, 1, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [5, 30, 50, 30, 5]),
    ([2, 1, 1, 1], [3, 12, 18, 12, 3]),
    ([2, 2, 1], [1, 2, 2, 2, 1]),
    ([3, 1, 1], [2, 6, 8, 6, 2]),
    ([4, 1], [1, 2, 2, 2, 1]),
    ])

hess_right[(0, 0, 1, 2, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 17, 42, 42, 17, 1]),
    ([2, 1, 1, 1], [1, 7, 4, -4, -7, -1]),
    ([2, 2, 1], [1, 1, -2, -2, 1, 1]),
    ([3, 1, 1], [1, 2, -3, -3, 2, 1]),
    ([3, 2], [1, -2, 1, -1, 2, -1]),
    ([4, 1], [1, -1, -2, 2, 1, -1]),
    ([5], [1, -3, 2, 2, -3, 1]),
    ])

hess_left[(0, 0, 1, 2, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 17, 42, 42, 17, 1]),
    ([2, 1, 1, 1], [1, 9, 20, 20, 9, 1]),
    ([2, 2, 1], [1, 5, 10, 10, 5, 1]),
    ([3, 1, 1], [1, 5, 9, 9, 5, 1]),
    ([3, 2], [1, 3, 5, 5, 3, 1]),
    ([4, 1], [1, 3, 4, 4, 3, 1]),
    ([5], [1, 2, 2, 2, 2, 1]),
    ])

hess_right[(0, 0, 1, 2, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 26, 66, 26, 1]),
    ([2, 1, 1, 1], [1, 10, 0, -10, -1]),
    ([2, 2, 1], [1, 2, -6, 2, 1]),
    ([3, 1, 1], [1, 2, -6, 2, 1]),
    ([3, 2], [1, -2, 0, 2, -1]),
    ([4, 1], [1, -2, 0, 2, -1]),
    ([5], [1, -4, 6, -4, 1]),
    ])

hess_left[(0, 0, 1, 2, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 26, 66, 26, 1]),
    ([2, 1, 1, 1], [1, 12, 22, 12, 1]),
    ([2, 2, 1], [1, 6, 10, 6, 1]),
    ([3, 1, 1], [1, 5, 6, 5, 1]),
    ([3, 2], [1, 3, 4, 3, 1]),
    ([4, 1], [1, 2, 2, 2, 1]),
    ([5], [1, 1, 1, 1, 1]),
    ])

hess_right[(0, 0, 1, 2, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [5, 55, 55, 5]),
    ([2, 1, 1, 1], [3, 9, -9, -3]),
    ([2, 2, 1], [1, -1, -1, 1]),
    ([3, 1, 1], [2, -2, -2, 2]),
    ([4, 1], [1, -3, 3, -1]),
    ])

hess_left[(0, 0, 1, 2, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [5, 55, 55, 5]),
    ([2, 1, 1, 1], [3, 15, 15, 3]),
    ([2, 2, 1], [1, 3, 3, 1]),
    ([3, 1, 1], [2, 4, 4, 2]),
    ([4, 1], [1, 1, 1, 1]),
    ])

hess_right[(0, 0, 1, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [10, 50, 50, 10]),
    ([2, 1, 1, 1], [4, 6, -6, -4]),
    ([2, 2, 1], [2, -2, -2, 2]),
    ([3, 1, 1], [1, -1, -1, 1]),
    ([3, 2], [1, -3, 3, -1]),
    ])

hess_left[(0, 0, 1, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [10, 50, 50, 10]),
    ([2, 1, 1, 1], [4, 14, 14, 4]),
    ([2, 2, 1], [2, 6, 6, 2]),
    ([3, 1, 1], [1, 2, 2, 1]),
    ([3, 2], [1, 2, 2, 1]),
    ])

hess_right[(0, 0, 1, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [20, 80, 20]),
    ([2, 1, 1, 1], [6, 0, -6]),
    ([3, 1, 1], [2, -4, 2]),
    ])

hess_left[(0, 0, 1, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [20, 80, 20]),
    ([2, 1, 1, 1], [6, 12, 6]),
    ([3, 1, 1], [2, 2, 2]),
    ])

hess_right[(0, 0, 2, 2, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [10, 30, 40, 30, 10]),
    ([2, 1, 1, 1], [4, 4, 0, -4, -4]),
    ([2, 2, 1], [2, -2, 0, -2, 2]),
    ([3, 1, 1], [1, 0, -2, 0, 1]),
    ([3, 2], [1, -2, 0, 2, -1]),
    ])

hess_left[(0, 0, 2, 2, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [10, 30, 40, 30, 10]),
    ([2, 1, 1, 1], [4, 12, 16, 12, 4]),
    ([2, 2, 1], [2, 6, 8, 6, 2]),
    ([3, 1, 1], [1, 3, 4, 3, 1]),
    ([3, 2], [1, 3, 4, 3, 1]),
    ])

hess_right[(0, 0, 2, 2, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [10, 50, 50, 10]),
    ([2, 1, 1, 1], [4, 6, -6, -4]),
    ([2, 2, 1], [2, -2, -2, 2]),
    ([3, 1, 1], [1, -1, -1, 1]),
    ([3, 2], [1, -3, 3, -1]),
    ])

hess_left[(0, 0, 2, 2, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [10, 50, 50, 10]),
    ([2, 1, 1, 1], [4, 14, 14, 4]),
    ([2, 2, 1], [2, 6, 6, 2]),
    ([3, 1, 1], [1, 2, 2, 1]),
    ([3, 2], [1, 2, 2, 1]),
    ])

hess_right[(0, 0, 2, 2, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [30, 60, 30]),
    ([2, 1, 1, 1], [6, 0, -6]),
    ([2, 2, 1], [2, -4, 2]),
    ])

hess_left[(0, 0, 2, 2, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [30, 60, 30]),
    ([2, 1, 1, 1], [6, 12, 6]),
    ([2, 2, 1], [2, 4, 2]),
    ])

hess_right[(0, 0, 2, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [30, 60, 30]),
    ([2, 1, 1, 1], [6, 0, -6]),
    ([2, 2, 1], [2, -4, 2]),
    ])

hess_left[(0, 0, 2, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [30, 60, 30]),
    ([2, 1, 1, 1], [6, 12, 6]),
    ([2, 2, 1], [2, 4, 2]),
    ])

hess_right[(0, 0, 2, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [60, 60]),
    ([2, 1, 1, 1], [6, -6]),
    ])

hess_left[(0, 0, 2, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [60, 60]),
    ([2, 1, 1, 1], [6, 6]),
    ])

hess_right[(0, 1, 1, 1, 1)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [5, 15, 25, 30, 25, 15, 5]),
    ([2, 1, 1, 1], [3, 3, 3, 0, -3, -3, -3]),
    ([2, 2, 1], [1, -1, 1, -2, 1, -1, 1]),
    ([3, 1, 1], [2, 0, -2, 0, -2, 0, 2]),
    ([4, 1], [1, -1, -1, 0, 1, 1, -1]),
    ])

hess_left[(0, 1, 1, 1, 1)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [5, 15, 25, 30, 25, 15, 5]),
    ([2, 1, 1, 1], [3, 9, 15, 18, 15, 9, 3]),
    ([2, 2, 1], [1, 3, 5, 6, 5, 3, 1]),
    ([3, 1, 1], [2, 6, 10, 12, 10, 6, 2]),
    ([4, 1], [1, 3, 5, 6, 5, 3, 1]),
    ])

hess_right[(0, 1, 1, 1, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [5, 15, 40, 40, 15, 5]),
    ([2, 1, 1, 1], [3, 3, 6, -6, -3, -3]),
    ([2, 2, 1], [1, -1, 0, 0, -1, 1]),
    ([3, 1, 1], [2, 0, -2, -2, 0, 2]),
    ([4, 1], [1, -1, -2, 2, 1, -1]),
    ])

hess_left[(0, 1, 1, 1, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [5, 15, 40, 40, 15, 5]),
    ([2, 1, 1, 1], [3, 9, 18, 18, 9, 3]),
    ([2, 2, 1], [1, 3, 4, 4, 3, 1]),
    ([3, 1, 1], [2, 6, 10, 10, 6, 2]),
    ([4, 1], [1, 3, 4, 4, 3, 1]),
    ])

hess_right[(0, 1, 1, 1, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [5, 30, 50, 30, 5]),
    ([2, 1, 1, 1], [3, 6, 0, -6, -3]),
    ([2, 2, 1], [1, -2, 2, -2, 1]),
    ([3, 1, 1], [2, 0, -4, 0, 2]),
    ([4, 1], [1, -2, 0, 2, -1]),
    ])

hess_left[(0, 1, 1, 1, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [5, 30, 50, 30, 5]),
    ([2, 1, 1, 1], [3, 12, 18, 12, 3]),
    ([2, 2, 1], [1, 2, 2, 2, 1]),
    ([3, 1, 1], [2, 6, 8, 6, 2]),
    ([4, 1], [1, 2, 2, 2, 1]),
    ])

hess_right[(0, 1, 1, 1, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [20, 40, 40, 20]),
    ([2, 1, 1, 1], [6, 0, 0, -6]),
    ([3, 1, 1], [2, -2, -2, 2]),
    ])

hess_left[(0, 1, 1, 1, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [20, 40, 40, 20]),
    ([2, 1, 1, 1], [6, 12, 12, 6]),
    ([3, 1, 1], [2, 4, 4, 2]),
    ])

hess_right[(0, 1, 1, 2, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [5, 30, 50, 30, 5]),
    ([2, 1, 1, 1], [3, 6, 0, -6, -3]),
    ([2, 2, 1], [1, -2, 2, -2, 1]),
    ([3, 1, 1], [2, 0, -4, 0, 2]),
    ([4, 1], [1, -2, 0, 2, -1]),
    ])

hess_left[(0, 1, 1, 2, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [5, 30, 50, 30, 5]),
    ([2, 1, 1, 1], [3, 12, 18, 12, 3]),
    ([2, 2, 1], [1, 2, 2, 2, 1]),
    ([3, 1, 1], [2, 6, 8, 6, 2]),
    ([4, 1], [1, 2, 2, 2, 1]),
    ])

hess_right[(0, 1, 1, 2, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [5, 55, 55, 5]),
    ([2, 1, 1, 1], [3, 9, -9, -3]),
    ([2, 2, 1], [1, -1, -1, 1]),
    ([3, 1, 1], [2, -2, -2, 2]),
    ([4, 1], [1, -3, 3, -1]),
    ])

hess_left[(0, 1, 1, 2, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [5, 55, 55, 5]),
    ([2, 1, 1, 1], [3, 15, 15, 3]),
    ([2, 2, 1], [1, 3, 3, 1]),
    ([3, 1, 1], [2, 4, 4, 2]),
    ([4, 1], [1, 1, 1, 1]),
    ])

hess_right[(0, 1, 1, 2, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [20, 80, 20]),
    ([2, 1, 1, 1], [6, 0, -6]),
    ([3, 1, 1], [2, -4, 2]),
    ])

hess_left[(0, 1, 1, 2, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [20, 80, 20]),
    ([2, 1, 1, 1], [6, 12, 6]),
    ([3, 1, 1], [2, 2, 2]),
    ])

hess_right[(0, 1, 1, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [30, 60, 30]),
    ([2, 1, 1, 1], [6, 0, -6]),
    ([2, 2, 1], [2, -4, 2]),
    ])

hess_left[(0, 1, 1, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [30, 60, 30]),
    ([2, 1, 1, 1], [6, 12, 6]),
    ([2, 2, 1], [2, 4, 2]),
    ])

hess_right[(0, 1, 1, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [60, 60]),
    ([2, 1, 1, 1], [6, -6]),
    ])

hess_left[(0, 1, 1, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [60, 60]),
    ([2, 1, 1, 1], [6, 6]),
    ])

hess_right[(0, 1, 2, 2, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [20, 40, 40, 20]),
    ([2, 1, 1, 1], [6, 0, 0, -6]),
    ([3, 1, 1], [2, -2, -2, 2]),
    ])

hess_left[(0, 1, 2, 2, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [20, 40, 40, 20]),
    ([2, 1, 1, 1], [6, 12, 12, 6]),
    ([3, 1, 1], [2, 4, 4, 2]),
    ])

hess_right[(0, 1, 2, 2, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [20, 80, 20]),
    ([2, 1, 1, 1], [6, 0, -6]),
    ([3, 1, 1], [2, -4, 2]),
    ])

hess_left[(0, 1, 2, 2, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [20, 80, 20]),
    ([2, 1, 1, 1], [6, 12, 6]),
    ([3, 1, 1], [2, 2, 2]),
    ])

hess_right[(0, 1, 2, 2, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [60, 60]),
    ([2, 1, 1, 1], [6, -6]),
    ])

hess_left[(0, 1, 2, 2, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [60, 60]),
    ([2, 1, 1, 1], [6, 6]),
    ])

hess_right[(0, 1, 2, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [60, 60]),
    ([2, 1, 1, 1], [6, -6]),
    ])

hess_left[(0, 1, 2, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [60, 60]),
    ([2, 1, 1, 1], [6, 6]),
    ])

hess_right[(0, 1, 2, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [120]),
    ])

hess_left[(0, 1, 2, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [120]),
    ])

