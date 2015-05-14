from sage.combinat.sf.sfa import zee
R, q = QQ['q'].objgen()
sym = SymmetricFunctions(R)
p = sym.p()
s = sym.s()
h = sym.h()

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

hess_right[(0, 0, 0, 0, 0, 0)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 5, 14, 29, 49, 71, 90, 101, 101, 90, 71, 49, 29, 14, 5, 1]),
    ([2, 1, 1, 1, 1], [1, 3, 6, 9, 11, 11, 8, 3, -3, -8, -11, -11, -9, -6, -3, -1]),
    ([2, 2, 1, 1], [1, 1, 2, 1, 1, -1, -2, -3, -3, -2, -1, 1, 1, 2, 1, 1]),
    ([2, 2, 2], [1, -1, 2, -3, 3, -5, 4, -5, 5, -4, 5, -3, 3, -2, 1, -1]),
    ([3, 1, 1, 1], [1, 2, 2, 2, 1, -1, -3, -4, -4, -3, -1, 1, 2, 2, 2, 1]),
    ([3, 2, 1], [1, 0, 0, 0, -1, -1, -1, 0, 0, 1, 1, 1, 0, 0, 0, -1]),
    ([3, 3], [1, -1, -1, 2, -2, -1, 3, -1, -1, 3, -1, -2, 2, -1, -1, 1]),
    ([4, 1, 1], [1, 1, 0, -1, -1, -1, -2, -1, 1, 2, 1, 1, 1, 0, -1, -1]),
    ([4, 2], [1, -1, 0, -1, 1, -1, 0, 1, 1, 0, -1, 1, -1, 0, -1, 1]),
    ([5, 1], [1, 0, -1, -1, -1, 1, 0, 1, 1, 0, 1, -1, -1, -1, 0, 1]),
    ([6], [1, -1, -1, 0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 1, 1, -1]),
    ])

hess_left[(0, 0, 0, 0, 0, 0)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 5, 14, 29, 49, 71, 90, 101, 101, 90, 71, 49, 29, 14, 5, 1]),
    ([2, 1, 1, 1, 1], [1, 5, 14, 29, 49, 71, 90, 101, 101, 90, 71, 49, 29, 14, 5, 1]),
    ([2, 2, 1, 1], [1, 5, 14, 29, 49, 71, 90, 101, 101, 90, 71, 49, 29, 14, 5, 1]),
    ([2, 2, 2], [1, 5, 14, 29, 49, 71, 90, 101, 101, 90, 71, 49, 29, 14, 5, 1]),
    ([3, 1, 1, 1], [1, 5, 14, 29, 49, 71, 90, 101, 101, 90, 71, 49, 29, 14, 5, 1]),
    ([3, 2, 1], [1, 5, 14, 29, 49, 71, 90, 101, 101, 90, 71, 49, 29, 14, 5, 1]),
    ([3, 3], [1, 5, 14, 29, 49, 71, 90, 101, 101, 90, 71, 49, 29, 14, 5, 1]),
    ([4, 1, 1], [1, 5, 14, 29, 49, 71, 90, 101, 101, 90, 71, 49, 29, 14, 5, 1]),
    ([4, 2], [1, 5, 14, 29, 49, 71, 90, 101, 101, 90, 71, 49, 29, 14, 5, 1]),
    ([5, 1], [1, 5, 14, 29, 49, 71, 90, 101, 101, 90, 71, 49, 29, 14, 5, 1]),
    ([6], [1, 5, 14, 29, 49, 71, 90, 101, 101, 90, 71, 49, 29, 14, 5, 1]),
    ])

hess_right[(0, 0, 0, 0, 0, 1)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 5, 14, 29, 54, 85, 111, 122, 111, 85, 54, 29, 14, 5, 1]),
    ([2, 1, 1, 1, 1], [1, 3, 6, 9, 14, 13, 9, 0, -9, -13, -14, -9, -6, -3, -1]),
    ([2, 2, 1, 1], [1, 1, 2, 1, 2, -3, -1, -6, -1, -3, 2, 1, 2, 1, 1]),
    ([2, 2, 2], [1, -1, 2, -3, 2, -3, 1, 0, -1, 3, -2, 3, -2, 1, -1]),
    ([3, 1, 1, 1], [1, 2, 2, 2, 3, -2, -6, -4, -6, -2, 3, 2, 2, 2, 1]),
    ([3, 2, 1], [1, 0, 0, 0, -1, -2, 0, 0, 0, 2, 1, 0, 0, 0, -1]),
    ([3, 3], [1, -1, -1, 2, -3, 1, 3, -4, 3, 1, -3, 2, -1, -1, 1]),
    ([4, 1, 1], [1, 1, 0, -1, 0, -3, -3, 0, 3, 3, 0, 1, 0, -1, -1]),
    ([4, 2], [1, -1, 0, -1, 0, 1, -1, 2, -1, 1, 0, -1, 0, -1, 1]),
    ([5, 1], [1, 0, -1, -1, -1, 0, 1, 2, 1, 0, -1, -1, -1, 0, 1]),
    ([6], [1, -1, -1, 0, -1, 3, 1, 0, -1, -3, 1, 0, 1, 1, -1]),
    ])

hess_left[(0, 0, 0, 0, 0, 1)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 5, 14, 29, 54, 85, 111, 122, 111, 85, 54, 29, 14, 5, 1]),
    ([2, 1, 1, 1, 1], [1, 5, 14, 29, 52, 79, 101, 110, 101, 79, 52, 29, 14, 5, 1]),
    ([2, 2, 1, 1], [1, 5, 14, 29, 50, 73, 91, 98, 91, 73, 50, 29, 14, 5, 1]),
    ([2, 2, 2], [1, 5, 14, 29, 48, 67, 81, 86, 81, 67, 48, 29, 14, 5, 1]),
    ([3, 1, 1, 1], [1, 5, 14, 29, 51, 76, 96, 104, 96, 76, 51, 29, 14, 5, 1]),
    ([3, 2, 1], [1, 5, 14, 29, 49, 70, 86, 92, 86, 70, 49, 29, 14, 5, 1]),
    ([3, 3], [1, 5, 14, 29, 48, 67, 81, 86, 81, 67, 48, 29, 14, 5, 1]),
    ([4, 1, 1], [1, 5, 14, 29, 50, 73, 91, 98, 91, 73, 50, 29, 14, 5, 1]),
    ([4, 2], [1, 5, 14, 29, 48, 67, 81, 86, 81, 67, 48, 29, 14, 5, 1]),
    ([5, 1], [1, 5, 14, 29, 49, 70, 86, 92, 86, 70, 49, 29, 14, 5, 1]),
    ([6], [1, 5, 14, 29, 48, 67, 81, 86, 81, 67, 48, 29, 14, 5, 1]),
    ])

hess_right[(0, 0, 0, 0, 0, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 5, 14, 34, 68, 106, 132, 132, 106, 68, 34, 14, 5, 1]),
    ([2, 1, 1, 1, 1], [1, 3, 6, 12, 16, 14, 6, -6, -14, -16, -12, -6, -3, -1]),
    ([2, 2, 1, 1], [1, 1, 2, 2, 0, -2, -4, -4, -2, 0, 2, 2, 1, 1]),
    ([2, 2, 2], [1, -1, 2, -4, 4, -6, 6, -6, 6, -4, 4, -2, 1, -1]),
    ([3, 1, 1, 1], [1, 2, 2, 4, 2, -5, -6, -6, -5, 2, 4, 2, 2, 1]),
    ([3, 2, 1], [1, 0, 0, 0, -2, -1, 0, 0, 1, 2, 0, 0, 0, -1]),
    ([3, 3], [1, -1, -1, 1, -1, 1, 0, 0, 1, -1, 1, -1, -1, 1]),
    ([4, 1, 1], [1, 1, 0, 0, -2, -4, -2, 2, 4, 2, 0, 0, -1, -1]),
    ([4, 2], [1, -1, 0, -2, 2, 0, 0, 0, 0, 2, -2, 0, -1, 1]),
    ([5, 1], [1, 0, -1, -1, -2, 1, 2, 2, 1, -2, -1, -1, 0, 1]),
    ([6], [1, -1, -1, -1, 1, 3, 0, 0, -3, -1, 1, 1, 1, -1]),
    ])

hess_left[(0, 0, 0, 0, 0, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 5, 14, 34, 68, 106, 132, 132, 106, 68, 34, 14, 5, 1]),
    ([2, 1, 1, 1, 1], [1, 5, 14, 32, 60, 90, 110, 110, 90, 60, 32, 14, 5, 1]),
    ([2, 2, 1, 1], [1, 5, 14, 30, 52, 74, 88, 88, 74, 52, 30, 14, 5, 1]),
    ([2, 2, 2], [1, 5, 14, 28, 44, 58, 66, 66, 58, 44, 28, 14, 5, 1]),
    ([3, 1, 1, 1], [1, 5, 14, 31, 56, 82, 99, 99, 82, 56, 31, 14, 5, 1]),
    ([3, 2, 1], [1, 5, 14, 29, 48, 66, 77, 77, 66, 48, 29, 14, 5, 1]),
    ([3, 3], [1, 5, 14, 28, 44, 58, 66, 66, 58, 44, 28, 14, 5, 1]),
    ([4, 1, 1], [1, 5, 14, 30, 52, 74, 88, 88, 74, 52, 30, 14, 5, 1]),
    ([4, 2], [1, 5, 14, 28, 44, 58, 66, 66, 58, 44, 28, 14, 5, 1]),
    ([5, 1], [1, 5, 14, 29, 48, 66, 77, 77, 66, 48, 29, 14, 5, 1]),
    ([6], [1, 5, 14, 28, 44, 58, 66, 66, 58, 44, 28, 14, 5, 1]),
    ])

hess_right[(0, 0, 0, 0, 0, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 5, 19, 48, 89, 127, 142, 127, 89, 48, 19, 5, 1]),
    ([2, 1, 1, 1, 1], [1, 3, 9, 14, 17, 11, 0, -11, -17, -14, -9, -3, -1]),
    ([2, 2, 1, 1], [1, 1, 3, 0, 1, -5, -2, -5, 1, 0, 3, 1, 1]),
    ([2, 2, 2], [1, -1, 1, -2, 1, -1, 0, 1, -1, 2, -1, 1, -1]),
    ([3, 1, 1, 1], [1, 2, 4, 3, -1, -5, -8, -5, -1, 3, 4, 2, 1]),
    ([3, 2, 1], [1, 0, 0, -1, -1, -1, 0, 1, 1, 1, 0, 0, -1]),
    ([3, 3], [1, -1, -2, 3, -1, -2, 4, -2, -1, 3, -2, -1, 1]),
    ([4, 1, 1], [1, 1, 1, -2, -3, -3, 0, 3, 3, 2, -1, -1, -1]),
    ([4, 2], [1, -1, -1, 0, 1, 1, -2, 1, 1, 0, -1, -1, 1]),
    ([5, 1], [1, 0, -1, -2, -1, 2, 2, 2, -1, -2, -1, 0, 1]),
    ([6], [1, -1, -2, 1, 1, 2, 0, -2, -1, -1, 2, 1, -1]),
    ])

hess_left[(0, 0, 0, 0, 0, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 5, 19, 48, 89, 127, 142, 127, 89, 48, 19, 5, 1]),
    ([2, 1, 1, 1, 1], [1, 5, 17, 40, 71, 99, 110, 99, 71, 40, 17, 5, 1]),
    ([2, 2, 1, 1], [1, 5, 15, 32, 53, 71, 78, 71, 53, 32, 15, 5, 1]),
    ([2, 2, 2], [1, 5, 13, 24, 35, 43, 46, 43, 35, 24, 13, 5, 1]),
    ([3, 1, 1, 1], [1, 5, 16, 36, 62, 85, 94, 85, 62, 36, 16, 5, 1]),
    ([3, 2, 1], [1, 5, 14, 28, 44, 57, 62, 57, 44, 28, 14, 5, 1]),
    ([3, 3], [1, 5, 13, 24, 35, 43, 46, 43, 35, 24, 13, 5, 1]),
    ([4, 1, 1], [1, 5, 15, 32, 53, 71, 78, 71, 53, 32, 15, 5, 1]),
    ([4, 2], [1, 5, 13, 24, 35, 43, 46, 43, 35, 24, 13, 5, 1]),
    ([5, 1], [1, 5, 14, 28, 44, 57, 62, 57, 44, 28, 14, 5, 1]),
    ([6], [1, 5, 13, 24, 35, 43, 46, 43, 35, 24, 13, 5, 1]),
    ])

hess_right[(0, 0, 0, 0, 0, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 10, 33, 69, 110, 137, 137, 110, 69, 33, 10, 1]),
    ([2, 1, 1, 1, 1], [1, 6, 11, 15, 14, 5, -5, -14, -15, -11, -6, -1]),
    ([2, 2, 1, 1], [1, 2, 1, 1, -2, -3, -3, -2, 1, 1, 2, 1]),
    ([2, 2, 2], [1, -2, 3, -5, 6, -7, 7, -6, 5, -3, 2, -1]),
    ([3, 1, 1, 1], [1, 4, 3, 0, -1, -7, -7, -1, 0, 3, 4, 1]),
    ([3, 2, 1], [1, 0, -1, 0, -1, -1, 1, 1, 0, 1, 0, -1]),
    ([3, 3], [1, -2, 0, 3, -4, 2, 2, -4, 3, 0, -2, 1]),
    ([4, 1, 1], [1, 2, -1, -3, -2, -1, 1, 2, 3, 1, -2, -1]),
    ([4, 2], [1, -2, 1, -1, 2, -1, -1, 2, -1, 1, -2, 1]),
    ([5, 1], [1, 0, -2, -1, 0, 2, 2, 0, -1, -2, 0, 1]),
    ([6], [1, -2, 0, 1, 0, 2, -2, 0, -1, 0, 2, -1]),
    ])

hess_left[(0, 0, 0, 0, 0, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 10, 33, 69, 110, 137, 137, 110, 69, 33, 10, 1]),
    ([2, 1, 1, 1, 1], [1, 8, 25, 51, 80, 99, 99, 80, 51, 25, 8, 1]),
    ([2, 2, 1, 1], [1, 6, 17, 33, 50, 61, 61, 50, 33, 17, 6, 1]),
    ([2, 2, 2], [1, 4, 9, 15, 20, 23, 23, 20, 15, 9, 4, 1]),
    ([3, 1, 1, 1], [1, 7, 21, 42, 65, 80, 80, 65, 42, 21, 7, 1]),
    ([3, 2, 1], [1, 5, 13, 24, 35, 42, 42, 35, 24, 13, 5, 1]),
    ([3, 3], [1, 4, 9, 15, 20, 23, 23, 20, 15, 9, 4, 1]),
    ([4, 1, 1], [1, 6, 17, 33, 50, 61, 61, 50, 33, 17, 6, 1]),
    ([4, 2], [1, 4, 9, 15, 20, 23, 23, 20, 15, 9, 4, 1]),
    ([5, 1], [1, 5, 13, 24, 35, 42, 42, 35, 24, 13, 5, 1]),
    ([6], [1, 4, 9, 15, 20, 23, 23, 20, 15, 9, 4, 1]),
    ])

hess_right[(0, 0, 0, 0, 0, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 24, 54, 90, 120, 132, 120, 90, 54, 24, 6]),
    ([2, 1, 1, 1, 1], [4, 8, 12, 12, 8, 0, -8, -12, -12, -8, -4]),
    ([2, 2, 1, 1], [2, 0, 2, -2, 0, -4, 0, -2, 2, 0, 2]),
    ([3, 1, 1, 1], [3, 3, 0, 0, -3, -6, -3, 0, 0, 3, 3]),
    ([3, 2, 1], [1, -1, 0, 0, -1, 0, 1, 0, 0, 1, -1]),
    ([4, 1, 1], [2, 0, -2, -2, 0, 0, 0, 2, 2, 0, -2]),
    ([5, 1], [1, -1, -1, 0, 0, 2, 0, 0, -1, -1, 1]),
    ])

hess_left[(0, 0, 0, 0, 0, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 24, 54, 90, 120, 132, 120, 90, 54, 24, 6]),
    ([2, 1, 1, 1, 1], [4, 16, 36, 60, 80, 88, 80, 60, 36, 16, 4]),
    ([2, 2, 1, 1], [2, 8, 18, 30, 40, 44, 40, 30, 18, 8, 2]),
    ([3, 1, 1, 1], [3, 12, 27, 45, 60, 66, 60, 45, 27, 12, 3]),
    ([3, 2, 1], [1, 4, 9, 15, 20, 22, 20, 15, 9, 4, 1]),
    ([4, 1, 1], [2, 8, 18, 30, 40, 44, 40, 30, 18, 8, 2]),
    ([5, 1], [1, 4, 9, 15, 20, 22, 20, 15, 9, 4, 1]),
    ])

hess_right[(0, 0, 0, 0, 1, 1)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 5, 14, 34, 68, 106, 132, 132, 106, 68, 34, 14, 5, 1]),
    ([2, 1, 1, 1, 1], [1, 3, 6, 12, 16, 14, 6, -6, -14, -16, -12, -6, -3, -1]),
    ([2, 2, 1, 1], [1, 1, 2, 2, 0, -2, -4, -4, -2, 0, 2, 2, 1, 1]),
    ([2, 2, 2], [1, -1, 2, -4, 4, -6, 6, -6, 6, -4, 4, -2, 1, -1]),
    ([3, 1, 1, 1], [1, 2, 2, 4, 2, -5, -6, -6, -5, 2, 4, 2, 2, 1]),
    ([3, 2, 1], [1, 0, 0, 0, -2, -1, 0, 0, 1, 2, 0, 0, 0, -1]),
    ([3, 3], [1, -1, -1, 1, -1, 1, 0, 0, 1, -1, 1, -1, -1, 1]),
    ([4, 1, 1], [1, 1, 0, 0, -2, -4, -2, 2, 4, 2, 0, 0, -1, -1]),
    ([4, 2], [1, -1, 0, -2, 2, 0, 0, 0, 0, 2, -2, 0, -1, 1]),
    ([5, 1], [1, 0, -1, -1, -2, 1, 2, 2, 1, -2, -1, -1, 0, 1]),
    ([6], [1, -1, -1, -1, 1, 3, 0, 0, -3, -1, 1, 1, 1, -1]),
    ])

hess_left[(0, 0, 0, 0, 1, 1)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 5, 14, 34, 68, 106, 132, 132, 106, 68, 34, 14, 5, 1]),
    ([2, 1, 1, 1, 1], [1, 5, 14, 32, 60, 90, 110, 110, 90, 60, 32, 14, 5, 1]),
    ([2, 2, 1, 1], [1, 5, 14, 30, 52, 74, 88, 88, 74, 52, 30, 14, 5, 1]),
    ([2, 2, 2], [1, 5, 14, 28, 44, 58, 66, 66, 58, 44, 28, 14, 5, 1]),
    ([3, 1, 1, 1], [1, 5, 14, 31, 56, 82, 99, 99, 82, 56, 31, 14, 5, 1]),
    ([3, 2, 1], [1, 5, 14, 29, 48, 66, 77, 77, 66, 48, 29, 14, 5, 1]),
    ([3, 3], [1, 5, 14, 28, 44, 58, 66, 66, 58, 44, 28, 14, 5, 1]),
    ([4, 1, 1], [1, 5, 14, 30, 52, 74, 88, 88, 74, 52, 30, 14, 5, 1]),
    ([4, 2], [1, 5, 14, 28, 44, 58, 66, 66, 58, 44, 28, 14, 5, 1]),
    ([5, 1], [1, 5, 14, 29, 48, 66, 77, 77, 66, 48, 29, 14, 5, 1]),
    ([6], [1, 5, 14, 28, 44, 58, 66, 66, 58, 44, 28, 14, 5, 1]),
    ])

hess_right[(0, 0, 0, 0, 1, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 5, 14, 39, 82, 136, 166, 136, 82, 39, 14, 5, 1]),
    ([2, 1, 1, 1, 1], [1, 3, 6, 15, 18, 18, 0, -18, -18, -15, -6, -3, -1]),
    ([2, 2, 1, 1], [1, 1, 2, 3, -2, 0, -10, 0, -2, 3, 2, 1, 1]),
    ([2, 2, 2], [1, -1, 2, -5, 6, -6, 0, 6, -6, 5, -2, 1, -1]),
    ([3, 1, 1, 1], [1, 2, 2, 6, 1, -8, -8, -8, 1, 6, 2, 2, 1]),
    ([3, 2, 1], [1, 0, 0, 0, -3, 0, 0, 0, 3, 0, 0, 0, -1]),
    ([3, 3], [1, -1, -1, 0, 1, 1, -2, 1, 1, 0, -1, -1, 1]),
    ([4, 1, 1], [1, 1, 0, 1, -4, -6, 0, 6, 4, -1, 0, -1, -1]),
    ([4, 2], [1, -1, 0, -3, 4, 0, -2, 0, 4, -3, 0, -1, 1]),
    ([5, 1], [1, 0, -1, -1, -3, 1, 6, 1, -3, -1, -1, 0, 1]),
    ([6], [1, -1, -1, -2, 3, 3, 0, -3, -3, 2, 1, 1, -1]),
    ])

hess_left[(0, 0, 0, 0, 1, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 5, 14, 39, 82, 136, 166, 136, 82, 39, 14, 5, 1]),
    ([2, 1, 1, 1, 1], [1, 5, 14, 35, 68, 104, 122, 104, 68, 35, 14, 5, 1]),
    ([2, 2, 1, 1], [1, 5, 14, 31, 54, 76, 86, 76, 54, 31, 14, 5, 1]),
    ([2, 2, 2], [1, 5, 14, 27, 40, 52, 58, 52, 40, 27, 14, 5, 1]),
    ([3, 1, 1, 1], [1, 5, 14, 33, 61, 88, 100, 88, 61, 33, 14, 5, 1]),
    ([3, 2, 1], [1, 5, 14, 29, 47, 62, 68, 62, 47, 29, 14, 5, 1]),
    ([3, 3], [1, 5, 14, 27, 40, 49, 52, 49, 40, 27, 14, 5, 1]),
    ([4, 1, 1], [1, 5, 14, 31, 54, 74, 82, 74, 54, 31, 14, 5, 1]),
    ([4, 2], [1, 5, 14, 27, 40, 50, 54, 50, 40, 27, 14, 5, 1]),
    ([5, 1], [1, 5, 14, 29, 47, 61, 66, 61, 47, 29, 14, 5, 1]),
    ([6], [1, 5, 14, 27, 40, 49, 52, 49, 40, 27, 14, 5, 1]),
    ])

hess_right[(0, 0, 0, 0, 1, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 5, 19, 53, 112, 170, 170, 112, 53, 19, 5, 1]),
    ([2, 1, 1, 1, 1], [1, 3, 9, 17, 22, 12, -12, -22, -17, -9, -3, -1]),
    ([2, 2, 1, 1], [1, 1, 3, 1, 0, -6, -6, 0, 1, 3, 1, 1]),
    ([2, 2, 2], [1, -1, 1, -3, 6, -12, 12, -6, 3, -1, 1, -1]),
    ([3, 1, 1, 1], [1, 2, 4, 5, -2, -10, -10, -2, 5, 4, 2, 1]),
    ([3, 2, 1], [1, 0, 0, -1, -2, 0, 0, 2, 1, 0, 0, -1]),
    ([3, 3], [1, -1, -2, 2, 1, -1, -1, 1, 2, -2, -1, 1]),
    ([4, 1, 1], [1, 1, 1, -1, -6, -4, 4, 6, 1, -1, -1, -1]),
    ([4, 2], [1, -1, -1, -1, 4, -2, -2, 4, -1, -1, -1, 1]),
    ([5, 1], [1, 0, -1, -2, -3, 5, 5, -3, -2, -1, 0, 1]),
    ([6], [1, -1, -2, 0, 3, 3, -3, -3, 0, 2, 1, -1]),
    ])

hess_left[(0, 0, 0, 0, 1, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 5, 19, 53, 112, 170, 170, 112, 53, 19, 5, 1]),
    ([2, 1, 1, 1, 1], [1, 5, 17, 43, 82, 116, 116, 82, 43, 17, 5, 1]),
    ([2, 2, 1, 1], [1, 5, 15, 33, 56, 74, 74, 56, 33, 15, 5, 1]),
    ([2, 2, 2], [1, 5, 13, 23, 34, 44, 44, 34, 23, 13, 5, 1]),
    ([3, 1, 1, 1], [1, 5, 16, 38, 67, 89, 89, 67, 38, 16, 5, 1]),
    ([3, 2, 1], [1, 5, 14, 28, 43, 53, 53, 43, 28, 14, 5, 1]),
    ([3, 3], [1, 5, 13, 23, 31, 35, 35, 31, 23, 13, 5, 1]),
    ([4, 1, 1], [1, 5, 15, 33, 54, 68, 68, 54, 33, 15, 5, 1]),
    ([4, 2], [1, 5, 13, 23, 32, 38, 38, 32, 23, 13, 5, 1]),
    ([5, 1], [1, 5, 14, 28, 42, 50, 50, 42, 28, 14, 5, 1]),
    ([6], [1, 5, 13, 23, 31, 35, 35, 31, 23, 13, 5, 1]),
    ])

hess_right[(0, 0, 0, 0, 1, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 10, 33, 83, 146, 174, 146, 83, 33, 10, 1]),
    ([2, 1, 1, 1, 1], [1, 6, 11, 21, 16, 0, -16, -21, -11, -6, -1]),
    ([2, 2, 1, 1], [1, 2, 1, 3, -6, -2, -6, 3, 1, 2, 1]),
    ([2, 2, 2], [1, -2, 3, -3, 0, 0, 0, 3, -3, 2, -1]),
    ([3, 1, 1, 1], [1, 4, 3, 2, -4, -12, -4, 2, 3, 4, 1]),
    ([3, 2, 1], [1, 0, -1, 0, -2, 0, 2, 0, 1, 0, -1]),
    ([3, 3], [1, -2, 0, 2, -1, 0, -1, 2, 0, -2, 1]),
    ([4, 1, 1], [1, 2, -1, -3, -4, 0, 4, 3, 1, -2, -1]),
    ([4, 2], [1, -2, 1, -1, 2, -2, 2, -1, 1, -2, 1]),
    ([5, 1], [1, 0, -2, -2, 1, 4, 1, -2, -2, 0, 1]),
    ([6], [1, -2, 0, 0, 3, 0, -3, 0, 0, 2, -1]),
    ])

hess_left[(0, 0, 0, 0, 1, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 10, 33, 83, 146, 174, 146, 83, 33, 10, 1]),
    ([2, 1, 1, 1, 1], [1, 8, 25, 57, 94, 110, 94, 57, 25, 8, 1]),
    ([2, 2, 1, 1], [1, 6, 17, 35, 54, 62, 54, 35, 17, 6, 1]),
    ([2, 2, 2], [1, 4, 9, 17, 26, 30, 26, 17, 9, 4, 1]),
    ([3, 1, 1, 1], [1, 7, 21, 44, 68, 78, 68, 44, 21, 7, 1]),
    ([3, 2, 1], [1, 5, 13, 24, 34, 38, 34, 24, 13, 5, 1]),
    ([3, 3], [1, 4, 9, 14, 17, 18, 17, 14, 9, 4, 1]),
    ([4, 1, 1], [1, 6, 17, 33, 48, 54, 48, 33, 17, 6, 1]),
    ([4, 2], [1, 4, 9, 15, 20, 22, 20, 15, 9, 4, 1]),
    ([5, 1], [1, 5, 13, 23, 31, 34, 31, 23, 13, 5, 1]),
    ([6], [1, 4, 9, 14, 17, 18, 17, 14, 9, 4, 1]),
    ])

hess_right[(0, 0, 0, 0, 1, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 24, 54, 114, 162, 162, 114, 54, 24, 6]),
    ([2, 1, 1, 1, 1], [4, 8, 12, 20, 4, -4, -20, -12, -8, -4]),
    ([2, 2, 1, 1], [2, 0, 2, -2, -2, -2, -2, 2, 0, 2]),
    ([3, 1, 1, 1], [3, 3, 0, 3, -9, -9, 3, 0, 3, 3]),
    ([3, 2, 1], [1, -1, 0, -1, 1, -1, 1, 0, 1, -1]),
    ([4, 1, 1], [2, 0, -2, -2, -2, 2, 2, 2, 0, -2]),
    ([5, 1], [1, -1, -1, -1, 2, 2, -1, -1, -1, 1]),
    ])

hess_left[(0, 0, 0, 0, 1, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 24, 54, 114, 162, 162, 114, 54, 24, 6]),
    ([2, 1, 1, 1, 1], [4, 16, 36, 68, 92, 92, 68, 36, 16, 4]),
    ([2, 2, 1, 1], [2, 8, 18, 30, 38, 38, 30, 18, 8, 2]),
    ([3, 1, 1, 1], [3, 12, 27, 48, 63, 63, 48, 27, 12, 3]),
    ([3, 2, 1], [1, 4, 9, 14, 17, 17, 14, 9, 4, 1]),
    ([4, 1, 1], [2, 8, 18, 30, 38, 38, 30, 18, 8, 2]),
    ([5, 1], [1, 4, 9, 14, 17, 17, 14, 9, 4, 1]),
    ])

hess_right[(0, 0, 0, 0, 2, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 5, 19, 53, 112, 170, 170, 112, 53, 19, 5, 1]),
    ([2, 1, 1, 1, 1], [1, 3, 9, 17, 22, 12, -12, -22, -17, -9, -3, -1]),
    ([2, 2, 1, 1], [1, 1, 3, 1, 0, -6, -6, 0, 1, 3, 1, 1]),
    ([2, 2, 2], [1, -1, 1, -3, 6, -12, 12, -6, 3, -1, 1, -1]),
    ([3, 1, 1, 1], [1, 2, 4, 5, -2, -10, -10, -2, 5, 4, 2, 1]),
    ([3, 2, 1], [1, 0, 0, -1, -2, 0, 0, 2, 1, 0, 0, -1]),
    ([3, 3], [1, -1, -2, 2, 1, -1, -1, 1, 2, -2, -1, 1]),
    ([4, 1, 1], [1, 1, 1, -1, -6, -4, 4, 6, 1, -1, -1, -1]),
    ([4, 2], [1, -1, -1, -1, 4, -2, -2, 4, -1, -1, -1, 1]),
    ([5, 1], [1, 0, -1, -2, -3, 5, 5, -3, -2, -1, 0, 1]),
    ([6], [1, -1, -2, 0, 3, 3, -3, -3, 0, 2, 1, -1]),
    ])

hess_left[(0, 0, 0, 0, 2, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 5, 19, 53, 112, 170, 170, 112, 53, 19, 5, 1]),
    ([2, 1, 1, 1, 1], [1, 5, 17, 43, 82, 116, 116, 82, 43, 17, 5, 1]),
    ([2, 2, 1, 1], [1, 5, 15, 33, 56, 74, 74, 56, 33, 15, 5, 1]),
    ([2, 2, 2], [1, 5, 13, 23, 34, 44, 44, 34, 23, 13, 5, 1]),
    ([3, 1, 1, 1], [1, 5, 16, 38, 67, 89, 89, 67, 38, 16, 5, 1]),
    ([3, 2, 1], [1, 5, 14, 28, 43, 53, 53, 43, 28, 14, 5, 1]),
    ([3, 3], [1, 5, 13, 23, 31, 35, 35, 31, 23, 13, 5, 1]),
    ([4, 1, 1], [1, 5, 15, 33, 54, 68, 68, 54, 33, 15, 5, 1]),
    ([4, 2], [1, 5, 13, 23, 32, 38, 38, 32, 23, 13, 5, 1]),
    ([5, 1], [1, 5, 14, 28, 42, 50, 50, 42, 28, 14, 5, 1]),
    ([6], [1, 5, 13, 23, 31, 35, 35, 31, 23, 13, 5, 1]),
    ])

hess_right[(0, 0, 0, 0, 2, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 5, 24, 76, 155, 198, 155, 76, 24, 5, 1]),
    ([2, 1, 1, 1, 1], [1, 3, 12, 22, 23, 0, -23, -22, -12, -3, -1]),
    ([2, 2, 1, 1], [1, 1, 4, 0, -1, -10, -1, 0, 4, 1, 1]),
    ([2, 2, 2], [1, -1, 0, 2, -5, 0, 5, -2, 0, 1, -1]),
    ([3, 1, 1, 1], [1, 2, 6, 4, -7, -12, -7, 4, 6, 2, 1]),
    ([3, 2, 1], [1, 0, 0, -2, -1, 0, 1, 2, 0, 0, -1]),
    ([3, 3], [1, -1, -3, 4, 2, -6, 2, 4, -3, -1, 1]),
    ([4, 1, 1], [1, 1, 2, -4, -7, 0, 7, 4, -2, -1, -1]),
    ([4, 2], [1, -1, -2, 2, 1, -2, 1, 2, -2, -1, 1]),
    ([5, 1], [1, 0, -1, -4, 0, 8, 0, -4, -1, 0, 1]),
    ([6], [1, -1, -3, 2, 4, 0, -4, -2, 3, 1, -1]),
    ])

hess_left[(0, 0, 0, 0, 2, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 5, 24, 76, 155, 198, 155, 76, 24, 5, 1]),
    ([2, 1, 1, 1, 1], [1, 5, 20, 54, 99, 122, 99, 54, 20, 5, 1]),
    ([2, 2, 1, 1], [1, 5, 16, 36, 59, 70, 59, 36, 16, 5, 1]),
    ([2, 2, 2], [1, 5, 12, 22, 35, 42, 35, 22, 12, 5, 1]),
    ([3, 1, 1, 1], [1, 5, 18, 43, 71, 84, 71, 43, 18, 5, 1]),
    ([3, 2, 1], [1, 5, 14, 27, 39, 44, 39, 27, 14, 5, 1]),
    ([3, 3], [1, 5, 12, 19, 23, 24, 23, 19, 12, 5, 1]),
    ([4, 1, 1], [1, 5, 16, 34, 51, 58, 51, 34, 16, 5, 1]),
    ([4, 2], [1, 5, 12, 20, 27, 30, 27, 20, 12, 5, 1]),
    ([5, 1], [1, 5, 14, 26, 35, 38, 35, 26, 14, 5, 1]),
    ([6], [1, 5, 12, 19, 23, 24, 23, 19, 12, 5, 1]),
    ])

hess_right[(0, 0, 0, 0, 2, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 10, 47, 119, 183, 183, 119, 47, 10, 1]),
    ([2, 1, 1, 1, 1], [1, 6, 17, 23, 11, -11, -23, -17, -6, -1]),
    ([2, 2, 1, 1], [1, 2, 3, -1, -5, -5, -1, 3, 2, 1]),
    ([2, 2, 2], [1, -2, 5, -9, 7, -7, 9, -5, 2, -1]),
    ([3, 1, 1, 1], [1, 4, 5, -1, -9, -9, -1, 5, 4, 1]),
    ([3, 2, 1], [1, 0, -1, -1, -1, 1, 1, 1, 0, -1]),
    ([3, 3], [1, -2, -1, 5, -3, -3, 5, -1, -2, 1]),
    ([4, 1, 1], [1, 2, -1, -5, -3, 3, 5, 1, -2, -1]),
    ([4, 2], [1, -2, 1, -1, 1, 1, -1, 1, -2, 1]),
    ([5, 1], [1, 0, -3, -1, 3, 3, -1, -3, 0, 1]),
    ([6], [1, -2, -1, 3, 1, -1, -3, 1, 2, -1]),
    ])

hess_left[(0, 0, 0, 0, 2, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 10, 47, 119, 183, 183, 119, 47, 10, 1]),
    ([2, 1, 1, 1, 1], [1, 8, 31, 71, 105, 105, 71, 31, 8, 1]),
    ([2, 2, 1, 1], [1, 6, 19, 39, 55, 55, 39, 19, 6, 1]),
    ([2, 2, 2], [1, 4, 11, 23, 33, 33, 23, 11, 4, 1]),
    ([3, 1, 1, 1], [1, 7, 23, 47, 66, 66, 47, 23, 7, 1]),
    ([3, 2, 1], [1, 5, 13, 23, 30, 30, 23, 13, 5, 1]),
    ([3, 3], [1, 4, 8, 11, 12, 12, 11, 8, 4, 1]),
    ([4, 1, 1], [1, 6, 17, 31, 41, 41, 31, 17, 6, 1]),
    ([4, 2], [1, 4, 9, 15, 19, 19, 15, 9, 4, 1]),
    ([5, 1], [1, 5, 12, 19, 23, 23, 19, 12, 5, 1]),
    ([6], [1, 4, 8, 11, 12, 12, 11, 8, 4, 1]),
    ])

hess_right[(0, 0, 0, 0, 2, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 24, 78, 156, 192, 156, 78, 24, 6]),
    ([2, 1, 1, 1, 1], [4, 8, 20, 16, 0, -16, -20, -8, -4]),
    ([2, 2, 1, 1], [2, 0, 2, -4, 0, -4, 2, 0, 2]),
    ([3, 1, 1, 1], [3, 3, 3, -3, -12, -3, 3, 3, 3]),
    ([3, 2, 1], [1, -1, -1, 1, 0, -1, 1, 1, -1]),
    ([4, 1, 1], [2, 0, -2, -4, 0, 4, 2, 0, -2]),
    ([5, 1], [1, -1, -2, 1, 2, 1, -2, -1, 1]),
    ])

hess_left[(0, 0, 0, 0, 2, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 24, 78, 156, 192, 156, 78, 24, 6]),
    ([2, 1, 1, 1, 1], [4, 16, 44, 80, 96, 80, 44, 16, 4]),
    ([2, 2, 1, 1], [2, 8, 18, 28, 32, 28, 18, 8, 2]),
    ([3, 1, 1, 1], [3, 12, 30, 51, 60, 51, 30, 12, 3]),
    ([3, 2, 1], [1, 4, 8, 11, 12, 11, 8, 4, 1]),
    ([4, 1, 1], [2, 8, 18, 28, 32, 28, 18, 8, 2]),
    ([5, 1], [1, 4, 8, 11, 12, 11, 8, 4, 1]),
    ])

hess_right[(0, 0, 0, 0, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 10, 47, 119, 183, 183, 119, 47, 10, 1]),
    ([2, 1, 1, 1, 1], [1, 6, 17, 23, 11, -11, -23, -17, -6, -1]),
    ([2, 2, 1, 1], [1, 2, 3, -1, -5, -5, -1, 3, 2, 1]),
    ([2, 2, 2], [1, -2, 5, -9, 7, -7, 9, -5, 2, -1]),
    ([3, 1, 1, 1], [1, 4, 5, -1, -9, -9, -1, 5, 4, 1]),
    ([3, 2, 1], [1, 0, -1, -1, -1, 1, 1, 1, 0, -1]),
    ([3, 3], [1, -2, -1, 5, -3, -3, 5, -1, -2, 1]),
    ([4, 1, 1], [1, 2, -1, -5, -3, 3, 5, 1, -2, -1]),
    ([4, 2], [1, -2, 1, -1, 1, 1, -1, 1, -2, 1]),
    ([5, 1], [1, 0, -3, -1, 3, 3, -1, -3, 0, 1]),
    ([6], [1, -2, -1, 3, 1, -1, -3, 1, 2, -1]),
    ])

hess_left[(0, 0, 0, 0, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 10, 47, 119, 183, 183, 119, 47, 10, 1]),
    ([2, 1, 1, 1, 1], [1, 8, 31, 71, 105, 105, 71, 31, 8, 1]),
    ([2, 2, 1, 1], [1, 6, 19, 39, 55, 55, 39, 19, 6, 1]),
    ([2, 2, 2], [1, 4, 11, 23, 33, 33, 23, 11, 4, 1]),
    ([3, 1, 1, 1], [1, 7, 23, 47, 66, 66, 47, 23, 7, 1]),
    ([3, 2, 1], [1, 5, 13, 23, 30, 30, 23, 13, 5, 1]),
    ([3, 3], [1, 4, 8, 11, 12, 12, 11, 8, 4, 1]),
    ([4, 1, 1], [1, 6, 17, 31, 41, 41, 31, 17, 6, 1]),
    ([4, 2], [1, 4, 9, 15, 19, 19, 15, 9, 4, 1]),
    ([5, 1], [1, 5, 12, 19, 23, 23, 19, 12, 5, 1]),
    ([6], [1, 4, 8, 11, 12, 12, 11, 8, 4, 1]),
    ])

hess_right[(0, 0, 0, 0, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 24, 83, 156, 192, 156, 83, 24, 1]),
    ([2, 1, 1, 1, 1], [1, 12, 19, 18, 0, -18, -19, -12, -1]),
    ([2, 2, 1, 1], [1, 4, -1, 0, -8, 0, -1, 4, 1]),
    ([2, 2, 2], [1, 0, -1, -2, 0, 2, 1, 0, -1]),
    ([3, 1, 1, 1], [1, 6, 2, -6, -6, -6, 2, 6, 1]),
    ([3, 2, 1], [1, 0, -2, 0, 0, 0, 2, 0, -1]),
    ([3, 3], [1, -3, 2, 3, -6, 3, 2, -3, 1]),
    ([4, 1, 1], [1, 2, -3, -4, 0, 4, 3, -2, -1]),
    ([4, 2], [1, -2, 1, -2, 4, -2, 1, -2, 1]),
    ([5, 1], [1, -1, -2, 1, 2, 1, -2, -1, 1]),
    ([6], [1, -3, 2, 1, 0, -1, -2, 3, -1]),
    ])

hess_left[(0, 0, 0, 0, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 24, 83, 156, 192, 156, 83, 24, 1]),
    ([2, 1, 1, 1, 1], [1, 14, 45, 82, 100, 82, 45, 14, 1]),
    ([2, 2, 1, 1], [1, 8, 23, 40, 48, 40, 23, 8, 1]),
    ([2, 2, 2], [1, 6, 17, 30, 36, 30, 17, 6, 1]),
    ([3, 1, 1, 1], [1, 9, 26, 45, 54, 45, 26, 9, 1]),
    ([3, 2, 1], [1, 5, 12, 19, 22, 19, 12, 5, 1]),
    ([3, 3], [1, 3, 5, 6, 6, 6, 5, 3, 1]),
    ([4, 1, 1], [1, 6, 15, 24, 28, 24, 15, 6, 1]),
    ([4, 2], [1, 4, 9, 14, 16, 14, 9, 4, 1]),
    ([5, 1], [1, 4, 8, 11, 12, 11, 8, 4, 1]),
    ([6], [1, 3, 5, 6, 6, 6, 5, 3, 1]),
    ])

hess_right[(0, 0, 0, 0, 3, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 48, 120, 186, 186, 120, 48, 6]),
    ([2, 1, 1, 1, 1], [4, 16, 16, 12, -12, -16, -16, -4]),
    ([2, 2, 1, 1], [2, 0, 0, -2, -2, 0, 0, 2]),
    ([3, 1, 1, 1], [3, 6, -3, -6, -6, -3, 6, 3]),
    ([3, 2, 1], [1, -2, 1, 0, 0, -1, 2, -1]),
    ([4, 1, 1], [2, 0, -4, -2, 2, 4, 0, -2]),
    ([5, 1], [1, -2, 0, 1, 1, 0, -2, 1]),
    ])

hess_left[(0, 0, 0, 0, 3, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 48, 120, 186, 186, 120, 48, 6]),
    ([2, 1, 1, 1, 1], [4, 24, 56, 84, 84, 56, 24, 4]),
    ([2, 2, 1, 1], [2, 8, 16, 22, 22, 16, 8, 2]),
    ([3, 1, 1, 1], [3, 15, 33, 48, 48, 33, 15, 3]),
    ([3, 2, 1], [1, 3, 5, 6, 6, 5, 3, 1]),
    ([4, 1, 1], [2, 8, 16, 22, 22, 16, 8, 2]),
    ([5, 1], [1, 3, 5, 6, 6, 5, 3, 1]),
    ])

hess_right[(0, 0, 0, 0, 4, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [15, 60, 120, 165, 165, 120, 60, 15]),
    ([2, 1, 1, 1, 1], [7, 14, 14, 7, -7, -14, -14, -7]),
    ([2, 2, 1, 1], [3, 0, 0, -3, -3, 0, 0, 3]),
    ([2, 2, 2], [3, -6, 6, -9, 9, -6, 6, -3]),
    ([3, 1, 1, 1], [3, 3, -3, -3, -3, -3, 3, 3]),
    ([3, 2, 1], [1, -1, -1, 1, -1, 1, 1, -1]),
    ([4, 1, 1], [1, 0, -2, -1, 1, 2, 0, -1]),
    ([4, 2], [1, -2, 0, 1, 1, 0, -2, 1]),
    ])

hess_left[(0, 0, 0, 0, 4, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [15, 60, 120, 165, 165, 120, 60, 15]),
    ([2, 1, 1, 1, 1], [7, 28, 56, 77, 77, 56, 28, 7]),
    ([2, 2, 1, 1], [3, 12, 24, 33, 33, 24, 12, 3]),
    ([2, 2, 2], [3, 12, 24, 33, 33, 24, 12, 3]),
    ([3, 1, 1, 1], [3, 12, 24, 33, 33, 24, 12, 3]),
    ([3, 2, 1], [1, 4, 8, 11, 11, 8, 4, 1]),
    ([4, 1, 1], [1, 4, 8, 11, 11, 8, 4, 1]),
    ([4, 2], [1, 4, 8, 11, 11, 8, 4, 1]),
    ])

hess_right[(0, 0, 0, 0, 4, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [30, 90, 150, 180, 150, 90, 30]),
    ([2, 1, 1, 1, 1], [12, 12, 12, 0, -12, -12, -12]),
    ([2, 2, 1, 1], [2, -2, 2, -4, 2, -2, 2]),
    ([3, 1, 1, 1], [6, 0, -6, 0, -6, 0, 6]),
    ([4, 1, 1], [2, -2, -2, 0, 2, 2, -2]),
    ])

hess_left[(0, 0, 0, 0, 4, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [30, 90, 150, 180, 150, 90, 30]),
    ([2, 1, 1, 1, 1], [12, 36, 60, 72, 60, 36, 12]),
    ([2, 2, 1, 1], [2, 6, 10, 12, 10, 6, 2]),
    ([3, 1, 1, 1], [6, 18, 30, 36, 30, 18, 6]),
    ([4, 1, 1], [2, 6, 10, 12, 10, 6, 2]),
    ])

hess_right[(0, 0, 0, 1, 1, 1)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 5, 19, 48, 89, 127, 142, 127, 89, 48, 19, 5, 1]),
    ([2, 1, 1, 1, 1], [1, 3, 9, 14, 17, 11, 0, -11, -17, -14, -9, -3, -1]),
    ([2, 2, 1, 1], [1, 1, 3, 0, 1, -5, -2, -5, 1, 0, 3, 1, 1]),
    ([2, 2, 2], [1, -1, 1, -2, 1, -1, 0, 1, -1, 2, -1, 1, -1]),
    ([3, 1, 1, 1], [1, 2, 4, 3, -1, -5, -8, -5, -1, 3, 4, 2, 1]),
    ([3, 2, 1], [1, 0, 0, -1, -1, -1, 0, 1, 1, 1, 0, 0, -1]),
    ([3, 3], [1, -1, -2, 3, -1, -2, 4, -2, -1, 3, -2, -1, 1]),
    ([4, 1, 1], [1, 1, 1, -2, -3, -3, 0, 3, 3, 2, -1, -1, -1]),
    ([4, 2], [1, -1, -1, 0, 1, 1, -2, 1, 1, 0, -1, -1, 1]),
    ([5, 1], [1, 0, -1, -2, -1, 2, 2, 2, -1, -2, -1, 0, 1]),
    ([6], [1, -1, -2, 1, 1, 2, 0, -2, -1, -1, 2, 1, -1]),
    ])

hess_left[(0, 0, 0, 1, 1, 1)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 5, 19, 48, 89, 127, 142, 127, 89, 48, 19, 5, 1]),
    ([2, 1, 1, 1, 1], [1, 5, 17, 40, 71, 99, 110, 99, 71, 40, 17, 5, 1]),
    ([2, 2, 1, 1], [1, 5, 15, 32, 53, 71, 78, 71, 53, 32, 15, 5, 1]),
    ([2, 2, 2], [1, 5, 13, 24, 35, 43, 46, 43, 35, 24, 13, 5, 1]),
    ([3, 1, 1, 1], [1, 5, 16, 36, 62, 85, 94, 85, 62, 36, 16, 5, 1]),
    ([3, 2, 1], [1, 5, 14, 28, 44, 57, 62, 57, 44, 28, 14, 5, 1]),
    ([3, 3], [1, 5, 13, 24, 35, 43, 46, 43, 35, 24, 13, 5, 1]),
    ([4, 1, 1], [1, 5, 15, 32, 53, 71, 78, 71, 53, 32, 15, 5, 1]),
    ([4, 2], [1, 5, 13, 24, 35, 43, 46, 43, 35, 24, 13, 5, 1]),
    ([5, 1], [1, 5, 14, 28, 44, 57, 62, 57, 44, 28, 14, 5, 1]),
    ([6], [1, 5, 13, 24, 35, 43, 46, 43, 35, 24, 13, 5, 1]),
    ])

hess_right[(0, 0, 0, 1, 1, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 5, 19, 53, 112, 170, 170, 112, 53, 19, 5, 1]),
    ([2, 1, 1, 1, 1], [1, 3, 9, 17, 22, 12, -12, -22, -17, -9, -3, -1]),
    ([2, 2, 1, 1], [1, 1, 3, 1, 0, -6, -6, 0, 1, 3, 1, 1]),
    ([2, 2, 2], [1, -1, 1, -3, 6, -12, 12, -6, 3, -1, 1, -1]),
    ([3, 1, 1, 1], [1, 2, 4, 5, -2, -10, -10, -2, 5, 4, 2, 1]),
    ([3, 2, 1], [1, 0, 0, -1, -2, 0, 0, 2, 1, 0, 0, -1]),
    ([3, 3], [1, -1, -2, 2, 1, -1, -1, 1, 2, -2, -1, 1]),
    ([4, 1, 1], [1, 1, 1, -1, -6, -4, 4, 6, 1, -1, -1, -1]),
    ([4, 2], [1, -1, -1, -1, 4, -2, -2, 4, -1, -1, -1, 1]),
    ([5, 1], [1, 0, -1, -2, -3, 5, 5, -3, -2, -1, 0, 1]),
    ([6], [1, -1, -2, 0, 3, 3, -3, -3, 0, 2, 1, -1]),
    ])

hess_left[(0, 0, 0, 1, 1, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 5, 19, 53, 112, 170, 170, 112, 53, 19, 5, 1]),
    ([2, 1, 1, 1, 1], [1, 5, 17, 43, 82, 116, 116, 82, 43, 17, 5, 1]),
    ([2, 2, 1, 1], [1, 5, 15, 33, 56, 74, 74, 56, 33, 15, 5, 1]),
    ([2, 2, 2], [1, 5, 13, 23, 34, 44, 44, 34, 23, 13, 5, 1]),
    ([3, 1, 1, 1], [1, 5, 16, 38, 67, 89, 89, 67, 38, 16, 5, 1]),
    ([3, 2, 1], [1, 5, 14, 28, 43, 53, 53, 43, 28, 14, 5, 1]),
    ([3, 3], [1, 5, 13, 23, 31, 35, 35, 31, 23, 13, 5, 1]),
    ([4, 1, 1], [1, 5, 15, 33, 54, 68, 68, 54, 33, 15, 5, 1]),
    ([4, 2], [1, 5, 13, 23, 32, 38, 38, 32, 23, 13, 5, 1]),
    ([5, 1], [1, 5, 14, 28, 42, 50, 50, 42, 28, 14, 5, 1]),
    ([6], [1, 5, 13, 23, 31, 35, 35, 31, 23, 13, 5, 1]),
    ])

hess_right[(0, 0, 0, 1, 1, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 5, 24, 76, 155, 198, 155, 76, 24, 5, 1]),
    ([2, 1, 1, 1, 1], [1, 3, 12, 22, 23, 0, -23, -22, -12, -3, -1]),
    ([2, 2, 1, 1], [1, 1, 4, 0, -1, -10, -1, 0, 4, 1, 1]),
    ([2, 2, 2], [1, -1, 0, 2, -5, 0, 5, -2, 0, 1, -1]),
    ([3, 1, 1, 1], [1, 2, 6, 4, -7, -12, -7, 4, 6, 2, 1]),
    ([3, 2, 1], [1, 0, 0, -2, -1, 0, 1, 2, 0, 0, -1]),
    ([3, 3], [1, -1, -3, 4, 2, -6, 2, 4, -3, -1, 1]),
    ([4, 1, 1], [1, 1, 2, -4, -7, 0, 7, 4, -2, -1, -1]),
    ([4, 2], [1, -1, -2, 2, 1, -2, 1, 2, -2, -1, 1]),
    ([5, 1], [1, 0, -1, -4, 0, 8, 0, -4, -1, 0, 1]),
    ([6], [1, -1, -3, 2, 4, 0, -4, -2, 3, 1, -1]),
    ])

hess_left[(0, 0, 0, 1, 1, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 5, 24, 76, 155, 198, 155, 76, 24, 5, 1]),
    ([2, 1, 1, 1, 1], [1, 5, 20, 54, 99, 122, 99, 54, 20, 5, 1]),
    ([2, 2, 1, 1], [1, 5, 16, 36, 59, 70, 59, 36, 16, 5, 1]),
    ([2, 2, 2], [1, 5, 12, 22, 35, 42, 35, 22, 12, 5, 1]),
    ([3, 1, 1, 1], [1, 5, 18, 43, 71, 84, 71, 43, 18, 5, 1]),
    ([3, 2, 1], [1, 5, 14, 27, 39, 44, 39, 27, 14, 5, 1]),
    ([3, 3], [1, 5, 12, 19, 23, 24, 23, 19, 12, 5, 1]),
    ([4, 1, 1], [1, 5, 16, 34, 51, 58, 51, 34, 16, 5, 1]),
    ([4, 2], [1, 5, 12, 20, 27, 30, 27, 20, 12, 5, 1]),
    ([5, 1], [1, 5, 14, 26, 35, 38, 35, 26, 14, 5, 1]),
    ([6], [1, 5, 12, 19, 23, 24, 23, 19, 12, 5, 1]),
    ])

hess_right[(0, 0, 0, 1, 1, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 10, 38, 116, 195, 195, 116, 38, 10, 1]),
    ([2, 1, 1, 1, 1], [1, 6, 14, 28, 11, -11, -28, -14, -6, -1]),
    ([2, 2, 1, 1], [1, 2, 2, 0, -5, -5, 0, 2, 2, 1]),
    ([2, 2, 2], [1, -2, 2, 0, -5, 5, 0, -2, 2, -1]),
    ([3, 1, 1, 1], [1, 4, 5, 2, -12, -12, 2, 5, 4, 1]),
    ([3, 2, 1], [1, 0, -1, -2, 2, -2, 2, 1, 0, -1]),
    ([3, 3], [1, -2, -1, 5, -3, -3, 5, -1, -2, 1]),
    ([4, 1, 1], [1, 2, 0, -6, -5, 5, 6, 0, -2, -1]),
    ([4, 2], [1, -2, 0, 2, -1, -1, 2, 0, -2, 1]),
    ([5, 1], [1, 0, -2, -4, 5, 5, -4, -2, 0, 1]),
    ([6], [1, -2, -1, 3, 1, -1, -3, 1, 2, -1]),
    ])

hess_left[(0, 0, 0, 1, 1, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 10, 38, 116, 195, 195, 116, 38, 10, 1]),
    ([2, 1, 1, 1, 1], [1, 8, 28, 70, 109, 109, 70, 28, 8, 1]),
    ([2, 2, 1, 1], [1, 6, 18, 36, 51, 51, 36, 18, 6, 1]),
    ([2, 2, 2], [1, 4, 8, 14, 21, 21, 14, 8, 4, 1]),
    ([3, 1, 1, 1], [1, 7, 23, 50, 72, 72, 50, 23, 7, 1]),
    ([3, 2, 1], [1, 5, 13, 22, 28, 28, 22, 13, 5, 1]),
    ([3, 3], [1, 4, 8, 11, 12, 12, 11, 8, 4, 1]),
    ([4, 1, 1], [1, 6, 18, 34, 45, 45, 34, 18, 6, 1]),
    ([4, 2], [1, 4, 8, 12, 15, 15, 12, 8, 4, 1]),
    ([5, 1], [1, 5, 13, 21, 25, 25, 21, 13, 5, 1]),
    ([6], [1, 4, 8, 11, 12, 12, 11, 8, 4, 1]),
    ])

hess_right[(0, 0, 0, 1, 1, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 24, 78, 156, 192, 156, 78, 24, 6]),
    ([2, 1, 1, 1, 1], [4, 8, 20, 16, 0, -16, -20, -8, -4]),
    ([2, 2, 1, 1], [2, 0, 2, -4, 0, -4, 2, 0, 2]),
    ([3, 1, 1, 1], [3, 3, 3, -3, -12, -3, 3, 3, 3]),
    ([3, 2, 1], [1, -1, -1, 1, 0, -1, 1, 1, -1]),
    ([4, 1, 1], [2, 0, -2, -4, 0, 4, 2, 0, -2]),
    ([5, 1], [1, -1, -2, 1, 2, 1, -2, -1, 1]),
    ])

hess_left[(0, 0, 0, 1, 1, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 24, 78, 156, 192, 156, 78, 24, 6]),
    ([2, 1, 1, 1, 1], [4, 16, 44, 80, 96, 80, 44, 16, 4]),
    ([2, 2, 1, 1], [2, 8, 18, 28, 32, 28, 18, 8, 2]),
    ([3, 1, 1, 1], [3, 12, 30, 51, 60, 51, 30, 12, 3]),
    ([3, 2, 1], [1, 4, 8, 11, 12, 11, 8, 4, 1]),
    ([4, 1, 1], [2, 8, 18, 28, 32, 28, 18, 8, 2]),
    ([5, 1], [1, 4, 8, 11, 12, 11, 8, 4, 1]),
    ])

hess_right[(0, 0, 0, 1, 2, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 5, 24, 76, 155, 198, 155, 76, 24, 5, 1]),
    ([2, 1, 1, 1, 1], [1, 3, 12, 22, 23, 0, -23, -22, -12, -3, -1]),
    ([2, 2, 1, 1], [1, 1, 4, 0, -1, -10, -1, 0, 4, 1, 1]),
    ([2, 2, 2], [1, -1, 0, 2, -5, 0, 5, -2, 0, 1, -1]),
    ([3, 1, 1, 1], [1, 2, 6, 4, -7, -12, -7, 4, 6, 2, 1]),
    ([3, 2, 1], [1, 0, 0, -2, -1, 0, 1, 2, 0, 0, -1]),
    ([3, 3], [1, -1, -3, 4, 2, -6, 2, 4, -3, -1, 1]),
    ([4, 1, 1], [1, 1, 2, -4, -7, 0, 7, 4, -2, -1, -1]),
    ([4, 2], [1, -1, -2, 2, 1, -2, 1, 2, -2, -1, 1]),
    ([5, 1], [1, 0, -1, -4, 0, 8, 0, -4, -1, 0, 1]),
    ([6], [1, -1, -3, 2, 4, 0, -4, -2, 3, 1, -1]),
    ])

hess_left[(0, 0, 0, 1, 2, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 5, 24, 76, 155, 198, 155, 76, 24, 5, 1]),
    ([2, 1, 1, 1, 1], [1, 5, 20, 54, 99, 122, 99, 54, 20, 5, 1]),
    ([2, 2, 1, 1], [1, 5, 16, 36, 59, 70, 59, 36, 16, 5, 1]),
    ([2, 2, 2], [1, 5, 12, 22, 35, 42, 35, 22, 12, 5, 1]),
    ([3, 1, 1, 1], [1, 5, 18, 43, 71, 84, 71, 43, 18, 5, 1]),
    ([3, 2, 1], [1, 5, 14, 27, 39, 44, 39, 27, 14, 5, 1]),
    ([3, 3], [1, 5, 12, 19, 23, 24, 23, 19, 12, 5, 1]),
    ([4, 1, 1], [1, 5, 16, 34, 51, 58, 51, 34, 16, 5, 1]),
    ([4, 2], [1, 5, 12, 20, 27, 30, 27, 20, 12, 5, 1]),
    ([5, 1], [1, 5, 14, 26, 35, 38, 35, 26, 14, 5, 1]),
    ([6], [1, 5, 12, 19, 23, 24, 23, 19, 12, 5, 1]),
    ])

hess_right[(0, 0, 0, 1, 2, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 5, 29, 113, 212, 212, 113, 29, 5, 1]),
    ([2, 1, 1, 1, 1], [1, 3, 15, 31, 18, -18, -31, -15, -3, -1]),
    ([2, 2, 1, 1], [1, 1, 5, 1, -8, -8, 1, 5, 1, 1]),
    ([2, 2, 2], [1, -1, -1, 7, -18, 18, -7, 1, 1, -1]),
    ([3, 1, 1, 1], [1, 2, 8, 2, -13, -13, 2, 8, 2, 1]),
    ([3, 2, 1], [1, 0, 0, -2, -3, 3, 2, 0, 0, -1]),
    ([3, 3], [1, -1, -4, 8, -4, -4, 8, -4, -1, 1]),
    ([4, 1, 1], [1, 1, 3, -9, -4, 4, 9, -3, -1, -1]),
    ([4, 2], [1, -1, -3, 5, -2, -2, 5, -3, -1, 1]),
    ([5, 1], [1, 0, -1, -7, 7, 7, -7, -1, 0, 1]),
    ([6], [1, -1, -4, 4, 6, -6, -4, 4, 1, -1]),
    ])

hess_left[(0, 0, 0, 1, 2, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 5, 29, 113, 212, 212, 113, 29, 5, 1]),
    ([2, 1, 1, 1, 1], [1, 5, 23, 69, 118, 118, 69, 23, 5, 1]),
    ([2, 2, 1, 1], [1, 5, 17, 41, 64, 64, 41, 17, 5, 1]),
    ([2, 2, 2], [1, 5, 11, 21, 34, 34, 21, 11, 5, 1]),
    ([3, 1, 1, 1], [1, 5, 20, 47, 71, 71, 47, 20, 5, 1]),
    ([3, 2, 1], [1, 5, 14, 27, 37, 37, 27, 14, 5, 1]),
    ([3, 3], [1, 5, 11, 17, 20, 20, 17, 11, 5, 1]),
    ([4, 1, 1], [1, 5, 17, 33, 44, 44, 33, 17, 5, 1]),
    ([4, 2], [1, 5, 11, 17, 22, 22, 17, 11, 5, 1]),
    ([5, 1], [1, 5, 14, 23, 27, 27, 23, 14, 5, 1]),
    ([6], [1, 5, 11, 15, 16, 16, 15, 11, 5, 1]),
    ])

hess_right[(0, 0, 0, 1, 2, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 10, 52, 166, 262, 166, 52, 10, 1]),
    ([2, 1, 1, 1, 1], [1, 6, 20, 34, 0, -34, -20, -6, -1]),
    ([2, 2, 1, 1], [1, 2, 4, -2, -10, -2, 4, 2, 1]),
    ([2, 2, 2], [1, -2, 4, -6, 0, 6, -4, 2, -1]),
    ([3, 1, 1, 1], [1, 4, 7, -2, -20, -2, 7, 4, 1]),
    ([3, 2, 1], [1, 0, -1, -2, 0, 2, 1, 0, -1]),
    ([3, 3], [1, -2, -2, 10, -14, 10, -2, -2, 1]),
    ([4, 1, 1], [1, 2, 0, -10, 0, 10, 0, -2, -1]),
    ([4, 2], [1, -2, 0, 2, -2, 2, 0, -2, 1]),
    ([5, 1], [1, 0, -3, -4, 12, -4, -3, 0, 1]),
    ([6], [1, -2, -2, 6, 0, -6, 2, 2, -1]),
    ])

hess_left[(0, 0, 0, 1, 2, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 10, 52, 166, 262, 166, 52, 10, 1]),
    ([2, 1, 1, 1, 1], [1, 8, 34, 88, 122, 88, 34, 8, 1]),
    ([2, 2, 1, 1], [1, 6, 20, 42, 54, 42, 20, 6, 1]),
    ([2, 2, 2], [1, 4, 10, 20, 26, 20, 10, 4, 1]),
    ([3, 1, 1, 1], [1, 7, 25, 52, 64, 52, 25, 7, 1]),
    ([3, 2, 1], [1, 5, 13, 22, 26, 22, 13, 5, 1]),
    ([3, 3], [1, 4, 7, 10, 10, 10, 7, 4, 1]),
    ([4, 1, 1], [1, 6, 18, 30, 34, 30, 18, 6, 1]),
    ([4, 2], [1, 4, 8, 12, 14, 12, 8, 4, 1]),
    ([5, 1], [1, 5, 12, 16, 17, 16, 12, 5, 1]),
    ([6], [1, 4, 7, 8, 8, 8, 7, 4, 1]),
    ])

hess_right[(0, 0, 0, 1, 2, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 24, 102, 228, 228, 102, 24, 6]),
    ([2, 1, 1, 1, 1], [4, 8, 28, 16, -16, -28, -8, -4]),
    ([2, 2, 1, 1], [2, 0, 2, -4, -4, 2, 0, 2]),
    ([3, 1, 1, 1], [3, 3, 6, -12, -12, 6, 3, 3]),
    ([3, 2, 1], [1, -1, -2, 4, -4, 2, 1, -1]),
    ([4, 1, 1], [2, 0, -2, -8, 8, 2, 0, -2]),
    ([5, 1], [1, -1, -3, 3, 3, -3, -1, 1]),
    ])

hess_left[(0, 0, 0, 1, 2, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 24, 102, 228, 228, 102, 24, 6]),
    ([2, 1, 1, 1, 1], [4, 16, 52, 96, 96, 52, 16, 4]),
    ([2, 2, 1, 1], [2, 8, 18, 28, 28, 18, 8, 2]),
    ([3, 1, 1, 1], [3, 12, 33, 51, 51, 33, 12, 3]),
    ([3, 2, 1], [1, 4, 7, 9, 9, 7, 4, 1]),
    ([4, 1, 1], [2, 8, 18, 24, 24, 18, 8, 2]),
    ([5, 1], [1, 4, 7, 8, 8, 7, 4, 1]),
    ])

hess_right[(0, 0, 0, 1, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 10, 66, 170, 226, 170, 66, 10, 1]),
    ([2, 1, 1, 1, 1], [1, 6, 24, 26, 0, -26, -24, -6, -1]),
    ([2, 2, 1, 1], [1, 2, 6, -6, -6, -6, 6, 2, 1]),
    ([2, 2, 2], [1, -2, 4, -6, 0, 6, -4, 2, -1]),
    ([3, 1, 1, 1], [1, 4, 6, -4, -14, -4, 6, 4, 1]),
    ([3, 2, 1], [1, 0, 0, -4, 0, 4, 0, 0, -1]),
    ([3, 3], [1, -2, 0, 2, -2, 2, 0, -2, 1]),
    ([4, 1, 1], [1, 2, -2, -6, 0, 6, 2, -2, -1]),
    ([4, 2], [1, -2, 0, 2, -2, 2, 0, -2, 1]),
    ([5, 1], [1, 0, -4, 0, 6, 0, -4, 0, 1]),
    ([6], [1, -2, -2, 6, 0, -6, 2, 2, -1]),
    ])

hess_left[(0, 0, 0, 1, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 10, 66, 170, 226, 170, 66, 10, 1]),
    ([2, 1, 1, 1, 1], [1, 8, 38, 88, 114, 88, 38, 8, 1]),
    ([2, 2, 1, 1], [1, 6, 22, 46, 58, 46, 22, 6, 1]),
    ([2, 2, 2], [1, 4, 10, 20, 26, 20, 10, 4, 1]),
    ([3, 1, 1, 1], [1, 7, 24, 47, 58, 47, 24, 7, 1]),
    ([3, 2, 1], [1, 5, 14, 25, 30, 25, 14, 5, 1]),
    ([3, 3], [1, 4, 9, 14, 16, 14, 9, 4, 1]),
    ([4, 1, 1], [1, 6, 16, 26, 30, 26, 16, 6, 1]),
    ([4, 2], [1, 4, 8, 12, 14, 12, 8, 4, 1]),
    ([5, 1], [1, 5, 11, 15, 16, 15, 11, 5, 1]),
    ([6], [1, 4, 7, 8, 8, 8, 7, 4, 1]),
    ])

hess_right[(0, 0, 0, 1, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 24, 102, 233, 233, 102, 24, 1]),
    ([2, 1, 1, 1, 1], [1, 12, 26, 23, -23, -26, -12, -1]),
    ([2, 2, 1, 1], [1, 4, 2, -7, -7, 2, 4, 1]),
    ([2, 2, 2], [1, 0, -2, -1, 1, 2, 0, -1]),
    ([3, 1, 1, 1], [1, 6, 3, -10, -10, 3, 6, 1]),
    ([3, 2, 1], [1, 0, -1, -4, 4, 1, 0, -1]),
    ([3, 3], [1, -3, 3, -1, -1, 3, -3, 1]),
    ([4, 1, 1], [1, 2, -4, -5, 5, 4, -2, -1]),
    ([4, 2], [1, -2, 0, 1, 1, 0, -2, 1]),
    ([5, 1], [1, -1, -3, 3, 3, -3, -1, 1]),
    ([6], [1, -3, 1, 5, -5, -1, 3, -1]),
    ])

hess_left[(0, 0, 0, 1, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 24, 102, 233, 233, 102, 24, 1]),
    ([2, 1, 1, 1, 1], [1, 14, 52, 101, 101, 52, 14, 1]),
    ([2, 2, 1, 1], [1, 8, 26, 45, 45, 26, 8, 1]),
    ([2, 2, 2], [1, 6, 16, 25, 25, 16, 6, 1]),
    ([3, 1, 1, 1], [1, 9, 27, 44, 44, 27, 9, 1]),
    ([3, 2, 1], [1, 5, 13, 20, 20, 13, 5, 1]),
    ([3, 3], [1, 3, 6, 8, 8, 6, 3, 1]),
    ([4, 1, 1], [1, 6, 14, 19, 19, 14, 6, 1]),
    ([4, 2], [1, 4, 8, 11, 11, 8, 4, 1]),
    ([5, 1], [1, 4, 7, 8, 8, 7, 4, 1]),
    ([6], [1, 3, 4, 4, 4, 4, 3, 1]),
    ])

hess_right[(0, 0, 0, 1, 3, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 48, 174, 264, 174, 48, 6]),
    ([2, 1, 1, 1, 1], [4, 16, 28, 0, -28, -16, -4]),
    ([2, 2, 1, 1], [2, 0, 2, -8, 2, 0, 2]),
    ([3, 1, 1, 1], [3, 6, -3, -12, -3, 6, 3]),
    ([3, 2, 1], [1, -2, 1, 0, -1, 2, -1]),
    ([4, 1, 1], [2, 0, -6, 0, 6, 0, -2]),
    ([5, 1], [1, -2, -1, 4, -1, -2, 1]),
    ])

hess_left[(0, 0, 0, 1, 3, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 48, 174, 264, 174, 48, 6]),
    ([2, 1, 1, 1, 1], [4, 24, 68, 96, 68, 24, 4]),
    ([2, 2, 1, 1], [2, 8, 18, 24, 18, 8, 2]),
    ([3, 1, 1, 1], [3, 15, 33, 42, 33, 15, 3]),
    ([3, 2, 1], [1, 3, 5, 6, 5, 3, 1]),
    ([4, 1, 1], [2, 8, 14, 16, 14, 8, 2]),
    ([5, 1], [1, 3, 4, 4, 4, 3, 1]),
    ])

hess_right[(0, 0, 0, 1, 4, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [15, 60, 165, 240, 165, 60, 15]),
    ([2, 1, 1, 1, 1], [7, 14, 23, 0, -23, -14, -7]),
    ([2, 2, 1, 1], [3, 0, 1, -8, 1, 0, 3]),
    ([2, 2, 2], [3, -6, 3, 0, -3, 6, -3]),
    ([3, 1, 1, 1], [3, 3, -3, -6, -3, 3, 3]),
    ([3, 2, 1], [1, -1, -1, 0, 1, 1, -1]),
    ([4, 1, 1], [1, 0, -3, 0, 3, 0, -1]),
    ([4, 2], [1, -2, -1, 4, -1, -2, 1]),
    ])

hess_left[(0, 0, 0, 1, 4, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [15, 60, 165, 240, 165, 60, 15]),
    ([2, 1, 1, 1, 1], [7, 28, 65, 88, 65, 28, 7]),
    ([2, 2, 1, 1], [3, 12, 25, 32, 25, 12, 3]),
    ([2, 2, 2], [3, 12, 21, 24, 21, 12, 3]),
    ([3, 1, 1, 1], [3, 12, 24, 30, 24, 12, 3]),
    ([3, 2, 1], [1, 4, 8, 10, 8, 4, 1]),
    ([4, 1, 1], [1, 4, 7, 8, 7, 4, 1]),
    ([4, 2], [1, 4, 7, 8, 7, 4, 1]),
    ])

hess_right[(0, 0, 0, 1, 4, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [30, 90, 240, 240, 90, 30]),
    ([2, 1, 1, 1, 1], [12, 12, 24, -24, -12, -12]),
    ([2, 2, 1, 1], [2, -2, 0, 0, -2, 2]),
    ([3, 1, 1, 1], [6, 0, -6, -6, 0, 6]),
    ([4, 1, 1], [2, -2, -4, 4, 2, -2]),
    ])

hess_left[(0, 0, 0, 1, 4, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [30, 90, 240, 240, 90, 30]),
    ([2, 1, 1, 1, 1], [12, 36, 72, 72, 36, 12]),
    ([2, 2, 1, 1], [2, 6, 8, 8, 6, 2]),
    ([3, 1, 1, 1], [6, 18, 30, 30, 18, 6]),
    ([4, 1, 1], [2, 6, 8, 8, 6, 2]),
    ])

hess_right[(0, 0, 0, 2, 2, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 10, 47, 119, 183, 183, 119, 47, 10, 1]),
    ([2, 1, 1, 1, 1], [1, 6, 17, 23, 11, -11, -23, -17, -6, -1]),
    ([2, 2, 1, 1], [1, 2, 3, -1, -5, -5, -1, 3, 2, 1]),
    ([2, 2, 2], [1, -2, 5, -9, 7, -7, 9, -5, 2, -1]),
    ([3, 1, 1, 1], [1, 4, 5, -1, -9, -9, -1, 5, 4, 1]),
    ([3, 2, 1], [1, 0, -1, -1, -1, 1, 1, 1, 0, -1]),
    ([3, 3], [1, -2, -1, 5, -3, -3, 5, -1, -2, 1]),
    ([4, 1, 1], [1, 2, -1, -5, -3, 3, 5, 1, -2, -1]),
    ([4, 2], [1, -2, 1, -1, 1, 1, -1, 1, -2, 1]),
    ([5, 1], [1, 0, -3, -1, 3, 3, -1, -3, 0, 1]),
    ([6], [1, -2, -1, 3, 1, -1, -3, 1, 2, -1]),
    ])

hess_left[(0, 0, 0, 2, 2, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 10, 47, 119, 183, 183, 119, 47, 10, 1]),
    ([2, 1, 1, 1, 1], [1, 8, 31, 71, 105, 105, 71, 31, 8, 1]),
    ([2, 2, 1, 1], [1, 6, 19, 39, 55, 55, 39, 19, 6, 1]),
    ([2, 2, 2], [1, 4, 11, 23, 33, 33, 23, 11, 4, 1]),
    ([3, 1, 1, 1], [1, 7, 23, 47, 66, 66, 47, 23, 7, 1]),
    ([3, 2, 1], [1, 5, 13, 23, 30, 30, 23, 13, 5, 1]),
    ([3, 3], [1, 4, 8, 11, 12, 12, 11, 8, 4, 1]),
    ([4, 1, 1], [1, 6, 17, 31, 41, 41, 31, 17, 6, 1]),
    ([4, 2], [1, 4, 9, 15, 19, 19, 15, 9, 4, 1]),
    ([5, 1], [1, 5, 12, 19, 23, 23, 19, 12, 5, 1]),
    ([6], [1, 4, 8, 11, 12, 12, 11, 8, 4, 1]),
    ])

hess_right[(0, 0, 0, 2, 2, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 10, 66, 170, 226, 170, 66, 10, 1]),
    ([2, 1, 1, 1, 1], [1, 6, 24, 26, 0, -26, -24, -6, -1]),
    ([2, 2, 1, 1], [1, 2, 6, -6, -6, -6, 6, 2, 1]),
    ([2, 2, 2], [1, -2, 4, -6, 0, 6, -4, 2, -1]),
    ([3, 1, 1, 1], [1, 4, 6, -4, -14, -4, 6, 4, 1]),
    ([3, 2, 1], [1, 0, 0, -4, 0, 4, 0, 0, -1]),
    ([3, 3], [1, -2, 0, 2, -2, 2, 0, -2, 1]),
    ([4, 1, 1], [1, 2, -2, -6, 0, 6, 2, -2, -1]),
    ([4, 2], [1, -2, 0, 2, -2, 2, 0, -2, 1]),
    ([5, 1], [1, 0, -4, 0, 6, 0, -4, 0, 1]),
    ([6], [1, -2, -2, 6, 0, -6, 2, 2, -1]),
    ])

hess_left[(0, 0, 0, 2, 2, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 10, 66, 170, 226, 170, 66, 10, 1]),
    ([2, 1, 1, 1, 1], [1, 8, 38, 88, 114, 88, 38, 8, 1]),
    ([2, 2, 1, 1], [1, 6, 22, 46, 58, 46, 22, 6, 1]),
    ([2, 2, 2], [1, 4, 10, 20, 26, 20, 10, 4, 1]),
    ([3, 1, 1, 1], [1, 7, 24, 47, 58, 47, 24, 7, 1]),
    ([3, 2, 1], [1, 5, 14, 25, 30, 25, 14, 5, 1]),
    ([3, 3], [1, 4, 9, 14, 16, 14, 9, 4, 1]),
    ([4, 1, 1], [1, 6, 16, 26, 30, 26, 16, 6, 1]),
    ([4, 2], [1, 4, 8, 12, 14, 12, 8, 4, 1]),
    ([5, 1], [1, 5, 11, 15, 16, 15, 11, 5, 1]),
    ([6], [1, 4, 7, 8, 8, 8, 7, 4, 1]),
    ])

hess_right[(0, 0, 0, 2, 2, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 15, 99, 245, 245, 99, 15, 1]),
    ([2, 1, 1, 1, 1], [1, 9, 31, 23, -23, -31, -9, -1]),
    ([2, 2, 1, 1], [1, 3, 3, -7, -7, 3, 3, 1]),
    ([2, 2, 2], [1, -3, 7, -13, 13, -7, 3, -1]),
    ([3, 1, 1, 1], [1, 6, 6, -13, -13, 6, 6, 1]),
    ([3, 2, 1], [1, 0, -2, -1, 1, 2, 0, -1]),
    ([3, 3], [1, -3, 3, -1, -1, 3, -3, 1]),
    ([4, 1, 1], [1, 3, -5, -7, 7, 5, -3, -1]),
    ([4, 2], [1, -3, 3, -1, -1, 3, -3, 1]),
    ([5, 1], [1, 0, -6, 5, 5, -6, 0, 1]),
    ([6], [1, -3, 1, 5, -5, -1, 3, -1]),
    ])

hess_left[(0, 0, 0, 2, 2, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 15, 99, 245, 245, 99, 15, 1]),
    ([2, 1, 1, 1, 1], [1, 11, 51, 105, 105, 51, 11, 1]),
    ([2, 2, 1, 1], [1, 7, 23, 41, 41, 23, 7, 1]),
    ([2, 2, 2], [1, 3, 7, 13, 13, 7, 3, 1]),
    ([3, 1, 1, 1], [1, 9, 30, 50, 50, 30, 9, 1]),
    ([3, 2, 1], [1, 5, 12, 18, 18, 12, 5, 1]),
    ([3, 3], [1, 3, 6, 8, 8, 6, 3, 1]),
    ([4, 1, 1], [1, 7, 17, 23, 23, 17, 7, 1]),
    ([4, 2], [1, 3, 5, 7, 7, 5, 3, 1]),
    ([5, 1], [1, 5, 9, 10, 10, 9, 5, 1]),
    ([6], [1, 3, 4, 4, 4, 4, 3, 1]),
    ])

hess_right[(0, 0, 0, 2, 2, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 48, 174, 264, 174, 48, 6]),
    ([2, 1, 1, 1, 1], [4, 16, 28, 0, -28, -16, -4]),
    ([2, 2, 1, 1], [2, 0, 2, -8, 2, 0, 2]),
    ([3, 1, 1, 1], [3, 6, -3, -12, -3, 6, 3]),
    ([3, 2, 1], [1, -2, 1, 0, -1, 2, -1]),
    ([4, 1, 1], [2, 0, -6, 0, 6, 0, -2]),
    ([5, 1], [1, -2, -1, 4, -1, -2, 1]),
    ])

hess_left[(0, 0, 0, 2, 2, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 48, 174, 264, 174, 48, 6]),
    ([2, 1, 1, 1, 1], [4, 24, 68, 96, 68, 24, 4]),
    ([2, 2, 1, 1], [2, 8, 18, 24, 18, 8, 2]),
    ([3, 1, 1, 1], [3, 15, 33, 42, 33, 15, 3]),
    ([3, 2, 1], [1, 3, 5, 6, 5, 3, 1]),
    ([4, 1, 1], [2, 8, 14, 16, 14, 8, 2]),
    ([5, 1], [1, 3, 4, 4, 4, 3, 1]),
    ])

hess_right[(0, 0, 0, 2, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 29, 117, 213, 213, 117, 29, 1]),
    ([2, 1, 1, 1, 1], [1, 13, 27, 15, -15, -27, -13, -1]),
    ([2, 2, 1, 1], [1, 5, 1, -7, -7, 1, 5, 1]),
    ([2, 2, 2], [1, -3, 7, -13, 13, -7, 3, -1]),
    ([3, 1, 1, 1], [1, 5, 3, -9, -9, 3, 5, 1]),
    ([3, 2, 1], [1, 1, -3, -3, 3, 3, -1, -1]),
    ([3, 3], [1, -1, -3, 3, 3, -3, -1, 1]),
    ([4, 1, 1], [1, 1, -3, -3, 3, 3, -1, -1]),
    ([4, 2], [1, -3, 3, -1, -1, 3, -3, 1]),
    ([5, 1], [1, -1, -3, 3, 3, -3, -1, 1]),
    ([6], [1, -3, 1, 5, -5, -1, 3, -1]),
    ])

hess_left[(0, 0, 0, 2, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 29, 117, 213, 213, 117, 29, 1]),
    ([2, 1, 1, 1, 1], [1, 15, 55, 97, 97, 55, 15, 1]),
    ([2, 2, 1, 1], [1, 9, 29, 49, 49, 29, 9, 1]),
    ([2, 2, 2], [1, 3, 7, 13, 13, 7, 3, 1]),
    ([3, 1, 1, 1], [1, 8, 24, 39, 39, 24, 8, 1]),
    ([3, 2, 1], [1, 6, 16, 25, 25, 16, 6, 1]),
    ([3, 3], [1, 5, 12, 18, 18, 12, 5, 1]),
    ([4, 1, 1], [1, 5, 11, 15, 15, 11, 5, 1]),
    ([4, 2], [1, 3, 5, 7, 7, 5, 3, 1]),
    ([5, 1], [1, 4, 7, 8, 8, 7, 4, 1]),
    ([6], [1, 3, 4, 4, 4, 4, 3, 1]),
    ])

hess_right[(0, 0, 0, 2, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 43, 179, 274, 179, 43, 1]),
    ([2, 1, 1, 1, 1], [1, 19, 31, 0, -31, -19, -1]),
    ([2, 2, 1, 1], [1, 7, -5, -6, -5, 7, 1]),
    ([2, 2, 2], [1, -1, -1, 0, 1, 1, -1]),
    ([3, 1, 1, 1], [1, 7, -1, -14, -1, 7, 1]),
    ([3, 2, 1], [1, 1, -5, 0, 5, -1, -1]),
    ([3, 3], [1, -2, -1, 4, -1, -2, 1]),
    ([4, 1, 1], [1, 1, -5, 0, 5, -1, -1]),
    ([4, 2], [1, -3, 3, -2, 3, -3, 1]),
    ([5, 1], [1, -2, -1, 4, -1, -2, 1]),
    ([6], [1, -4, 5, 0, -5, 4, -1]),
    ])

hess_left[(0, 0, 0, 2, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 43, 179, 274, 179, 43, 1]),
    ([2, 1, 1, 1, 1], [1, 21, 71, 102, 71, 21, 1]),
    ([2, 2, 1, 1], [1, 11, 31, 42, 31, 11, 1]),
    ([2, 2, 2], [1, 5, 11, 14, 11, 5, 1]),
    ([3, 1, 1, 1], [1, 10, 26, 34, 26, 10, 1]),
    ([3, 2, 1], [1, 6, 14, 18, 14, 6, 1]),
    ([3, 3], [1, 4, 8, 10, 8, 4, 1]),
    ([4, 1, 1], [1, 5, 9, 10, 9, 5, 1]),
    ([4, 2], [1, 3, 5, 6, 5, 3, 1]),
    ([5, 1], [1, 3, 4, 4, 4, 3, 1]),
    ([6], [1, 2, 2, 2, 2, 2, 1]),
    ])

hess_right[(0, 0, 0, 2, 3, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 102, 252, 252, 102, 6]),
    ([2, 1, 1, 1, 1], [4, 28, 16, -16, -28, -4]),
    ([2, 2, 1, 1], [2, 2, -4, -4, 2, 2]),
    ([3, 1, 1, 1], [3, 6, -9, -9, 6, 3]),
    ([3, 2, 1], [1, -2, 1, -1, 2, -1]),
    ([4, 1, 1], [2, -2, -4, 4, 2, -2]),
    ([5, 1], [1, -3, 2, 2, -3, 1]),
    ])

hess_left[(0, 0, 0, 2, 3, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 102, 252, 252, 102, 6]),
    ([2, 1, 1, 1, 1], [4, 36, 80, 80, 36, 4]),
    ([2, 2, 1, 1], [2, 10, 20, 20, 10, 2]),
    ([3, 1, 1, 1], [3, 15, 27, 27, 15, 3]),
    ([3, 2, 1], [1, 3, 5, 5, 3, 1]),
    ([4, 1, 1], [2, 6, 8, 8, 6, 2]),
    ([5, 1], [1, 2, 2, 2, 2, 1]),
    ])

hess_right[(0, 0, 0, 2, 4, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [15, 105, 240, 240, 105, 15]),
    ([2, 1, 1, 1, 1], [7, 23, 16, -16, -23, -7]),
    ([2, 2, 1, 1], [3, 1, -4, -4, 1, 3]),
    ([2, 2, 2], [3, -9, 12, -12, 9, -3]),
    ([3, 1, 1, 1], [3, 3, -6, -6, 3, 3]),
    ([3, 2, 1], [1, -1, -2, 2, 1, -1]),
    ([4, 1, 1], [1, -1, -2, 2, 1, -1]),
    ([4, 2], [1, -3, 2, 2, -3, 1]),
    ])

hess_left[(0, 0, 0, 2, 4, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [15, 105, 240, 240, 105, 15]),
    ([2, 1, 1, 1, 1], [7, 37, 76, 76, 37, 7]),
    ([2, 2, 1, 1], [3, 13, 24, 24, 13, 3]),
    ([2, 2, 2], [3, 9, 12, 12, 9, 3]),
    ([3, 1, 1, 1], [3, 12, 21, 21, 12, 3]),
    ([3, 2, 1], [1, 4, 7, 7, 4, 1]),
    ([4, 1, 1], [1, 3, 4, 4, 3, 1]),
    ([4, 2], [1, 3, 4, 4, 3, 1]),
    ])

hess_right[(0, 0, 0, 2, 4, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [30, 180, 300, 180, 30]),
    ([2, 1, 1, 1, 1], [12, 24, 0, -24, -12]),
    ([2, 2, 1, 1], [2, -4, 4, -4, 2]),
    ([3, 1, 1, 1], [6, 0, -12, 0, 6]),
    ([4, 1, 1], [2, -4, 0, 4, -2]),
    ])

hess_left[(0, 0, 0, 2, 4, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [30, 180, 300, 180, 30]),
    ([2, 1, 1, 1, 1], [12, 48, 72, 48, 12]),
    ([2, 2, 1, 1], [2, 4, 4, 4, 2]),
    ([3, 1, 1, 1], [6, 18, 24, 18, 6]),
    ([4, 1, 1], [2, 4, 4, 4, 2]),
    ])

hess_right[(0, 0, 0, 3, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [20, 80, 160, 200, 160, 80, 20]),
    ([2, 1, 1, 1, 1], [8, 16, 16, 0, -16, -16, -8]),
    ([2, 2, 1, 1], [4, 0, 0, -8, 0, 0, 4]),
    ([3, 1, 1, 1], [2, 2, -2, -4, -2, 2, 2]),
    ([3, 2, 1], [2, -2, -2, 0, 2, 2, -2]),
    ([3, 3], [2, -4, -2, 8, -2, -4, 2]),
    ])

hess_left[(0, 0, 0, 3, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [20, 80, 160, 200, 160, 80, 20]),
    ([2, 1, 1, 1, 1], [8, 32, 64, 80, 64, 32, 8]),
    ([2, 2, 1, 1], [4, 16, 32, 40, 32, 16, 4]),
    ([3, 1, 1, 1], [2, 8, 16, 20, 16, 8, 2]),
    ([3, 2, 1], [2, 8, 16, 20, 16, 8, 2]),
    ([3, 3], [2, 8, 16, 20, 16, 8, 2]),
    ])

hess_right[(0, 0, 0, 3, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [20, 120, 220, 220, 120, 20]),
    ([2, 1, 1, 1, 1], [8, 24, 8, -8, -24, -8]),
    ([2, 2, 1, 1], [4, 0, -4, -4, 0, 4]),
    ([3, 1, 1, 1], [2, 3, -5, -5, 3, 2]),
    ([3, 2, 1], [2, -3, -1, 1, 3, -2]),
    ([3, 3], [2, -6, 4, 4, -6, 2]),
    ])

hess_left[(0, 0, 0, 3, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [20, 120, 220, 220, 120, 20]),
    ([2, 1, 1, 1, 1], [8, 40, 72, 72, 40, 8]),
    ([2, 2, 1, 1], [4, 16, 28, 28, 16, 4]),
    ([3, 1, 1, 1], [2, 9, 16, 16, 9, 2]),
    ([3, 2, 1], [2, 7, 12, 12, 7, 2]),
    ([3, 3], [2, 6, 10, 10, 6, 2]),
    ])

hess_right[(0, 0, 0, 3, 3, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [60, 180, 240, 180, 60]),
    ([2, 1, 1, 1, 1], [16, 16, 0, -16, -16]),
    ([2, 2, 1, 1], [4, -4, 0, -4, 4]),
    ([3, 1, 1, 1], [3, 0, -6, 0, 3]),
    ([3, 2, 1], [1, -2, 0, 2, -1]),
    ])

hess_left[(0, 0, 0, 3, 3, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [60, 180, 240, 180, 60]),
    ([2, 1, 1, 1, 1], [16, 48, 64, 48, 16]),
    ([2, 2, 1, 1], [4, 12, 16, 12, 4]),
    ([3, 1, 1, 1], [3, 9, 12, 9, 3]),
    ([3, 2, 1], [1, 3, 4, 3, 1]),
    ])

hess_right[(0, 0, 0, 3, 4, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [60, 180, 240, 180, 60]),
    ([2, 1, 1, 1, 1], [16, 16, 0, -16, -16]),
    ([2, 2, 1, 1], [4, -4, 0, -4, 4]),
    ([3, 1, 1, 1], [3, 0, -6, 0, 3]),
    ([3, 2, 1], [1, -2, 0, 2, -1]),
    ])

hess_left[(0, 0, 0, 3, 4, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [60, 180, 240, 180, 60]),
    ([2, 1, 1, 1, 1], [16, 48, 64, 48, 16]),
    ([2, 2, 1, 1], [4, 12, 16, 12, 4]),
    ([3, 1, 1, 1], [3, 9, 12, 9, 3]),
    ([3, 2, 1], [1, 3, 4, 3, 1]),
    ])

hess_right[(0, 0, 0, 3, 4, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [120, 240, 240, 120]),
    ([2, 1, 1, 1, 1], [24, 0, 0, -24]),
    ([3, 1, 1, 1], [6, -6, -6, 6]),
    ])

hess_left[(0, 0, 0, 3, 4, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [120, 240, 240, 120]),
    ([2, 1, 1, 1, 1], [24, 48, 48, 24]),
    ([3, 1, 1, 1], [6, 12, 12, 6]),
    ])

hess_right[(0, 0, 1, 1, 1, 1)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 10, 33, 69, 110, 137, 137, 110, 69, 33, 10, 1]),
    ([2, 1, 1, 1, 1], [1, 6, 11, 15, 14, 5, -5, -14, -15, -11, -6, -1]),
    ([2, 2, 1, 1], [1, 2, 1, 1, -2, -3, -3, -2, 1, 1, 2, 1]),
    ([2, 2, 2], [1, -2, 3, -5, 6, -7, 7, -6, 5, -3, 2, -1]),
    ([3, 1, 1, 1], [1, 4, 3, 0, -1, -7, -7, -1, 0, 3, 4, 1]),
    ([3, 2, 1], [1, 0, -1, 0, -1, -1, 1, 1, 0, 1, 0, -1]),
    ([3, 3], [1, -2, 0, 3, -4, 2, 2, -4, 3, 0, -2, 1]),
    ([4, 1, 1], [1, 2, -1, -3, -2, -1, 1, 2, 3, 1, -2, -1]),
    ([4, 2], [1, -2, 1, -1, 2, -1, -1, 2, -1, 1, -2, 1]),
    ([5, 1], [1, 0, -2, -1, 0, 2, 2, 0, -1, -2, 0, 1]),
    ([6], [1, -2, 0, 1, 0, 2, -2, 0, -1, 0, 2, -1]),
    ])

hess_left[(0, 0, 1, 1, 1, 1)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 10, 33, 69, 110, 137, 137, 110, 69, 33, 10, 1]),
    ([2, 1, 1, 1, 1], [1, 8, 25, 51, 80, 99, 99, 80, 51, 25, 8, 1]),
    ([2, 2, 1, 1], [1, 6, 17, 33, 50, 61, 61, 50, 33, 17, 6, 1]),
    ([2, 2, 2], [1, 4, 9, 15, 20, 23, 23, 20, 15, 9, 4, 1]),
    ([3, 1, 1, 1], [1, 7, 21, 42, 65, 80, 80, 65, 42, 21, 7, 1]),
    ([3, 2, 1], [1, 5, 13, 24, 35, 42, 42, 35, 24, 13, 5, 1]),
    ([3, 3], [1, 4, 9, 15, 20, 23, 23, 20, 15, 9, 4, 1]),
    ([4, 1, 1], [1, 6, 17, 33, 50, 61, 61, 50, 33, 17, 6, 1]),
    ([4, 2], [1, 4, 9, 15, 20, 23, 23, 20, 15, 9, 4, 1]),
    ([5, 1], [1, 5, 13, 24, 35, 42, 42, 35, 24, 13, 5, 1]),
    ([6], [1, 4, 9, 15, 20, 23, 23, 20, 15, 9, 4, 1]),
    ])

hess_right[(0, 0, 1, 1, 1, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 10, 33, 83, 146, 174, 146, 83, 33, 10, 1]),
    ([2, 1, 1, 1, 1], [1, 6, 11, 21, 16, 0, -16, -21, -11, -6, -1]),
    ([2, 2, 1, 1], [1, 2, 1, 3, -6, -2, -6, 3, 1, 2, 1]),
    ([2, 2, 2], [1, -2, 3, -3, 0, 0, 0, 3, -3, 2, -1]),
    ([3, 1, 1, 1], [1, 4, 3, 2, -4, -12, -4, 2, 3, 4, 1]),
    ([3, 2, 1], [1, 0, -1, 0, -2, 0, 2, 0, 1, 0, -1]),
    ([3, 3], [1, -2, 0, 2, -1, 0, -1, 2, 0, -2, 1]),
    ([4, 1, 1], [1, 2, -1, -3, -4, 0, 4, 3, 1, -2, -1]),
    ([4, 2], [1, -2, 1, -1, 2, -2, 2, -1, 1, -2, 1]),
    ([5, 1], [1, 0, -2, -2, 1, 4, 1, -2, -2, 0, 1]),
    ([6], [1, -2, 0, 0, 3, 0, -3, 0, 0, 2, -1]),
    ])

hess_left[(0, 0, 1, 1, 1, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 10, 33, 83, 146, 174, 146, 83, 33, 10, 1]),
    ([2, 1, 1, 1, 1], [1, 8, 25, 57, 94, 110, 94, 57, 25, 8, 1]),
    ([2, 2, 1, 1], [1, 6, 17, 35, 54, 62, 54, 35, 17, 6, 1]),
    ([2, 2, 2], [1, 4, 9, 17, 26, 30, 26, 17, 9, 4, 1]),
    ([3, 1, 1, 1], [1, 7, 21, 44, 68, 78, 68, 44, 21, 7, 1]),
    ([3, 2, 1], [1, 5, 13, 24, 34, 38, 34, 24, 13, 5, 1]),
    ([3, 3], [1, 4, 9, 14, 17, 18, 17, 14, 9, 4, 1]),
    ([4, 1, 1], [1, 6, 17, 33, 48, 54, 48, 33, 17, 6, 1]),
    ([4, 2], [1, 4, 9, 15, 20, 22, 20, 15, 9, 4, 1]),
    ([5, 1], [1, 5, 13, 23, 31, 34, 31, 23, 13, 5, 1]),
    ([6], [1, 4, 9, 14, 17, 18, 17, 14, 9, 4, 1]),
    ])

hess_right[(0, 0, 1, 1, 1, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 10, 38, 116, 195, 195, 116, 38, 10, 1]),
    ([2, 1, 1, 1, 1], [1, 6, 14, 28, 11, -11, -28, -14, -6, -1]),
    ([2, 2, 1, 1], [1, 2, 2, 0, -5, -5, 0, 2, 2, 1]),
    ([2, 2, 2], [1, -2, 2, 0, -5, 5, 0, -2, 2, -1]),
    ([3, 1, 1, 1], [1, 4, 5, 2, -12, -12, 2, 5, 4, 1]),
    ([3, 2, 1], [1, 0, -1, -2, 2, -2, 2, 1, 0, -1]),
    ([3, 3], [1, -2, -1, 5, -3, -3, 5, -1, -2, 1]),
    ([4, 1, 1], [1, 2, 0, -6, -5, 5, 6, 0, -2, -1]),
    ([4, 2], [1, -2, 0, 2, -1, -1, 2, 0, -2, 1]),
    ([5, 1], [1, 0, -2, -4, 5, 5, -4, -2, 0, 1]),
    ([6], [1, -2, -1, 3, 1, -1, -3, 1, 2, -1]),
    ])

hess_left[(0, 0, 1, 1, 1, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 10, 38, 116, 195, 195, 116, 38, 10, 1]),
    ([2, 1, 1, 1, 1], [1, 8, 28, 70, 109, 109, 70, 28, 8, 1]),
    ([2, 2, 1, 1], [1, 6, 18, 36, 51, 51, 36, 18, 6, 1]),
    ([2, 2, 2], [1, 4, 8, 14, 21, 21, 14, 8, 4, 1]),
    ([3, 1, 1, 1], [1, 7, 23, 50, 72, 72, 50, 23, 7, 1]),
    ([3, 2, 1], [1, 5, 13, 22, 28, 28, 22, 13, 5, 1]),
    ([3, 3], [1, 4, 8, 11, 12, 12, 11, 8, 4, 1]),
    ([4, 1, 1], [1, 6, 18, 34, 45, 45, 34, 18, 6, 1]),
    ([4, 2], [1, 4, 8, 12, 15, 15, 12, 8, 4, 1]),
    ([5, 1], [1, 5, 13, 21, 25, 25, 21, 13, 5, 1]),
    ([6], [1, 4, 8, 11, 12, 12, 11, 8, 4, 1]),
    ])

hess_right[(0, 0, 1, 1, 1, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 15, 71, 165, 216, 165, 71, 15, 1]),
    ([2, 1, 1, 1, 1], [1, 9, 21, 23, 0, -23, -21, -9, -1]),
    ([2, 2, 1, 1], [1, 3, -1, 1, -8, 1, -1, 3, 1]),
    ([2, 2, 2], [1, -3, 5, -5, 0, 5, -5, 3, -1]),
    ([3, 1, 1, 1], [1, 6, 5, -6, -12, -6, 5, 6, 1]),
    ([3, 2, 1], [1, 0, -3, 2, 0, -2, 3, 0, -1]),
    ([3, 3], [1, -3, 2, 3, -6, 3, 2, -3, 1]),
    ([4, 1, 1], [1, 3, -3, -7, 0, 7, 3, -3, -1]),
    ([4, 2], [1, -3, 3, -1, 0, -1, 3, -3, 1]),
    ([5, 1], [1, 0, -4, 0, 6, 0, -4, 0, 1]),
    ([6], [1, -3, 2, 1, 0, -1, -2, 3, -1]),
    ])

hess_left[(0, 0, 1, 1, 1, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 15, 71, 165, 216, 165, 71, 15, 1]),
    ([2, 1, 1, 1, 1], [1, 11, 41, 85, 108, 85, 41, 11, 1]),
    ([2, 2, 1, 1], [1, 7, 19, 33, 40, 33, 19, 7, 1]),
    ([2, 2, 2], [1, 3, 5, 9, 12, 9, 5, 3, 1]),
    ([3, 1, 1, 1], [1, 9, 29, 54, 66, 54, 29, 9, 1]),
    ([3, 2, 1], [1, 5, 11, 16, 18, 16, 11, 5, 1]),
    ([3, 3], [1, 3, 5, 6, 6, 6, 5, 3, 1]),
    ([4, 1, 1], [1, 7, 19, 31, 36, 31, 19, 7, 1]),
    ([4, 2], [1, 3, 5, 7, 8, 7, 5, 3, 1]),
    ([5, 1], [1, 5, 11, 15, 16, 15, 11, 5, 1]),
    ([6], [1, 3, 5, 6, 6, 6, 5, 3, 1]),
    ])

hess_right[(0, 0, 1, 1, 1, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 48, 120, 186, 186, 120, 48, 6]),
    ([2, 1, 1, 1, 1], [4, 16, 16, 12, -12, -16, -16, -4]),
    ([2, 2, 1, 1], [2, 0, 0, -2, -2, 0, 0, 2]),
    ([3, 1, 1, 1], [3, 6, -3, -6, -6, -3, 6, 3]),
    ([3, 2, 1], [1, -2, 1, 0, 0, -1, 2, -1]),
    ([4, 1, 1], [2, 0, -4, -2, 2, 4, 0, -2]),
    ([5, 1], [1, -2, 0, 1, 1, 0, -2, 1]),
    ])

hess_left[(0, 0, 1, 1, 1, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 48, 120, 186, 186, 120, 48, 6]),
    ([2, 1, 1, 1, 1], [4, 24, 56, 84, 84, 56, 24, 4]),
    ([2, 2, 1, 1], [2, 8, 16, 22, 22, 16, 8, 2]),
    ([3, 1, 1, 1], [3, 15, 33, 48, 48, 33, 15, 3]),
    ([3, 2, 1], [1, 3, 5, 6, 6, 5, 3, 1]),
    ([4, 1, 1], [2, 8, 16, 22, 22, 16, 8, 2]),
    ([5, 1], [1, 3, 5, 6, 6, 5, 3, 1]),
    ])

hess_right[(0, 0, 1, 1, 2, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 10, 47, 119, 183, 183, 119, 47, 10, 1]),
    ([2, 1, 1, 1, 1], [1, 6, 17, 23, 11, -11, -23, -17, -6, -1]),
    ([2, 2, 1, 1], [1, 2, 3, -1, -5, -5, -1, 3, 2, 1]),
    ([2, 2, 2], [1, -2, 5, -9, 7, -7, 9, -5, 2, -1]),
    ([3, 1, 1, 1], [1, 4, 5, -1, -9, -9, -1, 5, 4, 1]),
    ([3, 2, 1], [1, 0, -1, -1, -1, 1, 1, 1, 0, -1]),
    ([3, 3], [1, -2, -1, 5, -3, -3, 5, -1, -2, 1]),
    ([4, 1, 1], [1, 2, -1, -5, -3, 3, 5, 1, -2, -1]),
    ([4, 2], [1, -2, 1, -1, 1, 1, -1, 1, -2, 1]),
    ([5, 1], [1, 0, -3, -1, 3, 3, -1, -3, 0, 1]),
    ([6], [1, -2, -1, 3, 1, -1, -3, 1, 2, -1]),
    ])

hess_left[(0, 0, 1, 1, 2, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 10, 47, 119, 183, 183, 119, 47, 10, 1]),
    ([2, 1, 1, 1, 1], [1, 8, 31, 71, 105, 105, 71, 31, 8, 1]),
    ([2, 2, 1, 1], [1, 6, 19, 39, 55, 55, 39, 19, 6, 1]),
    ([2, 2, 2], [1, 4, 11, 23, 33, 33, 23, 11, 4, 1]),
    ([3, 1, 1, 1], [1, 7, 23, 47, 66, 66, 47, 23, 7, 1]),
    ([3, 2, 1], [1, 5, 13, 23, 30, 30, 23, 13, 5, 1]),
    ([3, 3], [1, 4, 8, 11, 12, 12, 11, 8, 4, 1]),
    ([4, 1, 1], [1, 6, 17, 31, 41, 41, 31, 17, 6, 1]),
    ([4, 2], [1, 4, 9, 15, 19, 19, 15, 9, 4, 1]),
    ([5, 1], [1, 5, 12, 19, 23, 23, 19, 12, 5, 1]),
    ([6], [1, 4, 8, 11, 12, 12, 11, 8, 4, 1]),
    ])

hess_right[(0, 0, 1, 1, 2, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 10, 52, 166, 262, 166, 52, 10, 1]),
    ([2, 1, 1, 1, 1], [1, 6, 20, 34, 0, -34, -20, -6, -1]),
    ([2, 2, 1, 1], [1, 2, 4, -2, -10, -2, 4, 2, 1]),
    ([2, 2, 2], [1, -2, 4, -6, 0, 6, -4, 2, -1]),
    ([3, 1, 1, 1], [1, 4, 7, -2, -20, -2, 7, 4, 1]),
    ([3, 2, 1], [1, 0, -1, -2, 0, 2, 1, 0, -1]),
    ([3, 3], [1, -2, -2, 10, -14, 10, -2, -2, 1]),
    ([4, 1, 1], [1, 2, 0, -10, 0, 10, 0, -2, -1]),
    ([4, 2], [1, -2, 0, 2, -2, 2, 0, -2, 1]),
    ([5, 1], [1, 0, -3, -4, 12, -4, -3, 0, 1]),
    ([6], [1, -2, -2, 6, 0, -6, 2, 2, -1]),
    ])

hess_left[(0, 0, 1, 1, 2, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 10, 52, 166, 262, 166, 52, 10, 1]),
    ([2, 1, 1, 1, 1], [1, 8, 34, 88, 122, 88, 34, 8, 1]),
    ([2, 2, 1, 1], [1, 6, 20, 42, 54, 42, 20, 6, 1]),
    ([2, 2, 2], [1, 4, 10, 20, 26, 20, 10, 4, 1]),
    ([3, 1, 1, 1], [1, 7, 25, 52, 64, 52, 25, 7, 1]),
    ([3, 2, 1], [1, 5, 13, 22, 26, 22, 13, 5, 1]),
    ([3, 3], [1, 4, 7, 10, 10, 10, 7, 4, 1]),
    ([4, 1, 1], [1, 6, 18, 30, 34, 30, 18, 6, 1]),
    ([4, 2], [1, 4, 8, 12, 14, 12, 8, 4, 1]),
    ([5, 1], [1, 5, 12, 16, 17, 16, 12, 5, 1]),
    ([6], [1, 4, 7, 8, 8, 8, 7, 4, 1]),
    ])

hess_right[(0, 0, 1, 1, 2, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 15, 99, 245, 245, 99, 15, 1]),
    ([2, 1, 1, 1, 1], [1, 9, 31, 23, -23, -31, -9, -1]),
    ([2, 2, 1, 1], [1, 3, 3, -7, -7, 3, 3, 1]),
    ([2, 2, 2], [1, -3, 7, -13, 13, -7, 3, -1]),
    ([3, 1, 1, 1], [1, 6, 6, -13, -13, 6, 6, 1]),
    ([3, 2, 1], [1, 0, -2, -1, 1, 2, 0, -1]),
    ([3, 3], [1, -3, 3, -1, -1, 3, -3, 1]),
    ([4, 1, 1], [1, 3, -5, -7, 7, 5, -3, -1]),
    ([4, 2], [1, -3, 3, -1, -1, 3, -3, 1]),
    ([5, 1], [1, 0, -6, 5, 5, -6, 0, 1]),
    ([6], [1, -3, 1, 5, -5, -1, 3, -1]),
    ])

hess_left[(0, 0, 1, 1, 2, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 15, 99, 245, 245, 99, 15, 1]),
    ([2, 1, 1, 1, 1], [1, 11, 51, 105, 105, 51, 11, 1]),
    ([2, 2, 1, 1], [1, 7, 23, 41, 41, 23, 7, 1]),
    ([2, 2, 2], [1, 3, 7, 13, 13, 7, 3, 1]),
    ([3, 1, 1, 1], [1, 9, 30, 50, 50, 30, 9, 1]),
    ([3, 2, 1], [1, 5, 12, 18, 18, 12, 5, 1]),
    ([3, 3], [1, 3, 6, 8, 8, 6, 3, 1]),
    ([4, 1, 1], [1, 7, 17, 23, 23, 17, 7, 1]),
    ([4, 2], [1, 3, 5, 7, 7, 5, 3, 1]),
    ([5, 1], [1, 5, 9, 10, 10, 9, 5, 1]),
    ([6], [1, 3, 4, 4, 4, 4, 3, 1]),
    ])

hess_right[(0, 0, 1, 1, 2, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 48, 174, 264, 174, 48, 6]),
    ([2, 1, 1, 1, 1], [4, 16, 28, 0, -28, -16, -4]),
    ([2, 2, 1, 1], [2, 0, 2, -8, 2, 0, 2]),
    ([3, 1, 1, 1], [3, 6, -3, -12, -3, 6, 3]),
    ([3, 2, 1], [1, -2, 1, 0, -1, 2, -1]),
    ([4, 1, 1], [2, 0, -6, 0, 6, 0, -2]),
    ([5, 1], [1, -2, -1, 4, -1, -2, 1]),
    ])

hess_left[(0, 0, 1, 1, 2, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 48, 174, 264, 174, 48, 6]),
    ([2, 1, 1, 1, 1], [4, 24, 68, 96, 68, 24, 4]),
    ([2, 2, 1, 1], [2, 8, 18, 24, 18, 8, 2]),
    ([3, 1, 1, 1], [3, 15, 33, 42, 33, 15, 3]),
    ([3, 2, 1], [1, 3, 5, 6, 5, 3, 1]),
    ([4, 1, 1], [2, 8, 14, 16, 14, 8, 2]),
    ([5, 1], [1, 3, 4, 4, 4, 3, 1]),
    ])

hess_right[(0, 0, 1, 1, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 15, 99, 245, 245, 99, 15, 1]),
    ([2, 1, 1, 1, 1], [1, 9, 31, 23, -23, -31, -9, -1]),
    ([2, 2, 1, 1], [1, 3, 3, -7, -7, 3, 3, 1]),
    ([2, 2, 2], [1, -3, 7, -13, 13, -7, 3, -1]),
    ([3, 1, 1, 1], [1, 6, 6, -13, -13, 6, 6, 1]),
    ([3, 2, 1], [1, 0, -2, -1, 1, 2, 0, -1]),
    ([3, 3], [1, -3, 3, -1, -1, 3, -3, 1]),
    ([4, 1, 1], [1, 3, -5, -7, 7, 5, -3, -1]),
    ([4, 2], [1, -3, 3, -1, -1, 3, -3, 1]),
    ([5, 1], [1, 0, -6, 5, 5, -6, 0, 1]),
    ([6], [1, -3, 1, 5, -5, -1, 3, -1]),
    ])

hess_left[(0, 0, 1, 1, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 15, 99, 245, 245, 99, 15, 1]),
    ([2, 1, 1, 1, 1], [1, 11, 51, 105, 105, 51, 11, 1]),
    ([2, 2, 1, 1], [1, 7, 23, 41, 41, 23, 7, 1]),
    ([2, 2, 2], [1, 3, 7, 13, 13, 7, 3, 1]),
    ([3, 1, 1, 1], [1, 9, 30, 50, 50, 30, 9, 1]),
    ([3, 2, 1], [1, 5, 12, 18, 18, 12, 5, 1]),
    ([3, 3], [1, 3, 6, 8, 8, 6, 3, 1]),
    ([4, 1, 1], [1, 7, 17, 23, 23, 17, 7, 1]),
    ([4, 2], [1, 3, 5, 7, 7, 5, 3, 1]),
    ([5, 1], [1, 5, 9, 10, 10, 9, 5, 1]),
    ([6], [1, 3, 4, 4, 4, 4, 3, 1]),
    ])

hess_right[(0, 0, 1, 1, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 29, 175, 310, 175, 29, 1]),
    ([2, 1, 1, 1, 1], [1, 15, 39, 0, -39, -15, -1]),
    ([2, 2, 1, 1], [1, 5, -1, -10, -1, 5, 1]),
    ([2, 2, 2], [1, -1, -1, 0, 1, 1, -1]),
    ([3, 1, 1, 1], [1, 8, 1, -20, 1, 8, 1]),
    ([3, 2, 1], [1, 0, -3, 0, 3, 0, -1]),
    ([3, 3], [1, -4, 7, -8, 7, -4, 1]),
    ([4, 1, 1], [1, 3, -9, 0, 9, -3, -1]),
    ([4, 2], [1, -3, 3, -2, 3, -3, 1]),
    ([5, 1], [1, -1, -5, 10, -5, -1, 1]),
    ([6], [1, -4, 5, 0, -5, 4, -1]),
    ])

hess_left[(0, 0, 1, 1, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 29, 175, 310, 175, 29, 1]),
    ([2, 1, 1, 1, 1], [1, 17, 71, 110, 71, 17, 1]),
    ([2, 2, 1, 1], [1, 9, 27, 38, 27, 9, 1]),
    ([2, 2, 2], [1, 5, 11, 14, 11, 5, 1]),
    ([3, 1, 1, 1], [1, 11, 31, 40, 31, 11, 1]),
    ([3, 2, 1], [1, 5, 11, 14, 11, 5, 1]),
    ([3, 3], [1, 2, 4, 4, 4, 2, 1]),
    ([4, 1, 1], [1, 7, 13, 14, 13, 7, 1]),
    ([4, 2], [1, 3, 5, 6, 5, 3, 1]),
    ([5, 1], [1, 4, 5, 5, 5, 4, 1]),
    ([6], [1, 2, 2, 2, 2, 2, 1]),
    ])

hess_right[(0, 0, 1, 1, 3, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 72, 282, 282, 72, 6]),
    ([2, 1, 1, 1, 1], [4, 24, 28, -28, -24, -4]),
    ([2, 2, 1, 1], [2, 0, -2, -2, 0, 2]),
    ([3, 1, 1, 1], [3, 9, -12, -12, 9, 3]),
    ([3, 2, 1], [1, -3, 4, -4, 3, -1]),
    ([4, 1, 1], [2, 0, -10, 10, 0, -2]),
    ([5, 1], [1, -3, 2, 2, -3, 1]),
    ])

hess_left[(0, 0, 1, 1, 3, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 72, 282, 282, 72, 6]),
    ([2, 1, 1, 1, 1], [4, 32, 84, 84, 32, 4]),
    ([2, 2, 1, 1], [2, 8, 14, 14, 8, 2]),
    ([3, 1, 1, 1], [3, 18, 33, 33, 18, 3]),
    ([3, 2, 1], [1, 2, 3, 3, 2, 1]),
    ([4, 1, 1], [2, 8, 10, 10, 8, 2]),
    ([5, 1], [1, 2, 2, 2, 2, 1]),
    ])

hess_right[(0, 0, 1, 1, 4, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [15, 105, 240, 240, 105, 15]),
    ([2, 1, 1, 1, 1], [7, 23, 16, -16, -23, -7]),
    ([2, 2, 1, 1], [3, 1, -4, -4, 1, 3]),
    ([2, 2, 2], [3, -9, 12, -12, 9, -3]),
    ([3, 1, 1, 1], [3, 3, -6, -6, 3, 3]),
    ([3, 2, 1], [1, -1, -2, 2, 1, -1]),
    ([4, 1, 1], [1, -1, -2, 2, 1, -1]),
    ([4, 2], [1, -3, 2, 2, -3, 1]),
    ])

hess_left[(0, 0, 1, 1, 4, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [15, 105, 240, 240, 105, 15]),
    ([2, 1, 1, 1, 1], [7, 37, 76, 76, 37, 7]),
    ([2, 2, 1, 1], [3, 13, 24, 24, 13, 3]),
    ([2, 2, 2], [3, 9, 12, 12, 9, 3]),
    ([3, 1, 1, 1], [3, 12, 21, 21, 12, 3]),
    ([3, 2, 1], [1, 4, 7, 7, 4, 1]),
    ([4, 1, 1], [1, 3, 4, 4, 3, 1]),
    ([4, 2], [1, 3, 4, 4, 3, 1]),
    ])

hess_right[(0, 0, 1, 1, 4, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [30, 180, 300, 180, 30]),
    ([2, 1, 1, 1, 1], [12, 24, 0, -24, -12]),
    ([2, 2, 1, 1], [2, -4, 4, -4, 2]),
    ([3, 1, 1, 1], [6, 0, -12, 0, 6]),
    ([4, 1, 1], [2, -4, 0, 4, -2]),
    ])

hess_left[(0, 0, 1, 1, 4, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [30, 180, 300, 180, 30]),
    ([2, 1, 1, 1, 1], [12, 48, 72, 48, 12]),
    ([2, 2, 1, 1], [2, 4, 4, 4, 2]),
    ([3, 1, 1, 1], [6, 18, 24, 18, 6]),
    ([4, 1, 1], [2, 4, 4, 4, 2]),
    ])

hess_right[(0, 0, 1, 2, 2, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 24, 83, 156, 192, 156, 83, 24, 1]),
    ([2, 1, 1, 1, 1], [1, 12, 19, 18, 0, -18, -19, -12, -1]),
    ([2, 2, 1, 1], [1, 4, -1, 0, -8, 0, -1, 4, 1]),
    ([2, 2, 2], [1, 0, -1, -2, 0, 2, 1, 0, -1]),
    ([3, 1, 1, 1], [1, 6, 2, -6, -6, -6, 2, 6, 1]),
    ([3, 2, 1], [1, 0, -2, 0, 0, 0, 2, 0, -1]),
    ([3, 3], [1, -3, 2, 3, -6, 3, 2, -3, 1]),
    ([4, 1, 1], [1, 2, -3, -4, 0, 4, 3, -2, -1]),
    ([4, 2], [1, -2, 1, -2, 4, -2, 1, -2, 1]),
    ([5, 1], [1, -1, -2, 1, 2, 1, -2, -1, 1]),
    ([6], [1, -3, 2, 1, 0, -1, -2, 3, -1]),
    ])

hess_left[(0, 0, 1, 2, 2, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 24, 83, 156, 192, 156, 83, 24, 1]),
    ([2, 1, 1, 1, 1], [1, 14, 45, 82, 100, 82, 45, 14, 1]),
    ([2, 2, 1, 1], [1, 8, 23, 40, 48, 40, 23, 8, 1]),
    ([2, 2, 2], [1, 6, 17, 30, 36, 30, 17, 6, 1]),
    ([3, 1, 1, 1], [1, 9, 26, 45, 54, 45, 26, 9, 1]),
    ([3, 2, 1], [1, 5, 12, 19, 22, 19, 12, 5, 1]),
    ([3, 3], [1, 3, 5, 6, 6, 6, 5, 3, 1]),
    ([4, 1, 1], [1, 6, 15, 24, 28, 24, 15, 6, 1]),
    ([4, 2], [1, 4, 9, 14, 16, 14, 9, 4, 1]),
    ([5, 1], [1, 4, 8, 11, 12, 11, 8, 4, 1]),
    ([6], [1, 3, 5, 6, 6, 6, 5, 3, 1]),
    ])

hess_right[(0, 0, 1, 2, 2, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 24, 102, 233, 233, 102, 24, 1]),
    ([2, 1, 1, 1, 1], [1, 12, 26, 23, -23, -26, -12, -1]),
    ([2, 2, 1, 1], [1, 4, 2, -7, -7, 2, 4, 1]),
    ([2, 2, 2], [1, 0, -2, -1, 1, 2, 0, -1]),
    ([3, 1, 1, 1], [1, 6, 3, -10, -10, 3, 6, 1]),
    ([3, 2, 1], [1, 0, -1, -4, 4, 1, 0, -1]),
    ([3, 3], [1, -3, 3, -1, -1, 3, -3, 1]),
    ([4, 1, 1], [1, 2, -4, -5, 5, 4, -2, -1]),
    ([4, 2], [1, -2, 0, 1, 1, 0, -2, 1]),
    ([5, 1], [1, -1, -3, 3, 3, -3, -1, 1]),
    ([6], [1, -3, 1, 5, -5, -1, 3, -1]),
    ])

hess_left[(0, 0, 1, 2, 2, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 24, 102, 233, 233, 102, 24, 1]),
    ([2, 1, 1, 1, 1], [1, 14, 52, 101, 101, 52, 14, 1]),
    ([2, 2, 1, 1], [1, 8, 26, 45, 45, 26, 8, 1]),
    ([2, 2, 2], [1, 6, 16, 25, 25, 16, 6, 1]),
    ([3, 1, 1, 1], [1, 9, 27, 44, 44, 27, 9, 1]),
    ([3, 2, 1], [1, 5, 13, 20, 20, 13, 5, 1]),
    ([3, 3], [1, 3, 6, 8, 8, 6, 3, 1]),
    ([4, 1, 1], [1, 6, 14, 19, 19, 14, 6, 1]),
    ([4, 2], [1, 4, 8, 11, 11, 8, 4, 1]),
    ([5, 1], [1, 4, 7, 8, 8, 7, 4, 1]),
    ([6], [1, 3, 4, 4, 4, 4, 3, 1]),
    ])

hess_right[(0, 0, 1, 2, 2, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 29, 175, 310, 175, 29, 1]),
    ([2, 1, 1, 1, 1], [1, 15, 39, 0, -39, -15, -1]),
    ([2, 2, 1, 1], [1, 5, -1, -10, -1, 5, 1]),
    ([2, 2, 2], [1, -1, -1, 0, 1, 1, -1]),
    ([3, 1, 1, 1], [1, 8, 1, -20, 1, 8, 1]),
    ([3, 2, 1], [1, 0, -3, 0, 3, 0, -1]),
    ([3, 3], [1, -4, 7, -8, 7, -4, 1]),
    ([4, 1, 1], [1, 3, -9, 0, 9, -3, -1]),
    ([4, 2], [1, -3, 3, -2, 3, -3, 1]),
    ([5, 1], [1, -1, -5, 10, -5, -1, 1]),
    ([6], [1, -4, 5, 0, -5, 4, -1]),
    ])

hess_left[(0, 0, 1, 2, 2, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 29, 175, 310, 175, 29, 1]),
    ([2, 1, 1, 1, 1], [1, 17, 71, 110, 71, 17, 1]),
    ([2, 2, 1, 1], [1, 9, 27, 38, 27, 9, 1]),
    ([2, 2, 2], [1, 5, 11, 14, 11, 5, 1]),
    ([3, 1, 1, 1], [1, 11, 31, 40, 31, 11, 1]),
    ([3, 2, 1], [1, 5, 11, 14, 11, 5, 1]),
    ([3, 3], [1, 2, 4, 4, 4, 2, 1]),
    ([4, 1, 1], [1, 7, 13, 14, 13, 7, 1]),
    ([4, 2], [1, 3, 5, 6, 5, 3, 1]),
    ([5, 1], [1, 4, 5, 5, 5, 4, 1]),
    ([6], [1, 2, 2, 2, 2, 2, 1]),
    ])

hess_right[(0, 0, 1, 2, 2, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 102, 252, 252, 102, 6]),
    ([2, 1, 1, 1, 1], [4, 28, 16, -16, -28, -4]),
    ([2, 2, 1, 1], [2, 2, -4, -4, 2, 2]),
    ([3, 1, 1, 1], [3, 6, -9, -9, 6, 3]),
    ([3, 2, 1], [1, -2, 1, -1, 2, -1]),
    ([4, 1, 1], [2, -2, -4, 4, 2, -2]),
    ([5, 1], [1, -3, 2, 2, -3, 1]),
    ])

hess_left[(0, 0, 1, 2, 2, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 102, 252, 252, 102, 6]),
    ([2, 1, 1, 1, 1], [4, 36, 80, 80, 36, 4]),
    ([2, 2, 1, 1], [2, 10, 20, 20, 10, 2]),
    ([3, 1, 1, 1], [3, 15, 27, 27, 15, 3]),
    ([3, 2, 1], [1, 3, 5, 5, 3, 1]),
    ([4, 1, 1], [2, 6, 8, 8, 6, 2]),
    ([5, 1], [1, 2, 2, 2, 2, 1]),
    ])

hess_right[(0, 0, 1, 2, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 43, 179, 274, 179, 43, 1]),
    ([2, 1, 1, 1, 1], [1, 19, 31, 0, -31, -19, -1]),
    ([2, 2, 1, 1], [1, 7, -5, -6, -5, 7, 1]),
    ([2, 2, 2], [1, -1, -1, 0, 1, 1, -1]),
    ([3, 1, 1, 1], [1, 7, -1, -14, -1, 7, 1]),
    ([3, 2, 1], [1, 1, -5, 0, 5, -1, -1]),
    ([3, 3], [1, -2, -1, 4, -1, -2, 1]),
    ([4, 1, 1], [1, 1, -5, 0, 5, -1, -1]),
    ([4, 2], [1, -3, 3, -2, 3, -3, 1]),
    ([5, 1], [1, -2, -1, 4, -1, -2, 1]),
    ([6], [1, -4, 5, 0, -5, 4, -1]),
    ])

hess_left[(0, 0, 1, 2, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 43, 179, 274, 179, 43, 1]),
    ([2, 1, 1, 1, 1], [1, 21, 71, 102, 71, 21, 1]),
    ([2, 2, 1, 1], [1, 11, 31, 42, 31, 11, 1]),
    ([2, 2, 2], [1, 5, 11, 14, 11, 5, 1]),
    ([3, 1, 1, 1], [1, 10, 26, 34, 26, 10, 1]),
    ([3, 2, 1], [1, 6, 14, 18, 14, 6, 1]),
    ([3, 3], [1, 4, 8, 10, 8, 4, 1]),
    ([4, 1, 1], [1, 5, 9, 10, 9, 5, 1]),
    ([4, 2], [1, 3, 5, 6, 5, 3, 1]),
    ([5, 1], [1, 3, 4, 4, 4, 3, 1]),
    ([6], [1, 2, 2, 2, 2, 2, 1]),
    ])

hess_right[(0, 0, 1, 2, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 57, 302, 302, 57, 1]),
    ([2, 1, 1, 1, 1], [1, 25, 40, -40, -25, -1]),
    ([2, 2, 1, 1], [1, 9, -10, -10, 9, 1]),
    ([2, 2, 2], [1, 1, -8, 8, -1, -1]),
    ([3, 1, 1, 1], [1, 9, -10, -10, 9, 1]),
    ([3, 2, 1], [1, 1, -8, 8, -1, -1]),
    ([3, 3], [1, -3, 2, 2, -3, 1]),
    ([4, 1, 1], [1, 1, -8, 8, -1, -1]),
    ([4, 2], [1, -3, 2, 2, -3, 1]),
    ([5, 1], [1, -3, 2, 2, -3, 1]),
    ([6], [1, -5, 10, -10, 5, -1]),
    ])

hess_left[(0, 0, 1, 2, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [1, 57, 302, 302, 57, 1]),
    ([2, 1, 1, 1, 1], [1, 27, 92, 92, 27, 1]),
    ([2, 2, 1, 1], [1, 13, 34, 34, 13, 1]),
    ([2, 2, 2], [1, 7, 16, 16, 7, 1]),
    ([3, 1, 1, 1], [1, 12, 23, 23, 12, 1]),
    ([3, 2, 1], [1, 6, 11, 11, 6, 1]),
    ([3, 3], [1, 3, 5, 5, 3, 1]),
    ([4, 1, 1], [1, 5, 6, 6, 5, 1]),
    ([4, 2], [1, 3, 4, 4, 3, 1]),
    ([5, 1], [1, 2, 2, 2, 2, 1]),
    ([6], [1, 1, 1, 1, 1, 1]),
    ])

hess_right[(0, 0, 1, 2, 3, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 156, 396, 156, 6]),
    ([2, 1, 1, 1, 1], [4, 40, 0, -40, -4]),
    ([2, 2, 1, 1], [2, 4, -12, 4, 2]),
    ([3, 1, 1, 1], [3, 6, -18, 6, 3]),
    ([3, 2, 1], [1, -2, 0, 2, -1]),
    ([4, 1, 1], [2, -4, 0, 4, -2]),
    ([5, 1], [1, -4, 6, -4, 1]),
    ])

hess_left[(0, 0, 1, 2, 3, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 156, 396, 156, 6]),
    ([2, 1, 1, 1, 1], [4, 48, 88, 48, 4]),
    ([2, 2, 1, 1], [2, 12, 20, 12, 2]),
    ([3, 1, 1, 1], [3, 15, 18, 15, 3]),
    ([3, 2, 1], [1, 3, 4, 3, 1]),
    ([4, 1, 1], [2, 4, 4, 4, 2]),
    ([5, 1], [1, 1, 1, 1, 1]),
    ])

hess_right[(0, 0, 1, 2, 4, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [15, 180, 330, 180, 15]),
    ([2, 1, 1, 1, 1], [7, 34, 0, -34, -7]),
    ([2, 2, 1, 1], [3, 4, -14, 4, 3]),
    ([2, 2, 2], [3, -6, 0, 6, -3]),
    ([3, 1, 1, 1], [3, 0, -6, 0, 3]),
    ([3, 2, 1], [1, -2, 0, 2, -1]),
    ([4, 1, 1], [1, -2, 0, 2, -1]),
    ([4, 2], [1, -4, 6, -4, 1]),
    ])

hess_left[(0, 0, 1, 2, 4, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [15, 180, 330, 180, 15]),
    ([2, 1, 1, 1, 1], [7, 48, 82, 48, 7]),
    ([2, 2, 1, 1], [3, 16, 26, 16, 3]),
    ([2, 2, 2], [3, 12, 18, 12, 3]),
    ([3, 1, 1, 1], [3, 9, 12, 9, 3]),
    ([3, 2, 1], [1, 3, 4, 3, 1]),
    ([4, 1, 1], [1, 2, 2, 2, 1]),
    ([4, 2], [1, 2, 2, 2, 1]),
    ])

hess_right[(0, 0, 1, 2, 4, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [30, 330, 330, 30]),
    ([2, 1, 1, 1, 1], [12, 36, -36, -12]),
    ([2, 2, 1, 1], [2, -2, -2, 2]),
    ([3, 1, 1, 1], [6, -6, -6, 6]),
    ([4, 1, 1], [2, -6, 6, -2]),
    ])

hess_left[(0, 0, 1, 2, 4, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [30, 330, 330, 30]),
    ([2, 1, 1, 1, 1], [12, 60, 60, 12]),
    ([2, 2, 1, 1], [2, 6, 6, 2]),
    ([3, 1, 1, 1], [6, 12, 12, 6]),
    ([4, 1, 1], [2, 2, 2, 2]),
    ])

hess_right[(0, 0, 1, 3, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [20, 120, 220, 220, 120, 20]),
    ([2, 1, 1, 1, 1], [8, 24, 8, -8, -24, -8]),
    ([2, 2, 1, 1], [4, 0, -4, -4, 0, 4]),
    ([3, 1, 1, 1], [2, 3, -5, -5, 3, 2]),
    ([3, 2, 1], [2, -3, -1, 1, 3, -2]),
    ([3, 3], [2, -6, 4, 4, -6, 2]),
    ])

hess_left[(0, 0, 1, 3, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [20, 120, 220, 220, 120, 20]),
    ([2, 1, 1, 1, 1], [8, 40, 72, 72, 40, 8]),
    ([2, 2, 1, 1], [4, 16, 28, 28, 16, 4]),
    ([3, 1, 1, 1], [2, 9, 16, 16, 9, 2]),
    ([3, 2, 1], [2, 7, 12, 12, 7, 2]),
    ([3, 3], [2, 6, 10, 10, 6, 2]),
    ])

hess_right[(0, 0, 1, 3, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [20, 160, 360, 160, 20]),
    ([2, 1, 1, 1, 1], [8, 32, 0, -32, -8]),
    ([2, 2, 1, 1], [4, 0, -8, 0, 4]),
    ([3, 1, 1, 1], [2, 4, -12, 4, 2]),
    ([3, 2, 1], [2, -4, 0, 4, -2]),
    ([3, 3], [2, -8, 12, -8, 2]),
    ])

hess_left[(0, 0, 1, 3, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [20, 160, 360, 160, 20]),
    ([2, 1, 1, 1, 1], [8, 48, 80, 48, 8]),
    ([2, 2, 1, 1], [4, 16, 24, 16, 4]),
    ([3, 1, 1, 1], [2, 10, 12, 10, 2]),
    ([3, 2, 1], [2, 6, 8, 6, 2]),
    ([3, 3], [2, 4, 6, 4, 2]),
    ])

hess_right[(0, 0, 1, 3, 3, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [60, 300, 300, 60]),
    ([2, 1, 1, 1, 1], [16, 24, -24, -16]),
    ([2, 2, 1, 1], [4, -4, -4, 4]),
    ([3, 1, 1, 1], [3, -3, -3, 3]),
    ([3, 2, 1], [1, -3, 3, -1]),
    ])

hess_left[(0, 0, 1, 3, 3, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [60, 300, 300, 60]),
    ([2, 1, 1, 1, 1], [16, 56, 56, 16]),
    ([2, 2, 1, 1], [4, 12, 12, 4]),
    ([3, 1, 1, 1], [3, 6, 6, 3]),
    ([3, 2, 1], [1, 2, 2, 1]),
    ])

hess_right[(0, 0, 1, 3, 4, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [60, 300, 300, 60]),
    ([2, 1, 1, 1, 1], [16, 24, -24, -16]),
    ([2, 2, 1, 1], [4, -4, -4, 4]),
    ([3, 1, 1, 1], [3, -3, -3, 3]),
    ([3, 2, 1], [1, -3, 3, -1]),
    ])

hess_left[(0, 0, 1, 3, 4, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [60, 300, 300, 60]),
    ([2, 1, 1, 1, 1], [16, 56, 56, 16]),
    ([2, 2, 1, 1], [4, 12, 12, 4]),
    ([3, 1, 1, 1], [3, 6, 6, 3]),
    ([3, 2, 1], [1, 2, 2, 1]),
    ])

hess_right[(0, 0, 1, 3, 4, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [120, 480, 120]),
    ([2, 1, 1, 1, 1], [24, 0, -24]),
    ([3, 1, 1, 1], [6, -12, 6]),
    ])

hess_left[(0, 0, 1, 3, 4, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [120, 480, 120]),
    ([2, 1, 1, 1, 1], [24, 48, 24]),
    ([3, 1, 1, 1], [6, 6, 6]),
    ])

hess_right[(0, 0, 2, 2, 2, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [15, 60, 120, 165, 165, 120, 60, 15]),
    ([2, 1, 1, 1, 1], [7, 14, 14, 7, -7, -14, -14, -7]),
    ([2, 2, 1, 1], [3, 0, 0, -3, -3, 0, 0, 3]),
    ([2, 2, 2], [3, -6, 6, -9, 9, -6, 6, -3]),
    ([3, 1, 1, 1], [3, 3, -3, -3, -3, -3, 3, 3]),
    ([3, 2, 1], [1, -1, -1, 1, -1, 1, 1, -1]),
    ([4, 1, 1], [1, 0, -2, -1, 1, 2, 0, -1]),
    ([4, 2], [1, -2, 0, 1, 1, 0, -2, 1]),
    ])

hess_left[(0, 0, 2, 2, 2, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [15, 60, 120, 165, 165, 120, 60, 15]),
    ([2, 1, 1, 1, 1], [7, 28, 56, 77, 77, 56, 28, 7]),
    ([2, 2, 1, 1], [3, 12, 24, 33, 33, 24, 12, 3]),
    ([2, 2, 2], [3, 12, 24, 33, 33, 24, 12, 3]),
    ([3, 1, 1, 1], [3, 12, 24, 33, 33, 24, 12, 3]),
    ([3, 2, 1], [1, 4, 8, 11, 11, 8, 4, 1]),
    ([4, 1, 1], [1, 4, 8, 11, 11, 8, 4, 1]),
    ([4, 2], [1, 4, 8, 11, 11, 8, 4, 1]),
    ])

hess_right[(0, 0, 2, 2, 2, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [15, 60, 165, 240, 165, 60, 15]),
    ([2, 1, 1, 1, 1], [7, 14, 23, 0, -23, -14, -7]),
    ([2, 2, 1, 1], [3, 0, 1, -8, 1, 0, 3]),
    ([2, 2, 2], [3, -6, 3, 0, -3, 6, -3]),
    ([3, 1, 1, 1], [3, 3, -3, -6, -3, 3, 3]),
    ([3, 2, 1], [1, -1, -1, 0, 1, 1, -1]),
    ([4, 1, 1], [1, 0, -3, 0, 3, 0, -1]),
    ([4, 2], [1, -2, -1, 4, -1, -2, 1]),
    ])

hess_left[(0, 0, 2, 2, 2, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [15, 60, 165, 240, 165, 60, 15]),
    ([2, 1, 1, 1, 1], [7, 28, 65, 88, 65, 28, 7]),
    ([2, 2, 1, 1], [3, 12, 25, 32, 25, 12, 3]),
    ([2, 2, 2], [3, 12, 21, 24, 21, 12, 3]),
    ([3, 1, 1, 1], [3, 12, 24, 30, 24, 12, 3]),
    ([3, 2, 1], [1, 4, 8, 10, 8, 4, 1]),
    ([4, 1, 1], [1, 4, 7, 8, 7, 4, 1]),
    ([4, 2], [1, 4, 7, 8, 7, 4, 1]),
    ])

hess_right[(0, 0, 2, 2, 2, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [15, 105, 240, 240, 105, 15]),
    ([2, 1, 1, 1, 1], [7, 23, 16, -16, -23, -7]),
    ([2, 2, 1, 1], [3, 1, -4, -4, 1, 3]),
    ([2, 2, 2], [3, -9, 12, -12, 9, -3]),
    ([3, 1, 1, 1], [3, 3, -6, -6, 3, 3]),
    ([3, 2, 1], [1, -1, -2, 2, 1, -1]),
    ([4, 1, 1], [1, -1, -2, 2, 1, -1]),
    ([4, 2], [1, -3, 2, 2, -3, 1]),
    ])

hess_left[(0, 0, 2, 2, 2, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [15, 105, 240, 240, 105, 15]),
    ([2, 1, 1, 1, 1], [7, 37, 76, 76, 37, 7]),
    ([2, 2, 1, 1], [3, 13, 24, 24, 13, 3]),
    ([2, 2, 2], [3, 9, 12, 12, 9, 3]),
    ([3, 1, 1, 1], [3, 12, 21, 21, 12, 3]),
    ([3, 2, 1], [1, 4, 7, 7, 4, 1]),
    ([4, 1, 1], [1, 3, 4, 4, 3, 1]),
    ([4, 2], [1, 3, 4, 4, 3, 1]),
    ])

hess_right[(0, 0, 2, 2, 2, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [60, 180, 240, 180, 60]),
    ([2, 1, 1, 1, 1], [16, 16, 0, -16, -16]),
    ([2, 2, 1, 1], [4, -4, 0, -4, 4]),
    ([3, 1, 1, 1], [3, 0, -6, 0, 3]),
    ([3, 2, 1], [1, -2, 0, 2, -1]),
    ])

hess_left[(0, 0, 2, 2, 2, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [60, 180, 240, 180, 60]),
    ([2, 1, 1, 1, 1], [16, 48, 64, 48, 16]),
    ([2, 2, 1, 1], [4, 12, 16, 12, 4]),
    ([3, 1, 1, 1], [3, 9, 12, 9, 3]),
    ([3, 2, 1], [1, 3, 4, 3, 1]),
    ])

hess_right[(0, 0, 2, 2, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [15, 105, 240, 240, 105, 15]),
    ([2, 1, 1, 1, 1], [7, 23, 16, -16, -23, -7]),
    ([2, 2, 1, 1], [3, 1, -4, -4, 1, 3]),
    ([2, 2, 2], [3, -9, 12, -12, 9, -3]),
    ([3, 1, 1, 1], [3, 3, -6, -6, 3, 3]),
    ([3, 2, 1], [1, -1, -2, 2, 1, -1]),
    ([4, 1, 1], [1, -1, -2, 2, 1, -1]),
    ([4, 2], [1, -3, 2, 2, -3, 1]),
    ])

hess_left[(0, 0, 2, 2, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [15, 105, 240, 240, 105, 15]),
    ([2, 1, 1, 1, 1], [7, 37, 76, 76, 37, 7]),
    ([2, 2, 1, 1], [3, 13, 24, 24, 13, 3]),
    ([2, 2, 2], [3, 9, 12, 12, 9, 3]),
    ([3, 1, 1, 1], [3, 12, 21, 21, 12, 3]),
    ([3, 2, 1], [1, 4, 7, 7, 4, 1]),
    ([4, 1, 1], [1, 3, 4, 4, 3, 1]),
    ([4, 2], [1, 3, 4, 4, 3, 1]),
    ])

hess_right[(0, 0, 2, 2, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [15, 180, 330, 180, 15]),
    ([2, 1, 1, 1, 1], [7, 34, 0, -34, -7]),
    ([2, 2, 1, 1], [3, 4, -14, 4, 3]),
    ([2, 2, 2], [3, -6, 0, 6, -3]),
    ([3, 1, 1, 1], [3, 0, -6, 0, 3]),
    ([3, 2, 1], [1, -2, 0, 2, -1]),
    ([4, 1, 1], [1, -2, 0, 2, -1]),
    ([4, 2], [1, -4, 6, -4, 1]),
    ])

hess_left[(0, 0, 2, 2, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [15, 180, 330, 180, 15]),
    ([2, 1, 1, 1, 1], [7, 48, 82, 48, 7]),
    ([2, 2, 1, 1], [3, 16, 26, 16, 3]),
    ([2, 2, 2], [3, 12, 18, 12, 3]),
    ([3, 1, 1, 1], [3, 9, 12, 9, 3]),
    ([3, 2, 1], [1, 3, 4, 3, 1]),
    ([4, 1, 1], [1, 2, 2, 2, 1]),
    ([4, 2], [1, 2, 2, 2, 1]),
    ])

hess_right[(0, 0, 2, 2, 3, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [60, 300, 300, 60]),
    ([2, 1, 1, 1, 1], [16, 24, -24, -16]),
    ([2, 2, 1, 1], [4, -4, -4, 4]),
    ([3, 1, 1, 1], [3, -3, -3, 3]),
    ([3, 2, 1], [1, -3, 3, -1]),
    ])

hess_left[(0, 0, 2, 2, 3, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [60, 300, 300, 60]),
    ([2, 1, 1, 1, 1], [16, 56, 56, 16]),
    ([2, 2, 1, 1], [4, 12, 12, 4]),
    ([3, 1, 1, 1], [3, 6, 6, 3]),
    ([3, 2, 1], [1, 2, 2, 1]),
    ])

hess_right[(0, 0, 2, 2, 4, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [90, 270, 270, 90]),
    ([2, 1, 1, 1, 1], [18, 18, -18, -18]),
    ([2, 2, 1, 1], [6, -6, -6, 6]),
    ([2, 2, 2], [6, -18, 18, -6]),
    ])

hess_left[(0, 0, 2, 2, 4, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [90, 270, 270, 90]),
    ([2, 1, 1, 1, 1], [18, 54, 54, 18]),
    ([2, 2, 1, 1], [6, 18, 18, 6]),
    ([2, 2, 2], [6, 18, 18, 6]),
    ])

hess_right[(0, 0, 2, 2, 4, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [180, 360, 180]),
    ([2, 1, 1, 1, 1], [24, 0, -24]),
    ([2, 2, 1, 1], [4, -8, 4]),
    ])

hess_left[(0, 0, 2, 2, 4, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [180, 360, 180]),
    ([2, 1, 1, 1, 1], [24, 48, 24]),
    ([2, 2, 1, 1], [4, 8, 4]),
    ])

hess_right[(0, 0, 2, 3, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [60, 180, 240, 180, 60]),
    ([2, 1, 1, 1, 1], [16, 16, 0, -16, -16]),
    ([2, 2, 1, 1], [4, -4, 0, -4, 4]),
    ([3, 1, 1, 1], [3, 0, -6, 0, 3]),
    ([3, 2, 1], [1, -2, 0, 2, -1]),
    ])

hess_left[(0, 0, 2, 3, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [60, 180, 240, 180, 60]),
    ([2, 1, 1, 1, 1], [16, 48, 64, 48, 16]),
    ([2, 2, 1, 1], [4, 12, 16, 12, 4]),
    ([3, 1, 1, 1], [3, 9, 12, 9, 3]),
    ([3, 2, 1], [1, 3, 4, 3, 1]),
    ])

hess_right[(0, 0, 2, 3, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [60, 300, 300, 60]),
    ([2, 1, 1, 1, 1], [16, 24, -24, -16]),
    ([2, 2, 1, 1], [4, -4, -4, 4]),
    ([3, 1, 1, 1], [3, -3, -3, 3]),
    ([3, 2, 1], [1, -3, 3, -1]),
    ])

hess_left[(0, 0, 2, 3, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [60, 300, 300, 60]),
    ([2, 1, 1, 1, 1], [16, 56, 56, 16]),
    ([2, 2, 1, 1], [4, 12, 12, 4]),
    ([3, 1, 1, 1], [3, 6, 6, 3]),
    ([3, 2, 1], [1, 2, 2, 1]),
    ])

hess_right[(0, 0, 2, 3, 3, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [180, 360, 180]),
    ([2, 1, 1, 1, 1], [24, 0, -24]),
    ([2, 2, 1, 1], [4, -8, 4]),
    ])

hess_left[(0, 0, 2, 3, 3, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [180, 360, 180]),
    ([2, 1, 1, 1, 1], [24, 48, 24]),
    ([2, 2, 1, 1], [4, 8, 4]),
    ])

hess_right[(0, 0, 2, 3, 4, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [180, 360, 180]),
    ([2, 1, 1, 1, 1], [24, 0, -24]),
    ([2, 2, 1, 1], [4, -8, 4]),
    ])

hess_left[(0, 0, 2, 3, 4, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [180, 360, 180]),
    ([2, 1, 1, 1, 1], [24, 48, 24]),
    ([2, 2, 1, 1], [4, 8, 4]),
    ])

hess_right[(0, 0, 2, 3, 4, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [360, 360]),
    ([2, 1, 1, 1, 1], [24, -24]),
    ])

hess_left[(0, 0, 2, 3, 4, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [360, 360]),
    ([2, 1, 1, 1, 1], [24, 24]),
    ])

hess_right[(0, 1, 1, 1, 1, 1)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 24, 54, 90, 120, 132, 120, 90, 54, 24, 6]),
    ([2, 1, 1, 1, 1], [4, 8, 12, 12, 8, 0, -8, -12, -12, -8, -4]),
    ([2, 2, 1, 1], [2, 0, 2, -2, 0, -4, 0, -2, 2, 0, 2]),
    ([3, 1, 1, 1], [3, 3, 0, 0, -3, -6, -3, 0, 0, 3, 3]),
    ([3, 2, 1], [1, -1, 0, 0, -1, 0, 1, 0, 0, 1, -1]),
    ([4, 1, 1], [2, 0, -2, -2, 0, 0, 0, 2, 2, 0, -2]),
    ([5, 1], [1, -1, -1, 0, 0, 2, 0, 0, -1, -1, 1]),
    ])

hess_left[(0, 1, 1, 1, 1, 1)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 24, 54, 90, 120, 132, 120, 90, 54, 24, 6]),
    ([2, 1, 1, 1, 1], [4, 16, 36, 60, 80, 88, 80, 60, 36, 16, 4]),
    ([2, 2, 1, 1], [2, 8, 18, 30, 40, 44, 40, 30, 18, 8, 2]),
    ([3, 1, 1, 1], [3, 12, 27, 45, 60, 66, 60, 45, 27, 12, 3]),
    ([3, 2, 1], [1, 4, 9, 15, 20, 22, 20, 15, 9, 4, 1]),
    ([4, 1, 1], [2, 8, 18, 30, 40, 44, 40, 30, 18, 8, 2]),
    ([5, 1], [1, 4, 9, 15, 20, 22, 20, 15, 9, 4, 1]),
    ])

hess_right[(0, 1, 1, 1, 1, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 24, 54, 114, 162, 162, 114, 54, 24, 6]),
    ([2, 1, 1, 1, 1], [4, 8, 12, 20, 4, -4, -20, -12, -8, -4]),
    ([2, 2, 1, 1], [2, 0, 2, -2, -2, -2, -2, 2, 0, 2]),
    ([3, 1, 1, 1], [3, 3, 0, 3, -9, -9, 3, 0, 3, 3]),
    ([3, 2, 1], [1, -1, 0, -1, 1, -1, 1, 0, 1, -1]),
    ([4, 1, 1], [2, 0, -2, -2, -2, 2, 2, 2, 0, -2]),
    ([5, 1], [1, -1, -1, -1, 2, 2, -1, -1, -1, 1]),
    ])

hess_left[(0, 1, 1, 1, 1, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 24, 54, 114, 162, 162, 114, 54, 24, 6]),
    ([2, 1, 1, 1, 1], [4, 16, 36, 68, 92, 92, 68, 36, 16, 4]),
    ([2, 2, 1, 1], [2, 8, 18, 30, 38, 38, 30, 18, 8, 2]),
    ([3, 1, 1, 1], [3, 12, 27, 48, 63, 63, 48, 27, 12, 3]),
    ([3, 2, 1], [1, 4, 9, 14, 17, 17, 14, 9, 4, 1]),
    ([4, 1, 1], [2, 8, 18, 30, 38, 38, 30, 18, 8, 2]),
    ([5, 1], [1, 4, 9, 14, 17, 17, 14, 9, 4, 1]),
    ])

hess_right[(0, 1, 1, 1, 1, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 24, 78, 156, 192, 156, 78, 24, 6]),
    ([2, 1, 1, 1, 1], [4, 8, 20, 16, 0, -16, -20, -8, -4]),
    ([2, 2, 1, 1], [2, 0, 2, -4, 0, -4, 2, 0, 2]),
    ([3, 1, 1, 1], [3, 3, 3, -3, -12, -3, 3, 3, 3]),
    ([3, 2, 1], [1, -1, -1, 1, 0, -1, 1, 1, -1]),
    ([4, 1, 1], [2, 0, -2, -4, 0, 4, 2, 0, -2]),
    ([5, 1], [1, -1, -2, 1, 2, 1, -2, -1, 1]),
    ])

hess_left[(0, 1, 1, 1, 1, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 24, 78, 156, 192, 156, 78, 24, 6]),
    ([2, 1, 1, 1, 1], [4, 16, 44, 80, 96, 80, 44, 16, 4]),
    ([2, 2, 1, 1], [2, 8, 18, 28, 32, 28, 18, 8, 2]),
    ([3, 1, 1, 1], [3, 12, 30, 51, 60, 51, 30, 12, 3]),
    ([3, 2, 1], [1, 4, 8, 11, 12, 11, 8, 4, 1]),
    ([4, 1, 1], [2, 8, 18, 28, 32, 28, 18, 8, 2]),
    ([5, 1], [1, 4, 8, 11, 12, 11, 8, 4, 1]),
    ])

hess_right[(0, 1, 1, 1, 1, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 48, 120, 186, 186, 120, 48, 6]),
    ([2, 1, 1, 1, 1], [4, 16, 16, 12, -12, -16, -16, -4]),
    ([2, 2, 1, 1], [2, 0, 0, -2, -2, 0, 0, 2]),
    ([3, 1, 1, 1], [3, 6, -3, -6, -6, -3, 6, 3]),
    ([3, 2, 1], [1, -2, 1, 0, 0, -1, 2, -1]),
    ([4, 1, 1], [2, 0, -4, -2, 2, 4, 0, -2]),
    ([5, 1], [1, -2, 0, 1, 1, 0, -2, 1]),
    ])

hess_left[(0, 1, 1, 1, 1, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 48, 120, 186, 186, 120, 48, 6]),
    ([2, 1, 1, 1, 1], [4, 24, 56, 84, 84, 56, 24, 4]),
    ([2, 2, 1, 1], [2, 8, 16, 22, 22, 16, 8, 2]),
    ([3, 1, 1, 1], [3, 15, 33, 48, 48, 33, 15, 3]),
    ([3, 2, 1], [1, 3, 5, 6, 6, 5, 3, 1]),
    ([4, 1, 1], [2, 8, 16, 22, 22, 16, 8, 2]),
    ([5, 1], [1, 3, 5, 6, 6, 5, 3, 1]),
    ])

hess_right[(0, 1, 1, 1, 1, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [30, 90, 150, 180, 150, 90, 30]),
    ([2, 1, 1, 1, 1], [12, 12, 12, 0, -12, -12, -12]),
    ([2, 2, 1, 1], [2, -2, 2, -4, 2, -2, 2]),
    ([3, 1, 1, 1], [6, 0, -6, 0, -6, 0, 6]),
    ([4, 1, 1], [2, -2, -2, 0, 2, 2, -2]),
    ])

hess_left[(0, 1, 1, 1, 1, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [30, 90, 150, 180, 150, 90, 30]),
    ([2, 1, 1, 1, 1], [12, 36, 60, 72, 60, 36, 12]),
    ([2, 2, 1, 1], [2, 6, 10, 12, 10, 6, 2]),
    ([3, 1, 1, 1], [6, 18, 30, 36, 30, 18, 6]),
    ([4, 1, 1], [2, 6, 10, 12, 10, 6, 2]),
    ])

hess_right[(0, 1, 1, 1, 2, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 24, 78, 156, 192, 156, 78, 24, 6]),
    ([2, 1, 1, 1, 1], [4, 8, 20, 16, 0, -16, -20, -8, -4]),
    ([2, 2, 1, 1], [2, 0, 2, -4, 0, -4, 2, 0, 2]),
    ([3, 1, 1, 1], [3, 3, 3, -3, -12, -3, 3, 3, 3]),
    ([3, 2, 1], [1, -1, -1, 1, 0, -1, 1, 1, -1]),
    ([4, 1, 1], [2, 0, -2, -4, 0, 4, 2, 0, -2]),
    ([5, 1], [1, -1, -2, 1, 2, 1, -2, -1, 1]),
    ])

hess_left[(0, 1, 1, 1, 2, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 24, 78, 156, 192, 156, 78, 24, 6]),
    ([2, 1, 1, 1, 1], [4, 16, 44, 80, 96, 80, 44, 16, 4]),
    ([2, 2, 1, 1], [2, 8, 18, 28, 32, 28, 18, 8, 2]),
    ([3, 1, 1, 1], [3, 12, 30, 51, 60, 51, 30, 12, 3]),
    ([3, 2, 1], [1, 4, 8, 11, 12, 11, 8, 4, 1]),
    ([4, 1, 1], [2, 8, 18, 28, 32, 28, 18, 8, 2]),
    ([5, 1], [1, 4, 8, 11, 12, 11, 8, 4, 1]),
    ])

hess_right[(0, 1, 1, 1, 2, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 24, 102, 228, 228, 102, 24, 6]),
    ([2, 1, 1, 1, 1], [4, 8, 28, 16, -16, -28, -8, -4]),
    ([2, 2, 1, 1], [2, 0, 2, -4, -4, 2, 0, 2]),
    ([3, 1, 1, 1], [3, 3, 6, -12, -12, 6, 3, 3]),
    ([3, 2, 1], [1, -1, -2, 4, -4, 2, 1, -1]),
    ([4, 1, 1], [2, 0, -2, -8, 8, 2, 0, -2]),
    ([5, 1], [1, -1, -3, 3, 3, -3, -1, 1]),
    ])

hess_left[(0, 1, 1, 1, 2, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 24, 102, 228, 228, 102, 24, 6]),
    ([2, 1, 1, 1, 1], [4, 16, 52, 96, 96, 52, 16, 4]),
    ([2, 2, 1, 1], [2, 8, 18, 28, 28, 18, 8, 2]),
    ([3, 1, 1, 1], [3, 12, 33, 51, 51, 33, 12, 3]),
    ([3, 2, 1], [1, 4, 7, 9, 9, 7, 4, 1]),
    ([4, 1, 1], [2, 8, 18, 24, 24, 18, 8, 2]),
    ([5, 1], [1, 4, 7, 8, 8, 7, 4, 1]),
    ])

hess_right[(0, 1, 1, 1, 2, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 48, 174, 264, 174, 48, 6]),
    ([2, 1, 1, 1, 1], [4, 16, 28, 0, -28, -16, -4]),
    ([2, 2, 1, 1], [2, 0, 2, -8, 2, 0, 2]),
    ([3, 1, 1, 1], [3, 6, -3, -12, -3, 6, 3]),
    ([3, 2, 1], [1, -2, 1, 0, -1, 2, -1]),
    ([4, 1, 1], [2, 0, -6, 0, 6, 0, -2]),
    ([5, 1], [1, -2, -1, 4, -1, -2, 1]),
    ])

hess_left[(0, 1, 1, 1, 2, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 48, 174, 264, 174, 48, 6]),
    ([2, 1, 1, 1, 1], [4, 24, 68, 96, 68, 24, 4]),
    ([2, 2, 1, 1], [2, 8, 18, 24, 18, 8, 2]),
    ([3, 1, 1, 1], [3, 15, 33, 42, 33, 15, 3]),
    ([3, 2, 1], [1, 3, 5, 6, 5, 3, 1]),
    ([4, 1, 1], [2, 8, 14, 16, 14, 8, 2]),
    ([5, 1], [1, 3, 4, 4, 4, 3, 1]),
    ])

hess_right[(0, 1, 1, 1, 2, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [30, 90, 240, 240, 90, 30]),
    ([2, 1, 1, 1, 1], [12, 12, 24, -24, -12, -12]),
    ([2, 2, 1, 1], [2, -2, 0, 0, -2, 2]),
    ([3, 1, 1, 1], [6, 0, -6, -6, 0, 6]),
    ([4, 1, 1], [2, -2, -4, 4, 2, -2]),
    ])

hess_left[(0, 1, 1, 1, 2, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [30, 90, 240, 240, 90, 30]),
    ([2, 1, 1, 1, 1], [12, 36, 72, 72, 36, 12]),
    ([2, 2, 1, 1], [2, 6, 8, 8, 6, 2]),
    ([3, 1, 1, 1], [6, 18, 30, 30, 18, 6]),
    ([4, 1, 1], [2, 6, 8, 8, 6, 2]),
    ])

hess_right[(0, 1, 1, 1, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 48, 174, 264, 174, 48, 6]),
    ([2, 1, 1, 1, 1], [4, 16, 28, 0, -28, -16, -4]),
    ([2, 2, 1, 1], [2, 0, 2, -8, 2, 0, 2]),
    ([3, 1, 1, 1], [3, 6, -3, -12, -3, 6, 3]),
    ([3, 2, 1], [1, -2, 1, 0, -1, 2, -1]),
    ([4, 1, 1], [2, 0, -6, 0, 6, 0, -2]),
    ([5, 1], [1, -2, -1, 4, -1, -2, 1]),
    ])

hess_left[(0, 1, 1, 1, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 48, 174, 264, 174, 48, 6]),
    ([2, 1, 1, 1, 1], [4, 24, 68, 96, 68, 24, 4]),
    ([2, 2, 1, 1], [2, 8, 18, 24, 18, 8, 2]),
    ([3, 1, 1, 1], [3, 15, 33, 42, 33, 15, 3]),
    ([3, 2, 1], [1, 3, 5, 6, 5, 3, 1]),
    ([4, 1, 1], [2, 8, 14, 16, 14, 8, 2]),
    ([5, 1], [1, 3, 4, 4, 4, 3, 1]),
    ])

hess_right[(0, 1, 1, 1, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 102, 252, 252, 102, 6]),
    ([2, 1, 1, 1, 1], [4, 28, 16, -16, -28, -4]),
    ([2, 2, 1, 1], [2, 2, -4, -4, 2, 2]),
    ([3, 1, 1, 1], [3, 6, -9, -9, 6, 3]),
    ([3, 2, 1], [1, -2, 1, -1, 2, -1]),
    ([4, 1, 1], [2, -2, -4, 4, 2, -2]),
    ([5, 1], [1, -3, 2, 2, -3, 1]),
    ])

hess_left[(0, 1, 1, 1, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 102, 252, 252, 102, 6]),
    ([2, 1, 1, 1, 1], [4, 36, 80, 80, 36, 4]),
    ([2, 2, 1, 1], [2, 10, 20, 20, 10, 2]),
    ([3, 1, 1, 1], [3, 15, 27, 27, 15, 3]),
    ([3, 2, 1], [1, 3, 5, 5, 3, 1]),
    ([4, 1, 1], [2, 6, 8, 8, 6, 2]),
    ([5, 1], [1, 2, 2, 2, 2, 1]),
    ])

hess_right[(0, 1, 1, 1, 3, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [30, 180, 300, 180, 30]),
    ([2, 1, 1, 1, 1], [12, 24, 0, -24, -12]),
    ([2, 2, 1, 1], [2, -4, 4, -4, 2]),
    ([3, 1, 1, 1], [6, 0, -12, 0, 6]),
    ([4, 1, 1], [2, -4, 0, 4, -2]),
    ])

hess_left[(0, 1, 1, 1, 3, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [30, 180, 300, 180, 30]),
    ([2, 1, 1, 1, 1], [12, 48, 72, 48, 12]),
    ([2, 2, 1, 1], [2, 4, 4, 4, 2]),
    ([3, 1, 1, 1], [6, 18, 24, 18, 6]),
    ([4, 1, 1], [2, 4, 4, 4, 2]),
    ])

hess_right[(0, 1, 1, 1, 4, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [60, 180, 240, 180, 60]),
    ([2, 1, 1, 1, 1], [16, 16, 0, -16, -16]),
    ([2, 2, 1, 1], [4, -4, 0, -4, 4]),
    ([3, 1, 1, 1], [3, 0, -6, 0, 3]),
    ([3, 2, 1], [1, -2, 0, 2, -1]),
    ])

hess_left[(0, 1, 1, 1, 4, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [60, 180, 240, 180, 60]),
    ([2, 1, 1, 1, 1], [16, 48, 64, 48, 16]),
    ([2, 2, 1, 1], [4, 12, 16, 12, 4]),
    ([3, 1, 1, 1], [3, 9, 12, 9, 3]),
    ([3, 2, 1], [1, 3, 4, 3, 1]),
    ])

hess_right[(0, 1, 1, 1, 4, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [120, 240, 240, 120]),
    ([2, 1, 1, 1, 1], [24, 0, 0, -24]),
    ([3, 1, 1, 1], [6, -6, -6, 6]),
    ])

hess_left[(0, 1, 1, 1, 4, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [120, 240, 240, 120]),
    ([2, 1, 1, 1, 1], [24, 48, 48, 24]),
    ([3, 1, 1, 1], [6, 12, 12, 6]),
    ])

hess_right[(0, 1, 1, 2, 2, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 48, 120, 186, 186, 120, 48, 6]),
    ([2, 1, 1, 1, 1], [4, 16, 16, 12, -12, -16, -16, -4]),
    ([2, 2, 1, 1], [2, 0, 0, -2, -2, 0, 0, 2]),
    ([3, 1, 1, 1], [3, 6, -3, -6, -6, -3, 6, 3]),
    ([3, 2, 1], [1, -2, 1, 0, 0, -1, 2, -1]),
    ([4, 1, 1], [2, 0, -4, -2, 2, 4, 0, -2]),
    ([5, 1], [1, -2, 0, 1, 1, 0, -2, 1]),
    ])

hess_left[(0, 1, 1, 2, 2, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 48, 120, 186, 186, 120, 48, 6]),
    ([2, 1, 1, 1, 1], [4, 24, 56, 84, 84, 56, 24, 4]),
    ([2, 2, 1, 1], [2, 8, 16, 22, 22, 16, 8, 2]),
    ([3, 1, 1, 1], [3, 15, 33, 48, 48, 33, 15, 3]),
    ([3, 2, 1], [1, 3, 5, 6, 6, 5, 3, 1]),
    ([4, 1, 1], [2, 8, 16, 22, 22, 16, 8, 2]),
    ([5, 1], [1, 3, 5, 6, 6, 5, 3, 1]),
    ])

hess_right[(0, 1, 1, 2, 2, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 48, 174, 264, 174, 48, 6]),
    ([2, 1, 1, 1, 1], [4, 16, 28, 0, -28, -16, -4]),
    ([2, 2, 1, 1], [2, 0, 2, -8, 2, 0, 2]),
    ([3, 1, 1, 1], [3, 6, -3, -12, -3, 6, 3]),
    ([3, 2, 1], [1, -2, 1, 0, -1, 2, -1]),
    ([4, 1, 1], [2, 0, -6, 0, 6, 0, -2]),
    ([5, 1], [1, -2, -1, 4, -1, -2, 1]),
    ])

hess_left[(0, 1, 1, 2, 2, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 48, 174, 264, 174, 48, 6]),
    ([2, 1, 1, 1, 1], [4, 24, 68, 96, 68, 24, 4]),
    ([2, 2, 1, 1], [2, 8, 18, 24, 18, 8, 2]),
    ([3, 1, 1, 1], [3, 15, 33, 42, 33, 15, 3]),
    ([3, 2, 1], [1, 3, 5, 6, 5, 3, 1]),
    ([4, 1, 1], [2, 8, 14, 16, 14, 8, 2]),
    ([5, 1], [1, 3, 4, 4, 4, 3, 1]),
    ])

hess_right[(0, 1, 1, 2, 2, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 72, 282, 282, 72, 6]),
    ([2, 1, 1, 1, 1], [4, 24, 28, -28, -24, -4]),
    ([2, 2, 1, 1], [2, 0, -2, -2, 0, 2]),
    ([3, 1, 1, 1], [3, 9, -12, -12, 9, 3]),
    ([3, 2, 1], [1, -3, 4, -4, 3, -1]),
    ([4, 1, 1], [2, 0, -10, 10, 0, -2]),
    ([5, 1], [1, -3, 2, 2, -3, 1]),
    ])

hess_left[(0, 1, 1, 2, 2, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 72, 282, 282, 72, 6]),
    ([2, 1, 1, 1, 1], [4, 32, 84, 84, 32, 4]),
    ([2, 2, 1, 1], [2, 8, 14, 14, 8, 2]),
    ([3, 1, 1, 1], [3, 18, 33, 33, 18, 3]),
    ([3, 2, 1], [1, 2, 3, 3, 2, 1]),
    ([4, 1, 1], [2, 8, 10, 10, 8, 2]),
    ([5, 1], [1, 2, 2, 2, 2, 1]),
    ])

hess_right[(0, 1, 1, 2, 2, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [30, 180, 300, 180, 30]),
    ([2, 1, 1, 1, 1], [12, 24, 0, -24, -12]),
    ([2, 2, 1, 1], [2, -4, 4, -4, 2]),
    ([3, 1, 1, 1], [6, 0, -12, 0, 6]),
    ([4, 1, 1], [2, -4, 0, 4, -2]),
    ])

hess_left[(0, 1, 1, 2, 2, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [30, 180, 300, 180, 30]),
    ([2, 1, 1, 1, 1], [12, 48, 72, 48, 12]),
    ([2, 2, 1, 1], [2, 4, 4, 4, 2]),
    ([3, 1, 1, 1], [6, 18, 24, 18, 6]),
    ([4, 1, 1], [2, 4, 4, 4, 2]),
    ])

hess_right[(0, 1, 1, 2, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 102, 252, 252, 102, 6]),
    ([2, 1, 1, 1, 1], [4, 28, 16, -16, -28, -4]),
    ([2, 2, 1, 1], [2, 2, -4, -4, 2, 2]),
    ([3, 1, 1, 1], [3, 6, -9, -9, 6, 3]),
    ([3, 2, 1], [1, -2, 1, -1, 2, -1]),
    ([4, 1, 1], [2, -2, -4, 4, 2, -2]),
    ([5, 1], [1, -3, 2, 2, -3, 1]),
    ])

hess_left[(0, 1, 1, 2, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 102, 252, 252, 102, 6]),
    ([2, 1, 1, 1, 1], [4, 36, 80, 80, 36, 4]),
    ([2, 2, 1, 1], [2, 10, 20, 20, 10, 2]),
    ([3, 1, 1, 1], [3, 15, 27, 27, 15, 3]),
    ([3, 2, 1], [1, 3, 5, 5, 3, 1]),
    ([4, 1, 1], [2, 6, 8, 8, 6, 2]),
    ([5, 1], [1, 2, 2, 2, 2, 1]),
    ])

hess_right[(0, 1, 1, 2, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 156, 396, 156, 6]),
    ([2, 1, 1, 1, 1], [4, 40, 0, -40, -4]),
    ([2, 2, 1, 1], [2, 4, -12, 4, 2]),
    ([3, 1, 1, 1], [3, 6, -18, 6, 3]),
    ([3, 2, 1], [1, -2, 0, 2, -1]),
    ([4, 1, 1], [2, -4, 0, 4, -2]),
    ([5, 1], [1, -4, 6, -4, 1]),
    ])

hess_left[(0, 1, 1, 2, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [6, 156, 396, 156, 6]),
    ([2, 1, 1, 1, 1], [4, 48, 88, 48, 4]),
    ([2, 2, 1, 1], [2, 12, 20, 12, 2]),
    ([3, 1, 1, 1], [3, 15, 18, 15, 3]),
    ([3, 2, 1], [1, 3, 4, 3, 1]),
    ([4, 1, 1], [2, 4, 4, 4, 2]),
    ([5, 1], [1, 1, 1, 1, 1]),
    ])

hess_right[(0, 1, 1, 2, 3, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [30, 330, 330, 30]),
    ([2, 1, 1, 1, 1], [12, 36, -36, -12]),
    ([2, 2, 1, 1], [2, -2, -2, 2]),
    ([3, 1, 1, 1], [6, -6, -6, 6]),
    ([4, 1, 1], [2, -6, 6, -2]),
    ])

hess_left[(0, 1, 1, 2, 3, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [30, 330, 330, 30]),
    ([2, 1, 1, 1, 1], [12, 60, 60, 12]),
    ([2, 2, 1, 1], [2, 6, 6, 2]),
    ([3, 1, 1, 1], [6, 12, 12, 6]),
    ([4, 1, 1], [2, 2, 2, 2]),
    ])

hess_right[(0, 1, 1, 2, 4, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [60, 300, 300, 60]),
    ([2, 1, 1, 1, 1], [16, 24, -24, -16]),
    ([2, 2, 1, 1], [4, -4, -4, 4]),
    ([3, 1, 1, 1], [3, -3, -3, 3]),
    ([3, 2, 1], [1, -3, 3, -1]),
    ])

hess_left[(0, 1, 1, 2, 4, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [60, 300, 300, 60]),
    ([2, 1, 1, 1, 1], [16, 56, 56, 16]),
    ([2, 2, 1, 1], [4, 12, 12, 4]),
    ([3, 1, 1, 1], [3, 6, 6, 3]),
    ([3, 2, 1], [1, 2, 2, 1]),
    ])

hess_right[(0, 1, 1, 2, 4, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [120, 480, 120]),
    ([2, 1, 1, 1, 1], [24, 0, -24]),
    ([3, 1, 1, 1], [6, -12, 6]),
    ])

hess_left[(0, 1, 1, 2, 4, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [120, 480, 120]),
    ([2, 1, 1, 1, 1], [24, 48, 24]),
    ([3, 1, 1, 1], [6, 6, 6]),
    ])

hess_right[(0, 1, 1, 3, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [60, 180, 240, 180, 60]),
    ([2, 1, 1, 1, 1], [16, 16, 0, -16, -16]),
    ([2, 2, 1, 1], [4, -4, 0, -4, 4]),
    ([3, 1, 1, 1], [3, 0, -6, 0, 3]),
    ([3, 2, 1], [1, -2, 0, 2, -1]),
    ])

hess_left[(0, 1, 1, 3, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [60, 180, 240, 180, 60]),
    ([2, 1, 1, 1, 1], [16, 48, 64, 48, 16]),
    ([2, 2, 1, 1], [4, 12, 16, 12, 4]),
    ([3, 1, 1, 1], [3, 9, 12, 9, 3]),
    ([3, 2, 1], [1, 3, 4, 3, 1]),
    ])

hess_right[(0, 1, 1, 3, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [60, 300, 300, 60]),
    ([2, 1, 1, 1, 1], [16, 24, -24, -16]),
    ([2, 2, 1, 1], [4, -4, -4, 4]),
    ([3, 1, 1, 1], [3, -3, -3, 3]),
    ([3, 2, 1], [1, -3, 3, -1]),
    ])

hess_left[(0, 1, 1, 3, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [60, 300, 300, 60]),
    ([2, 1, 1, 1, 1], [16, 56, 56, 16]),
    ([2, 2, 1, 1], [4, 12, 12, 4]),
    ([3, 1, 1, 1], [3, 6, 6, 3]),
    ([3, 2, 1], [1, 2, 2, 1]),
    ])

hess_right[(0, 1, 1, 3, 3, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [180, 360, 180]),
    ([2, 1, 1, 1, 1], [24, 0, -24]),
    ([2, 2, 1, 1], [4, -8, 4]),
    ])

hess_left[(0, 1, 1, 3, 3, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [180, 360, 180]),
    ([2, 1, 1, 1, 1], [24, 48, 24]),
    ([2, 2, 1, 1], [4, 8, 4]),
    ])

hess_right[(0, 1, 1, 3, 4, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [180, 360, 180]),
    ([2, 1, 1, 1, 1], [24, 0, -24]),
    ([2, 2, 1, 1], [4, -8, 4]),
    ])

hess_left[(0, 1, 1, 3, 4, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [180, 360, 180]),
    ([2, 1, 1, 1, 1], [24, 48, 24]),
    ([2, 2, 1, 1], [4, 8, 4]),
    ])

hess_right[(0, 1, 1, 3, 4, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [360, 360]),
    ([2, 1, 1, 1, 1], [24, -24]),
    ])

hess_left[(0, 1, 1, 3, 4, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [360, 360]),
    ([2, 1, 1, 1, 1], [24, 24]),
    ])

hess_right[(0, 1, 2, 2, 2, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [30, 90, 150, 180, 150, 90, 30]),
    ([2, 1, 1, 1, 1], [12, 12, 12, 0, -12, -12, -12]),
    ([2, 2, 1, 1], [2, -2, 2, -4, 2, -2, 2]),
    ([3, 1, 1, 1], [6, 0, -6, 0, -6, 0, 6]),
    ([4, 1, 1], [2, -2, -2, 0, 2, 2, -2]),
    ])

hess_left[(0, 1, 2, 2, 2, 2)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [30, 90, 150, 180, 150, 90, 30]),
    ([2, 1, 1, 1, 1], [12, 36, 60, 72, 60, 36, 12]),
    ([2, 2, 1, 1], [2, 6, 10, 12, 10, 6, 2]),
    ([3, 1, 1, 1], [6, 18, 30, 36, 30, 18, 6]),
    ([4, 1, 1], [2, 6, 10, 12, 10, 6, 2]),
    ])

hess_right[(0, 1, 2, 2, 2, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [30, 90, 240, 240, 90, 30]),
    ([2, 1, 1, 1, 1], [12, 12, 24, -24, -12, -12]),
    ([2, 2, 1, 1], [2, -2, 0, 0, -2, 2]),
    ([3, 1, 1, 1], [6, 0, -6, -6, 0, 6]),
    ([4, 1, 1], [2, -2, -4, 4, 2, -2]),
    ])

hess_left[(0, 1, 2, 2, 2, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [30, 90, 240, 240, 90, 30]),
    ([2, 1, 1, 1, 1], [12, 36, 72, 72, 36, 12]),
    ([2, 2, 1, 1], [2, 6, 8, 8, 6, 2]),
    ([3, 1, 1, 1], [6, 18, 30, 30, 18, 6]),
    ([4, 1, 1], [2, 6, 8, 8, 6, 2]),
    ])

hess_right[(0, 1, 2, 2, 2, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [30, 180, 300, 180, 30]),
    ([2, 1, 1, 1, 1], [12, 24, 0, -24, -12]),
    ([2, 2, 1, 1], [2, -4, 4, -4, 2]),
    ([3, 1, 1, 1], [6, 0, -12, 0, 6]),
    ([4, 1, 1], [2, -4, 0, 4, -2]),
    ])

hess_left[(0, 1, 2, 2, 2, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [30, 180, 300, 180, 30]),
    ([2, 1, 1, 1, 1], [12, 48, 72, 48, 12]),
    ([2, 2, 1, 1], [2, 4, 4, 4, 2]),
    ([3, 1, 1, 1], [6, 18, 24, 18, 6]),
    ([4, 1, 1], [2, 4, 4, 4, 2]),
    ])

hess_right[(0, 1, 2, 2, 2, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [120, 240, 240, 120]),
    ([2, 1, 1, 1, 1], [24, 0, 0, -24]),
    ([3, 1, 1, 1], [6, -6, -6, 6]),
    ])

hess_left[(0, 1, 2, 2, 2, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [120, 240, 240, 120]),
    ([2, 1, 1, 1, 1], [24, 48, 48, 24]),
    ([3, 1, 1, 1], [6, 12, 12, 6]),
    ])

hess_right[(0, 1, 2, 2, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [30, 180, 300, 180, 30]),
    ([2, 1, 1, 1, 1], [12, 24, 0, -24, -12]),
    ([2, 2, 1, 1], [2, -4, 4, -4, 2]),
    ([3, 1, 1, 1], [6, 0, -12, 0, 6]),
    ([4, 1, 1], [2, -4, 0, 4, -2]),
    ])

hess_left[(0, 1, 2, 2, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [30, 180, 300, 180, 30]),
    ([2, 1, 1, 1, 1], [12, 48, 72, 48, 12]),
    ([2, 2, 1, 1], [2, 4, 4, 4, 2]),
    ([3, 1, 1, 1], [6, 18, 24, 18, 6]),
    ([4, 1, 1], [2, 4, 4, 4, 2]),
    ])

hess_right[(0, 1, 2, 2, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [30, 330, 330, 30]),
    ([2, 1, 1, 1, 1], [12, 36, -36, -12]),
    ([2, 2, 1, 1], [2, -2, -2, 2]),
    ([3, 1, 1, 1], [6, -6, -6, 6]),
    ([4, 1, 1], [2, -6, 6, -2]),
    ])

hess_left[(0, 1, 2, 2, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [30, 330, 330, 30]),
    ([2, 1, 1, 1, 1], [12, 60, 60, 12]),
    ([2, 2, 1, 1], [2, 6, 6, 2]),
    ([3, 1, 1, 1], [6, 12, 12, 6]),
    ([4, 1, 1], [2, 2, 2, 2]),
    ])

hess_right[(0, 1, 2, 2, 3, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [120, 480, 120]),
    ([2, 1, 1, 1, 1], [24, 0, -24]),
    ([3, 1, 1, 1], [6, -12, 6]),
    ])

hess_left[(0, 1, 2, 2, 3, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [120, 480, 120]),
    ([2, 1, 1, 1, 1], [24, 48, 24]),
    ([3, 1, 1, 1], [6, 6, 6]),
    ])

hess_right[(0, 1, 2, 2, 4, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [180, 360, 180]),
    ([2, 1, 1, 1, 1], [24, 0, -24]),
    ([2, 2, 1, 1], [4, -8, 4]),
    ])

hess_left[(0, 1, 2, 2, 4, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [180, 360, 180]),
    ([2, 1, 1, 1, 1], [24, 48, 24]),
    ([2, 2, 1, 1], [4, 8, 4]),
    ])

hess_right[(0, 1, 2, 2, 4, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [360, 360]),
    ([2, 1, 1, 1, 1], [24, -24]),
    ])

hess_left[(0, 1, 2, 2, 4, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [360, 360]),
    ([2, 1, 1, 1, 1], [24, 24]),
    ])

hess_right[(0, 1, 2, 3, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [120, 240, 240, 120]),
    ([2, 1, 1, 1, 1], [24, 0, 0, -24]),
    ([3, 1, 1, 1], [6, -6, -6, 6]),
    ])

hess_left[(0, 1, 2, 3, 3, 3)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [120, 240, 240, 120]),
    ([2, 1, 1, 1, 1], [24, 48, 48, 24]),
    ([3, 1, 1, 1], [6, 12, 12, 6]),
    ])

hess_right[(0, 1, 2, 3, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [120, 480, 120]),
    ([2, 1, 1, 1, 1], [24, 0, -24]),
    ([3, 1, 1, 1], [6, -12, 6]),
    ])

hess_left[(0, 1, 2, 3, 3, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [120, 480, 120]),
    ([2, 1, 1, 1, 1], [24, 48, 24]),
    ([3, 1, 1, 1], [6, 6, 6]),
    ])

hess_right[(0, 1, 2, 3, 3, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [360, 360]),
    ([2, 1, 1, 1, 1], [24, -24]),
    ])

hess_left[(0, 1, 2, 3, 3, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [360, 360]),
    ([2, 1, 1, 1, 1], [24, 24]),
    ])

hess_right[(0, 1, 2, 3, 4, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [360, 360]),
    ([2, 1, 1, 1, 1], [24, -24]),
    ])

hess_left[(0, 1, 2, 3, 4, 4)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [360, 360]),
    ([2, 1, 1, 1, 1], [24, 24]),
    ])

hess_right[(0, 1, 2, 3, 4, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [720]),
    ])

hess_left[(0, 1, 2, 3, 4, 5)] = p.sum(
    p.term(Partition(index), R(coeffs) / zee(index))
    for index, coeffs in [
    ([1, 1, 1, 1, 1, 1], [720]),
    ])

