
R, q = QQ['q'].objgen()
sym = SymmetricFunctions(R)
m = sym.m()
e = sym.e()

# q-chromatic symmetric functions for unit interval orders.
# There is one of each for each Dyck path.
# Dyck paths are represented as numbers of missing boxes in each row.
# For example, the following Dyck paths:

#   /\
#  /  \
# /    \
# >>> csf[0, 0, 0]

#  /\/\
# /    \
# >>> csf[0, 0, 1]

#  /\
# /  \/\
# >>> csf[0, 0, 2]

#    /\
# /\/  \
# >>> csf[0, 1, 1]

# /\/\/\
# >>> csf[0, 1, 2]

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
    return all(
        is_e_positive(csf[path])
        for path in csf
        )

#--------------------------------
# data
#--------------------------------

csf = {}

csf[(0, 0, 0, 0)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1], [1, 3, 5, 6, 5, 3, 1]),
    ])

csf[(0, 0, 0, 1)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1], [1, 3, 8, 8, 3, 1]),
    ([2, 1, 1], [0, 0, 1, 1]),
    ])

csf[(0, 0, 0, 2)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1], [1, 6, 10, 6, 1]),
    ([2, 1, 1], [0, 1, 2, 1]),
    ])

csf[(0, 0, 0, 3)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1], [4, 8, 8, 4]),
    ([2, 1, 1], [1, 2, 2, 1]),
    ])

csf[(0, 0, 1, 1)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1], [1, 6, 10, 6, 1]),
    ([2, 1, 1], [0, 1, 2, 1]),
    ])

csf[(0, 0, 1, 2)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1], [1, 11, 11, 1]),
    ([2, 1, 1], [0, 3, 3]),
    ([2, 2], [0, 1, 1]),
    ])

csf[(0, 0, 1, 3)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1], [4, 16, 4]),
    ([2, 1, 1], [1, 6, 1]),
    ([2, 2], [0, 2]),
    ([3, 1], [0, 1]),
    ])

csf[(0, 0, 2, 2)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1], [6, 12, 6]),
    ([2, 1, 1], [2, 4, 2]),
    ([2, 2], [1, 2, 1]),
    ])

csf[(0, 0, 2, 3)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1], [12, 12]),
    ([2, 1, 1], [5, 5]),
    ([2, 2], [2, 2]),
    ([3, 1], [1, 1]),
    ])

csf[(0, 1, 1, 1)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1], [4, 8, 8, 4]),
    ([2, 1, 1], [1, 2, 2, 1]),
    ])

csf[(0, 1, 1, 2)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1], [4, 16, 4]),
    ([2, 1, 1], [1, 6, 1]),
    ([2, 2], [0, 2]),
    ([3, 1], [0, 1]),
    ])

csf[(0, 1, 1, 3)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1], [12, 12]),
    ([2, 1, 1], [5, 5]),
    ([2, 2], [2, 2]),
    ([3, 1], [1, 1]),
    ])

csf[(0, 1, 2, 2)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1], [12, 12]),
    ([2, 1, 1], [5, 5]),
    ([2, 2], [2, 2]),
    ([3, 1], [1, 1]),
    ])

csf[(0, 1, 2, 3)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1], [24]),
    ([2, 1, 1], [12]),
    ([2, 2], [6]),
    ([3, 1], [4]),
    ([4], [1]),
    ])

