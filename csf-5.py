
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

csf[(0, 0, 0, 0, 0)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 4, 9, 15, 20, 22, 20, 15, 9, 4, 1]),
    ])

csf[(0, 0, 0, 0, 1)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 4, 9, 19, 27, 27, 19, 9, 4, 1]),
    ([2, 1, 1, 1], [0, 0, 0, 1, 2, 2, 1]),
    ])

csf[(0, 0, 0, 0, 2)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 4, 13, 26, 32, 26, 13, 4, 1]),
    ([2, 1, 1, 1], [0, 0, 1, 3, 4, 3, 1]),
    ])

csf[(0, 0, 0, 0, 3)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 8, 20, 31, 31, 20, 8, 1]),
    ([2, 1, 1, 1], [0, 1, 3, 5, 5, 3, 1]),
    ])

csf[(0, 0, 0, 0, 4)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [5, 15, 25, 30, 25, 15, 5]),
    ([2, 1, 1, 1], [1, 3, 5, 6, 5, 3, 1]),
    ])

csf[(0, 0, 0, 1, 1)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 4, 13, 26, 32, 26, 13, 4, 1]),
    ([2, 1, 1, 1], [0, 0, 1, 3, 4, 3, 1]),
    ])

csf[(0, 0, 0, 1, 2)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 4, 17, 38, 38, 17, 4, 1]),
    ([2, 1, 1, 1], [0, 0, 2, 7, 7, 2]),
    ([2, 2, 1], [0, 0, 0, 1, 1]),
    ])

csf[(0, 0, 0, 1, 3)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 8, 29, 44, 29, 8, 1]),
    ([2, 1, 1, 1], [0, 1, 6, 10, 6, 1]),
    ([2, 2, 1], [0, 0, 1, 2, 1]),
    ])

csf[(0, 0, 0, 1, 4)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [5, 15, 40, 40, 15, 5]),
    ([2, 1, 1, 1], [1, 3, 11, 11, 3, 1]),
    ([2, 2, 1], [0, 0, 2, 2]),
    ([3, 1, 1], [0, 0, 1, 1]),
    ])

csf[(0, 0, 0, 2, 2)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 8, 29, 44, 29, 8, 1]),
    ([2, 1, 1, 1], [0, 1, 6, 10, 6, 1]),
    ([2, 2, 1], [0, 0, 1, 2, 1]),
    ])

csf[(0, 0, 0, 2, 3)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 17, 42, 42, 17, 1]),
    ([2, 1, 1, 1], [0, 4, 11, 11, 4]),
    ([2, 2, 1], [0, 1, 3, 3, 1]),
    ])

csf[(0, 0, 0, 2, 4)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [5, 30, 50, 30, 5]),
    ([2, 1, 1, 1], [1, 9, 16, 9, 1]),
    ([2, 2, 1], [0, 2, 4, 2]),
    ([3, 1, 1], [0, 1, 2, 1]),
    ])

csf[(0, 0, 0, 3, 3)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [10, 30, 40, 30, 10]),
    ([2, 1, 1, 1], [3, 9, 12, 9, 3]),
    ([2, 2, 1], [1, 3, 4, 3, 1]),
    ])

csf[(0, 0, 0, 3, 4)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [20, 40, 40, 20]),
    ([2, 1, 1, 1], [7, 14, 14, 7]),
    ([2, 2, 1], [2, 4, 4, 2]),
    ([3, 1, 1], [1, 2, 2, 1]),
    ])

csf[(0, 0, 1, 1, 1)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 8, 20, 31, 31, 20, 8, 1]),
    ([2, 1, 1, 1], [0, 1, 3, 5, 5, 3, 1]),
    ])

csf[(0, 0, 1, 1, 2)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 8, 29, 44, 29, 8, 1]),
    ([2, 1, 1, 1], [0, 1, 6, 10, 6, 1]),
    ([2, 2, 1], [0, 0, 1, 2, 1]),
    ])

csf[(0, 0, 1, 1, 3)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 12, 47, 47, 12, 1]),
    ([2, 1, 1, 1], [0, 2, 13, 13, 2]),
    ([2, 2, 1], [0, 0, 3, 3]),
    ([3, 1, 1], [0, 0, 1, 1]),
    ])

csf[(0, 0, 1, 1, 4)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [5, 30, 50, 30, 5]),
    ([2, 1, 1, 1], [1, 9, 16, 9, 1]),
    ([2, 2, 1], [0, 2, 4, 2]),
    ([3, 1, 1], [0, 1, 2, 1]),
    ])

csf[(0, 0, 1, 2, 2)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 17, 42, 42, 17, 1]),
    ([2, 1, 1, 1], [0, 4, 11, 11, 4]),
    ([2, 2, 1], [0, 1, 3, 3, 1]),
    ])

csf[(0, 0, 1, 2, 3)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [1, 26, 66, 26, 1]),
    ([2, 1, 1, 1], [0, 7, 22, 7]),
    ([2, 2, 1], [0, 2, 8, 2]),
    ([3, 1, 1], [0, 0, 2]),
    ([3, 2], [0, 0, 1]),
    ])

csf[(0, 0, 1, 2, 4)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [5, 55, 55, 5]),
    ([2, 1, 1, 1], [1, 20, 20, 1]),
    ([2, 2, 1], [0, 7, 7]),
    ([3, 1, 1], [0, 3, 3]),
    ([3, 2], [0, 1, 1]),
    ])

csf[(0, 0, 1, 3, 3)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [10, 50, 50, 10]),
    ([2, 1, 1, 1], [3, 18, 18, 3]),
    ([2, 2, 1], [1, 7, 7, 1]),
    ([3, 1, 1], [0, 2, 2]),
    ([3, 2], [0, 1, 1]),
    ])

csf[(0, 0, 1, 3, 4)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [20, 80, 20]),
    ([2, 1, 1, 1], [7, 34, 7]),
    ([2, 2, 1], [2, 14, 2]),
    ([3, 1, 1], [1, 8, 1]),
    ([3, 2], [0, 3]),
    ([4, 1], [0, 1]),
    ])

csf[(0, 0, 2, 2, 2)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [10, 30, 40, 30, 10]),
    ([2, 1, 1, 1], [3, 9, 12, 9, 3]),
    ([2, 2, 1], [1, 3, 4, 3, 1]),
    ])

csf[(0, 0, 2, 2, 3)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [10, 50, 50, 10]),
    ([2, 1, 1, 1], [3, 18, 18, 3]),
    ([2, 2, 1], [1, 7, 7, 1]),
    ([3, 1, 1], [0, 2, 2]),
    ([3, 2], [0, 1, 1]),
    ])

csf[(0, 0, 2, 2, 4)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [30, 60, 30]),
    ([2, 1, 1, 1], [12, 24, 12]),
    ([2, 2, 1], [5, 10, 5]),
    ([3, 1, 1], [2, 4, 2]),
    ([3, 2], [1, 2, 1]),
    ])

csf[(0, 0, 2, 3, 3)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [30, 60, 30]),
    ([2, 1, 1, 1], [12, 24, 12]),
    ([2, 2, 1], [5, 10, 5]),
    ([3, 1, 1], [2, 4, 2]),
    ([3, 2], [1, 2, 1]),
    ])

csf[(0, 0, 2, 3, 4)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [60, 60]),
    ([2, 1, 1, 1], [27, 27]),
    ([2, 2, 1], [12, 12]),
    ([3, 1, 1], [7, 7]),
    ([3, 2], [3, 3]),
    ([4, 1], [1, 1]),
    ])

csf[(0, 1, 1, 1, 1)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [5, 15, 25, 30, 25, 15, 5]),
    ([2, 1, 1, 1], [1, 3, 5, 6, 5, 3, 1]),
    ])

csf[(0, 1, 1, 1, 2)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [5, 15, 40, 40, 15, 5]),
    ([2, 1, 1, 1], [1, 3, 11, 11, 3, 1]),
    ([2, 2, 1], [0, 0, 2, 2]),
    ([3, 1, 1], [0, 0, 1, 1]),
    ])

csf[(0, 1, 1, 1, 3)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [5, 30, 50, 30, 5]),
    ([2, 1, 1, 1], [1, 9, 16, 9, 1]),
    ([2, 2, 1], [0, 2, 4, 2]),
    ([3, 1, 1], [0, 1, 2, 1]),
    ])

csf[(0, 1, 1, 1, 4)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [20, 40, 40, 20]),
    ([2, 1, 1, 1], [7, 14, 14, 7]),
    ([2, 2, 1], [2, 4, 4, 2]),
    ([3, 1, 1], [1, 2, 2, 1]),
    ])

csf[(0, 1, 1, 2, 2)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [5, 30, 50, 30, 5]),
    ([2, 1, 1, 1], [1, 9, 16, 9, 1]),
    ([2, 2, 1], [0, 2, 4, 2]),
    ([3, 1, 1], [0, 1, 2, 1]),
    ])

csf[(0, 1, 1, 2, 3)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [5, 55, 55, 5]),
    ([2, 1, 1, 1], [1, 20, 20, 1]),
    ([2, 2, 1], [0, 7, 7]),
    ([3, 1, 1], [0, 3, 3]),
    ([3, 2], [0, 1, 1]),
    ])

csf[(0, 1, 1, 2, 4)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [20, 80, 20]),
    ([2, 1, 1, 1], [7, 34, 7]),
    ([2, 2, 1], [2, 14, 2]),
    ([3, 1, 1], [1, 8, 1]),
    ([3, 2], [0, 3]),
    ([4, 1], [0, 1]),
    ])

csf[(0, 1, 1, 3, 3)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [30, 60, 30]),
    ([2, 1, 1, 1], [12, 24, 12]),
    ([2, 2, 1], [5, 10, 5]),
    ([3, 1, 1], [2, 4, 2]),
    ([3, 2], [1, 2, 1]),
    ])

csf[(0, 1, 1, 3, 4)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [60, 60]),
    ([2, 1, 1, 1], [27, 27]),
    ([2, 2, 1], [12, 12]),
    ([3, 1, 1], [7, 7]),
    ([3, 2], [3, 3]),
    ([4, 1], [1, 1]),
    ])

csf[(0, 1, 2, 2, 2)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [20, 40, 40, 20]),
    ([2, 1, 1, 1], [7, 14, 14, 7]),
    ([2, 2, 1], [2, 4, 4, 2]),
    ([3, 1, 1], [1, 2, 2, 1]),
    ])

csf[(0, 1, 2, 2, 3)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [20, 80, 20]),
    ([2, 1, 1, 1], [7, 34, 7]),
    ([2, 2, 1], [2, 14, 2]),
    ([3, 1, 1], [1, 8, 1]),
    ([3, 2], [0, 3]),
    ([4, 1], [0, 1]),
    ])

csf[(0, 1, 2, 2, 4)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [60, 60]),
    ([2, 1, 1, 1], [27, 27]),
    ([2, 2, 1], [12, 12]),
    ([3, 1, 1], [7, 7]),
    ([3, 2], [3, 3]),
    ([4, 1], [1, 1]),
    ])

csf[(0, 1, 2, 3, 3)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [60, 60]),
    ([2, 1, 1, 1], [27, 27]),
    ([2, 2, 1], [12, 12]),
    ([3, 1, 1], [7, 7]),
    ([3, 2], [3, 3]),
    ([4, 1], [1, 1]),
    ])

csf[(0, 1, 2, 3, 4)] = m.sum(
    m.term(Partition(index), R(coeffs))
    for index, coeffs in [
    ([1, 1, 1, 1, 1], [120]),
    ([2, 1, 1, 1], [60]),
    ([2, 2, 1], [30]),
    ([3, 1, 1], [20]),
    ([3, 2], [10]),
    ([4, 1], [5]),
    ([5], [1]),
    ])

