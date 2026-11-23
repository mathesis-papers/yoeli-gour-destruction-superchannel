#!/usr/bin/env python3
"""Exact finite witnesses for dephasing-covariant fixed-channel rigidity."""
import json
from itertools import permutations
import sympy as algebra


def equal(actual, expected):
    assert actual.shape == expected.shape
    assert all(algebra.expand(entry) == 0 for entry in actual - expected)


def matrix_unit(dimension, row, column):
    result = algebra.zeros(dimension)
    result[row, column] = 1
    return result


def dephaser(dimension):
    return algebra.diag(*[int(row == column)
                          for row in range(dimension)
                          for column in range(dimension)])


def monomial(permutation, phase_index):
    dimension = len(permutation)
    return algebra.Matrix(dimension, dimension, lambda row, column:
                          (algebra.I if column == phase_index else 1)
                          if row == permutation[column] else 0)


def vectorize(matrix):
    return algebra.Matrix([matrix[row, column]
                           for column in range(matrix.cols)
                           for row in range(matrix.rows)])


def channel_matrix(unitary):
    return algebra.kronecker_product(unitary, unitary.conjugate())


def choi(channel, dimension):
    return algebra.Matrix(dimension ** 2, dimension ** 2, lambda row, column:
                          channel[(row % dimension) * dimension + column % dimension,
                                  (row // dimension) * dimension + column // dimension])


def check_family(dimension, permutation_list):
    delta = dephaser(dimension)
    vectors = []
    for permutation in permutation_list:
        base = vectorize(monomial(permutation, None))
        for phase_index in [None] + list(range(dimension)):
            unitary = monomial(permutation, phase_index)
            equal(unitary.H * unitary, algebra.eye(dimension))
            channel = channel_matrix(unitary)
            equal(channel * delta, delta * channel)
            vector = vectorize(unitary)
            equal(choi(channel, dimension), vector * vector.H)
            assert (vector.H * vector)[0] == dimension
            if phase_index is not None:
                expected = algebra.zeros(dimension ** 2, 1)
                expected[phase_index * dimension + permutation[phase_index]] = algebra.I - 1
                equal(vector - base, expected)
            vectors.append(vector)
    span_rank = algebra.Matrix.hstack(*vectors).rank()
    assert span_rank == dimension ** 2
    remaining = set(range(len(vectors)))
    components = []
    while remaining:
        frontier = [remaining.pop()]
        component_size = 0
        while frontier:
            vertex = frontier.pop()
            component_size += 1
            adjacent = {other for other in remaining
                        if algebra.expand((vectors[vertex].H * vectors[other])[0]) != 0}
            remaining -= adjacent
            frontier.extend(adjacent)
        components.append(component_size)
    assert len(components) == (2 if dimension == 2 else 1)
    rotation = algebra.eye(dimension)
    rotation[:2, :2] = algebra.Matrix([[algebra.Rational(3, 5), -algebra.Rational(4, 5)],
                                     [algebra.Rational(4, 5), algebra.Rational(3, 5)]])
    equal(rotation.H * rotation, algebra.eye(dimension))
    assert (rotation * matrix_unit(dimension, 0, 0) * rotation.H)[0, 1] == algebra.Rational(12, 25)
    assert channel_matrix(rotation) * delta != delta * channel_matrix(rotation)
    return {'dimension': dimension, 'unitary_vectors': len(vectors),
            'span_rank': span_rank, 'overlap_component_sizes': sorted(components),
            'rotation_off_diagonal': '12/25'}


def check_qubit():
    phase_flip = channel_matrix(algebra.diag(1, -1))
    delta = dephaser(2)
    equal(delta, (algebra.eye(4) + phase_flip) / 2)
    theta_mask = algebra.Matrix(4, 4, lambda row, column:
                               (1 + phase_flip[row, row] * phase_flip[column, column]) / 2)
    covariance_mask = algebra.Matrix(4, 4, lambda row, column:
                                    int(delta[row, row] == delta[column, column]))
    equal(theta_mask, covariance_mask)
    for row in range(4):
        for column in range(4):
            elementary = matrix_unit(4, row, column)
            image = (elementary + phase_flip * elementary * phase_flip) / 2
            equal((image + phase_flip * image * phase_flip) / 2, image)
    return {'fixed_space_mask': 'equal to covariance mask', 'idempotence_basis_units': 16}


def check_source_construction():
    for dimension in (1, 2, 3):
        delta = dephaser(dimension)
        equal(delta * delta, delta)
        identity_vector = vectorize(algebra.eye(dimension))
        depolarizing = identity_vector * identity_vector.T / dimension
        equal(choi(depolarizing, dimension), algebra.eye(dimension ** 2) / dimension)
        equal(delta * depolarizing, depolarizing * delta)
        for row in range(dimension ** 2):
            for column in range(dimension ** 2):
                elementary = matrix_unit(dimension ** 2, row, column)
                projected = delta * elementary * delta
                equal(delta * projected * delta, projected)
                equal(choi(projected, dimension), algebra.diag(*choi(elementary, dimension).diagonal()))
    return {'dimensions': [1, 2, 3], 'two_sided_choi_action': 'complete dephasing'}


def check_uniformity():
    dimension = 2
    uniform_vector = vectorize(algebra.eye(dimension))
    depolarizing = uniform_vector * uniform_vector.T / dimension
    unitaries = [monomial(permutation, phase_index)
                 for permutation in permutations(range(dimension))
                 for phase_index in [None, 0, 1]]
    rotation = algebra.Matrix([[algebra.Rational(3, 5), -algebra.Rational(4, 5)],
                               [algebra.Rational(4, 5), algebra.Rational(3, 5)]])
    unitaries.append(rotation)
    vectors = [vectorize(unitary) for unitary in unitaries]
    for unitary in unitaries:
        equal(channel_matrix(unitary) * depolarizing, depolarizing)
        equal(depolarizing * channel_matrix(unitary), depolarizing)
    reached = {0}
    while True:
        adjacent = {index for index, vector in enumerate(vectors)
                    if any(algebra.expand((vectors[other].H * vector)[0]) != 0
                           for other in reached)}
        if adjacent <= reached:
            break
        reached |= adjacent
    assert len(reached) == 7
    assert algebra.Matrix.hstack(*vectors).rank() == 4
    equal((vectors[0].H * vectors[-1]) / 2, algebra.Matrix([[algebra.Rational(3, 5)]]))
    equal((vectors[-1].H * vectors[4]) / 2,
          algebra.Matrix([[2 * (algebra.I - 1) / 5]]))
    reset = vectorize(matrix_unit(2, 0, 0)) * uniform_vector.T
    assert reset * depolarizing != depolarizing * reset
    equal(choi(reset, 2), algebra.kronecker_product(algebra.eye(2), matrix_unit(2, 0, 0)))
    return {'dimension': 2, 'unitary_vectors': 7, 'span_rank': 4,
            'overlap_component_sizes': [7], 'nonunital_witness': 'reset to E_00'}


def check_two_sector_pinching():
    cases = []
    for dimension, rank in [(2, 1), (3, 1), (4, 1), (4, 2)]:
        involution = algebra.diag(*([1] * rank + [-1] * (dimension - rank)))
        conjugation = channel_matrix(involution)
        pinching = (algebra.eye(dimension ** 2) + conjugation) / 2
        equal(pinching * pinching, pinching)
        for row in range(dimension ** 2):
            for column in range(dimension ** 2):
                equal_mask = int(pinching[row, row] == pinching[column, column])
                twirl_mask = (1 + conjugation[row, row] * conjugation[column, column]) / 2
                assert equal_mask == twirl_mask
                assert twirl_mask ** 2 == twirl_mask
        cases.append({'dimension': dimension, 'first_sector_rank': rank,
                      'coefficient_units': dimension ** 4})
    return {'fixed_space_equals_covariance': True, 'cases': cases}


def check_algebraic_projection():
    dimension = 3
    householder = algebra.Matrix([[7, -4, -4], [-4, 1, -8], [-4, -8, 1]]) / 9
    equal(householder.H * householder, algebra.eye(dimension))
    channel = channel_matrix(householder)
    delta = dephaser(dimension)
    complement = algebra.eye(dimension ** 2) - delta
    for row in range(dimension ** 2):
        for column in range(dimension ** 2):
            keep = delta[row, row] * delta[column, column] + complement[row, row] * complement[column, column]
            assert keep == int(delta[row, row] == delta[column, column])
            assert keep ** 2 == keep
    projected = delta * channel * delta + complement * channel * complement
    equal(delta * projected, projected * delta)
    equal(delta * projected * delta + complement * projected * complement, projected)
    trace_row = vectorize(algebra.eye(dimension)).T
    equal(trace_row * projected, trace_row)
    principal = choi(projected, dimension).extract([0, 4, 2], [0, 4, 2])
    determinant = algebra.factor(principal.det())
    assert determinant == -algebra.Rational(784, 531441)
    affine_projection = algebra.Matrix([[1, 1, 0], [1, 1, 1], [0, 1, 1]]) / 3
    assert affine_projection.det() == -algebra.Rational(1, 27)
    assert affine_projection.trace() == 1
    return {'householder_unitary': True, 'projection_idempotent': True,
            'projection_coefficient_units': 81,
            'output_covariant_and_trace_preserving': True,
            'choi_principal_indices': ['(0,0)', '(1,1)', '(0,2)'],
            'choi_principal_determinant': str(determinant),
            'affine_state_projection_determinant': '-1/27'}


def check_replacer_and_asymmetry():
    gamma = algebra.diag(algebra.Rational(1, 3), algebra.Rational(2, 3))
    trace_row = vectorize(algebra.eye(2)).T
    replacer = vectorize(gamma) * trace_row
    equal(replacer * replacer, replacer)
    for unitary in [algebra.eye(2), algebra.Matrix([[0, 1], [1, 0]])]:
        equal(replacer * channel_matrix(unitary), replacer)
    swap = channel_matrix(algebra.Matrix([[0, 1], [1, 0]]))
    phase = channel_matrix(algebra.diag(1, algebra.I))
    delta = dephaser(2)
    equal(swap * delta, delta * swap)
    assert swap * phase != phase * swap
    return {'full_rank_replacer': ['1/3', '2/3'],
            'swap_commutes_with_dephasing': True,
            'swap_does_not_commute_with_phase_1_i': True}


def main():
    qutrit_permutations = [(0, 1, 2), (1, 2, 0), (2, 0, 1), (1, 0, 2)]
    result = {'arithmetic': 'exact over Q(i)',
              'qutrit_family': check_family(3, qutrit_permutations),
              'qubit_family': check_family(2, list(permutations(range(2)))),
              'qubit_superchannel': check_qubit(),
              'two_sided_construction': check_source_construction(),
              'uniformity_rigidity_qubit_bridge': check_uniformity(),
              'two_sector_pinching': check_two_sector_pinching(),
              'algebraic_projection': check_algebraic_projection(),
              'replacer_and_asymmetry': check_replacer_and_asymmetry(),
              'scope': 'Finite witnesses; the dimension-uniform statements have written proofs in note.tex.',
              'status': 'all finite checks passed'}
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
