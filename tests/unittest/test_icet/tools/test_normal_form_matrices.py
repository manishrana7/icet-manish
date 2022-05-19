import pytest
from icet.tools.structure_enumeration_support.normal_form_matrices import \
    (HermiteNormalForm,
     SmithNormalForm, get_unique_snfs,
     yield_hermite_normal_forms,
     yield_reduced_hnfs)
import numpy as np

"""
@pytest.fixture
def hnf():
    '''A Hermite Normal Form object.'''
    H = [[1, 0, 0], [1, 2, 0], [0, 0, 2]]
    rotations = [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
    translations = [[0, 0, 0]]
    basis_shifts = [[0, 0, 0]]
    return HermiteNormalForm(H, rotations, translations, basis_shifts)
"""


@pytest.mark.parametrize("H, rotations, translations, basis_shifts", (
    (np.array([[1, 0, 0], [1, 2, 0], [0, 0, 2]]),
     [np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])],
     [[0, 0, 0]],
     [[0, 0, 0]]),),
)
def test_initialize_hnf(H, rotations, translations, basis_shifts):
    """Tests initialization of HermiteNormalForm objects."""
    hnf = HermiteNormalForm(H, rotations, translations, basis_shifts)
    assert np.allclose(hnf.H, H)
    assert isinstance(hnf.snf, SmithNormalForm)
    assert np.allclose(hnf.snf.S, (1, 2, 2))

    assert len(hnf.transformations) == 1
    assert len(hnf.transformations[0]) == 3
    assert np.allclose(hnf.transformations[0][0], [[-1, -1,  0],
                                                   [2,  1,  0],
                                                   [0,  0,  1]])
    assert np.allclose(hnf.transformations[0][1], [0, 0, 0])
    assert np.allclose(hnf.transformations[0][2], [0, 0, 0])


@pytest.mark.parametrize("det, pbc, target_hnfs", (
    (1, [True, True, True], [[[1, 0, 0], [0, 1, 0], [0, 0, 1]]]),
    (2, [True, True, True], [[[1, 0, 0], [0, 1, 0], [0, 0, 2]],
                             [[1, 0, 0], [0, 1, 0], [0, 1, 2]],
                             [[1, 0, 0], [0, 1, 0], [1, 0, 2]],
                             [[1, 0, 0], [0, 1, 0], [1, 1, 2]],
                             [[1, 0, 0], [0, 2, 0], [0, 0, 1]],
                             [[1, 0, 0], [1, 2, 0], [0, 0, 1]],
                             [[2, 0, 0], [0, 1, 0], [0, 0, 1]]]),
    (2, [True, True, False], [[[1, 0, 0], [0, 2, 0], [0, 0, 1]],
                              [[1, 0, 0], [1, 2, 0], [0, 0, 1]],
                              [[2, 0, 0], [0, 1, 0], [0, 0, 1]]]),
    (2, [False, True, True], [[[1, 0, 0], [0, 1, 0], [0, 0, 2]],
                              [[1, 0, 0], [0, 1, 0], [0, 1, 2]],
                              [[1, 0, 0], [0, 2, 0], [0, 0, 1]]]),
    (2, [False, False, True], [[[1, 0, 0], [0, 1, 0], [0, 0, 2]]]),
))
def test_yield_hermite_normal_forms(det, pbc, target_hnfs):
    """Tests that Hermite normal form matrices are properly generated."""
    hnfs = [hnf for hnf in yield_hermite_normal_forms(det, pbc)]
    assert len(hnfs) == len(target_hnfs)
    for hnf, target_hnf in zip(hnfs, target_hnfs):
        assert np.allclose(hnf, target_hnf)


@pytest.mark.parametrize("ncells, rotations, pbc, target_hnf_matrices", (
    (2, [[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
         [[-1, 0, 0], [0, 1, 0], [0, 0, -1]]],
     [True, True, True],
     [[[1, 0, 0], [0, 1, 0], [0, 0, 2]],
      [[1, 0, 0], [0, 1, 0], [0, 1, 2]],
      [[1, 0, 0], [0, 1, 0], [1, 0, 2]],
      [[1, 0, 0], [0, 1, 0], [1, 1, 2]],
      [[1, 0, 0], [0, 2, 0], [0, 0, 1]],
      [[1, 0, 0], [1, 2, 0], [0, 0, 1]],
      [[2, 0, 0], [0, 1, 0], [0, 0, 1]]]),
))
def test_yield_reduced_hnfs(ncells, rotations, pbc, target_hnf_matrices):
    """Tests that matrices are properly excluded based on symmetry."""
    symmetries = {'rotations': rotations,
                  'translations': [[0, 0, 1] for _ in rotations],
                  'basis_shifts': [[0, 0, 0] for _ in rotations]}
    for i, hnf in enumerate(yield_reduced_hnfs(ncells, symmetries, pbc)):
        assert np.allclose(hnf.H, target_hnf_matrices[i])


@pytest.mark.parametrize("H, target_S", (
                         ([[2, 0, 0], [0, 2, 0], [0, 0, 2]], [2, 2, 2]),
                         ([[3, 0, 0], [0, 2, 0], [0, 0, 2]], [1, 2, 6]),
                         ([[1, 0, 0], [0, 2, 0], [0, 0, 4]], [1, 2, 4]),
                         ([[1, 0, 0], [0, 2, 0], [0, 2, 4]], [1, 2, 4]),
                         ([[1, 0, 0], [0, 2, 0], [0, 1, 4]], [1, 1, 8]),
                         ))
def test_initialize_snf(H, target_S):
    """Tests initialization of Smith Normal Form objects."""
    snf = SmithNormalForm(np.array(H))
    assert np.allclose(snf.S, target_S)


def test_add_hnf():
    """Tests that addition of H matrix to SNF works."""
    snf = SmithNormalForm(np.array([[1, 0, 0], [0, 2, 0], [0, 0, 4]]))
    assert len(snf.hnfs) == 0
    snf.add_hnf(HermiteNormalForm(np.array([[1, 0, 0], [1, 2, 0], [0, 0, 2]]),
                                  [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                                  [[0, 0, 0]],
                                  [[0, 0, 0]]))
    assert len(snf.hnfs) == 1


def test_set_group_order():
    """Tests that calculation of "group order"."""
    snf = SmithNormalForm(np.array([[1, 0, 0], [0, 2, 0], [0, 0, 4]]))
    assert snf.group_order is None
    snf.set_group_order()
    assert np.allclose(snf.group_order,
                       [[0, 0, 0], [0, 0, 1], [0, 0, 2], [0, 0, 3],
                        [0, 1, 0], [0, 1, 1], [0, 1, 2], [0, 1, 3]])


def test_get_unique_snfs():
    """Tests generation of unique SNFs based on a list of HNFs"""
    Hs = ([[1, 0, 0], [0, 2, 0], [0, 0, 4]],
          [[1, 0, 0], [0, 2, 0], [0, 2, 4]],
          [[1, 0, 0], [0, 2, 0], [0, 1, 4]])
    hnfs = []
    for H in Hs:
        hnfs.append(HermiteNormalForm(np.array(H),
                    [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                    [[0, 0, 0]],
                    [[0, 0, 0]]))

    snfs = get_unique_snfs(hnfs)
    assert len(snfs) == 2
    assert np.allclose(snfs[0].S, [1, 2, 4])
    assert np.allclose(snfs[1].S, [1, 1, 8])
