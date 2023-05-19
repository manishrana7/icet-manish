#!/usr/bin/env python3
import unittest
import numpy as np
from icet.tools.structure_generation import occupy_structure_randomly
from icet.tools.training_set_generation import structure_annealing
from icet import ClusterSpace
from ase.build import bulk


class TestStructurePoolGenerationWithAnnealing(unittest.TestCase):
    """
    Container for tests of the class functionality
    """

    def __init__(self, *args, **kwargs):
        super(TestStructurePoolGenerationWithAnnealing, self).__init__(*args, **kwargs)
        prim = bulk('Au', a=4.0)
        self.cs = ClusterSpace(prim, [6.0], [['Au', 'Pd']])
        target_concentrations = []
        for i in range(1, 10):
            target_concentrations.append({'Au': i / 10, 'Pd': (10 - i) / 10})
        self.structure_pool = []
        for _ in range(int(5e2)):
            structure = prim.copy().repeat(10)
            indx = np.random.choice(range(9))
            concentration = target_concentrations[indx]
            occupy_structure_randomly(structure, self.cs,
                                      concentration)
            self.structure_pool.append(structure)

    def test_structure_annealing_correct_steps_and_number_of_structures(self):
        """
        Test that the number of structure is correct for a few cases
        """
        (inds, traj) = structure_annealing(self.cs,
                                           self.structure_pool,
                                           10,
                                           int(1e1))
        assert len(inds) == 10

        initial_indices = [f for f in range(10)]
        (inds, traj) = structure_annealing(self.cs,
                                           self.structure_pool,
                                           10,
                                           int(1e1),
                                           initial_indices=initial_indices)
        assert len(inds) == 10

        base = [self.structure_pool[i] for i in initial_indices]
        tmp_structure_pool = self.structure_pool[5:]
        (inds, traj) = structure_annealing(self.cs,
                                           tmp_structure_pool,
                                           10,
                                           int(1e1),
                                           base_structures=base,
                                           initial_indices=initial_indices)
        assert len(inds) == 10

        base = [self.structure_pool[i] for i in initial_indices]
        tmp_structure_pool = self.structure_pool[5:]
        (inds, traj) = structure_annealing(self.cs,
                                           tmp_structure_pool,
                                           10,
                                           int(1e1),
                                           base_structures=base)
        assert len(inds) == 10

    def test_structure_annealing_with_initial_indices(self):
        """
        Test that it start with the correct indices
        """
        initial_indices = [f for f in range(10)]
        (inds, _) = structure_annealing(self.cs,
                                        self.structure_pool,
                                        10,
                                        0,
                                        initial_indices=initial_indices)
        assert np.all(np.array(initial_indices) == inds)

    def test_structure_annealing_minimization_condition_number(self):
        """
        Test that the condition number is minimized
        """
        (_, traj) = structure_annealing(self.cs,
                                        self.structure_pool,
                                        10,
                                        int(1e3))
        assert traj[-1] < traj[0]

    def test_structure_annealing_minimization_correlation(self):
        """
        Test that the correlation is minimized
        """
        (_, traj) = structure_annealing(self.cs,
                                        self.structure_pool,
                                        10,
                                        int(1e3),
                                        metric_function='covariance')
        assert traj[-1] < traj[0]

    def test_structure_annealing_check_correlation_minimize_condition(self):
        """
        Test that the correlation also minimizes the condition number
        """
        initial_indices = [f for f in range(10)]
        (end_inds, traj) = structure_annealing(self.cs,
                                               self.structure_pool,
                                               10,
                                               int(1e3),
                                               initial_indices=initial_indices)
        A_start = []
        for ind in initial_indices:
            A_start.append(self.cs.get_cluster_vector(self.structure_pool[ind]))
        A_start = np.array(A_start)

        A_end = []
        for ind in end_inds:
            A_end.append(self.cs.get_cluster_vector(self.structure_pool[ind]))
        A_end = np.array(A_end)

        assert np.linalg.cond(A_end) < np.linalg.cond(A_start)

    def test_structure_annealing_gives_correct_last_indices(self):
        """
        Test that we return the correct structure indices
        """
        initial_indices = [f for f in range(10)]
        (end_inds, traj) = structure_annealing(self.cs,
                                               self.structure_pool,
                                               10,
                                               int(1e3),
                                               initial_indices=initial_indices)
        A_end = []
        for ind in end_inds:
            A_end.append(self.cs.get_cluster_vector(self.structure_pool[ind]))
        A_end = np.array(A_end)
        assert np.linalg.cond(A_end) == traj[-1]

    def test_structure_annealing_with_base_input(self):
        """
        Test that we return the correct structure indices
        """
        initial_indices = [i for i in range(5)]
        base = [self.structure_pool[i] for i in initial_indices]
        tmp_structure_pool = self.structure_pool[5:]
        (end_inds, traj) = structure_annealing(self.cs,
                                               tmp_structure_pool,
                                               10,
                                               int(1e3),
                                               base_structures=base)
        A = [self.cs.get_cluster_vector(base[i]) for i in range(5)]
        for ind in end_inds:
            A.append(self.cs.get_cluster_vector(tmp_structure_pool[ind]))
        A = np.array(A)
        assert np.linalg.cond(A) == traj[-1]


if __name__ == '__main__':
    unittest.main()
