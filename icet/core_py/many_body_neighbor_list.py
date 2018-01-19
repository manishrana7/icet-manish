import numpy as np

from icet.core_py.lattice_site import LatticeSite
import copy
class ManyBodyNeighborList(object):

    """
    This is a Python implementation of icet's many-body neighbor list.
    This is used both as a tester of the C++ version and also for
    trying out new functionalities.

    Functionalities are similar to ASE's neigbhor list but are extended
    for more connections than for pairs.

    `bothways=False` means that when collecting all neighbors over
    all indices there will not be any duplicates. In other words, in
    total the algorithm will generate `i,j,k` but not `i,k,j`. To make
    sure of this is the case the generated indices are sorted such
    that `i < j < k`.

    `bothways=True` causes the algorithm to return every neighbor of
    an index.  However, it will not return both `i,j,k` and
    `i,k,j`. In contrast to `bothways=False` the following is allowed:
    `i>j`, `i>k` etc.. (but always `j<k`).
    """

    def __init__(self):
        pass

    def build(self, neighbor_lists, index, bothways=False):
        """
        Will take the neighbor list object and combine the neighbors
        of `index` up to `order`.

        Parameters
        ----------
        neighbor_lists : list of ASE NeighborList objects
            ASE neighbor lists
        index : int
            index of site for which to return neighbors
        bothways : boolean
            `False` will return all indices that are bigger than `index`;
            `True` will also return indices that are smaller.
        """
        if not isinstance(neighbor_lists, list):
            neighbor_lists = [neighbor_lists]

        if not neighbor_lists:
            raise RuntimeError(
                'neighbor_lists empty in TestManyBodyNeighborList::build')

        many_body_neighbor_indices = []

        self.add_singlet(index, many_body_neighbor_indices)

        self.add_pairs(index, neighbor_lists[0],
                       many_body_neighbor_indices, bothways)

        """ Add neighbors of higher order (k>=2) """
        for k in range(2, len(neighbor_lists) + 2):

            """ Get neighbors of index in icet format """
            neighbor = self.get_neighbor_from_nl(neighbor_lists[k - 2], index)

            zero_vector = np.array([0., 0., 0., ])

            current_original_neighbors = [LatticeSite(index, zero_vector)]

            self.combine_to_higher_order(
                neighbor_lists[k - 2], many_body_neighbor_indices, neighbor,
                current_original_neighbors, bothways, k)

        return many_body_neighbor_indices

    def combine_to_higher_order(self, nl, many_body_neighbor_indices,
                                neighbor_i, current_original_neighbors,
                                bothways, order):
        """
        For each `j` in `neighbor` construct
        the intersect of `neighbor_j` and `neighbor`,
        call the intersect `neighbor_ij`. All
        neighbors in `neighbor_ij` are then
        neighbors of `i` and `j` What is saved
        then is `(i,j)` and `neighbor_ij` up
        to the desired `order`.

        Parameters
        ----------
        nl : ASE neighbor_list object
        many_body_neighbor_indices: list
            neighbor_lists, each inner list is made up.
        neighbor_i : list
            neighbors of a chosen index in the icet format [[index,offset]]
        current_original_neighbors :
        bothways : boolean
            False will return all indices that are bigger than site "index"
            True will return also return indices that are smaller.
        order : int
            highest order for many-body neighbor indices.
        """
        for j in neighbor_i:

            original_neighbor_copy = copy.deepcopy( current_original_neighbors )

            # explain this line here
            if (len(original_neighbor_copy) == 1 and
                (not bothways and
                 j < original_neighbor_copy[-1])):
                continue

            original_neighbor_copy.append(j)

            neighbor_j_offset = self.translate_all_neighbor(
                self.get_neighbor_from_nl(nl, j.index), j.unitcell_offset)

            if not bothways:
                neighbor_j_offset = self.filter_neighbor_from_smaller(
                    neighbor_j_offset, j)

            intersection_with_j = self.get_intersection(
                neighbor_i, neighbor_j_offset)

            if len(original_neighbor_copy) + 1 < order:
                self.combine_to_higher_order(
                    nl, many_body_neighbor_indices, intersection_with_j,
                    original_neighbor_copy, bothways, order)

            if (len(intersection_with_j) > 0 and
                    len(original_neighbor_copy) == (order - 1)):
                many_body_neighbor_indices.append(
                    [original_neighbor_copy, intersection_with_j])

    def add_singlet(self, index, mbn_indices):
        """
        Add singlet to many-body neighbor indices
        """
        offset = [0., 0., 0.]        
        mbn_indices.append(LatticeSite(index, offset))

    def add_pairs(self, index, neigbhorlist, mbn_indices, bothways):
        """
        Add pairs from neigbhorlist to many-body neighbor indices
        """
        offset = [0., 0., 0.]
        first_site = LatticeSite(index, offset)
        neighbor = self.get_neighbor_from_nl(neigbhorlist, index)
        if not bothways:
            neighbor = self.filter_neighbor_from_smaller(neighbor, first_site)
        if len(neighbor) == 0:
            return
        mbn_indices.append([first_site, neighbor])

    def get_intersection(self, neighbor_i, neighbor_j):
        """
        Return intersection of neighbor_i with neighbor_j
        using the is_j_in_neighbor bool method.
        """
        neighbor_ij = []
        for j in neighbor_i:
            if self.is_j_in_neighbor(j, neighbor_j):
                neighbor_ij.append(j)
        return neighbor_ij

    def is_j_in_neighbor(self, j, neighbor):
        """
        Returns true if there is an index in neighbor that is equal to j
        """
        return j in neighbor
            
        for k in neighbor:
            if k == j:
                return True
        return False

    def filter_neighbor_from_smaller(self, neighbor_i, j):
        """
        Returns all k in neighbor_i that are bigger than j
        """
        neighbor_j_filtered = []
        for k in neighbor_i:
             if j <  k:
                neighbor_j_filtered.append(k)
        return neighbor_j_filtered

    def translate_all_neighbor(self, neighbor, offset):
        """
        Make a copy of neighbor and returns neighbor
        but with all offset "offseted" an
        addition offset
        """
        neighbor_i_offset = neighbor.copy()
        for j in neighbor_i_offset:
            j.unitcell_offset += offset
        return neighbor_i_offset

    def compare_arrays(self, arr1, arr2):
        """
        Compares two arrays.
        Compares element by element in order.
        Returns true if arr1 < arr2
        """
        assert len(arr1) == len(arr2)
        for i in range(len(arr1)):
            if arr1[i] < arr2[i]:
                return True
            if arr1[i] > arr2[i]:
                return False
        return False


    def get_neighbor_from_nl(self, ase_nl, index):
        """
        Get the neighbors of index in the format used in icet
        """
        indices, offsets = ase_nl.get_neighbors(index)
        neighbor = []
        for ind, offs in zip(indices.copy(), offsets.copy()):
            neighbor.append(LatticeSite(ind, offs))
        return neighbor
