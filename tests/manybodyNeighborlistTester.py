from __future__ import print_function, division
from ase import Atoms
from ase.neighborlist import NeighborList
from ase.build import bulk
import numpy as np

"""
This is a python implementation of icet's  manybody neighborlist.
This is used both as a tester of the c++ version and also a venue
to try out new functionalities.


Functionalities are similar to ASE neigbhorlist but is extended
for more connections than for pairs.

Bothways = False will mean that if you collect all neighbors over all indices
then you will not have any duplicates. I.e. in total you will only generate
i,j,k and never i,k,j.
To make sure of this the indices you generate will always be sorted that is:
i < j < k  etc

Bothways = True will mean that you will return every neighbor of an indice.
However you will not return both i,j,k and i,k,j.
Contrary to bothways = False this is allowed: i>j, i>k etc.. (but always j < k)


"""


class manybodyNeighborlistTester():
    def __init__(self):
        self.initiated = True

    def build(self, nl, index, order, bothways=False):
        """
        Will take the neighborlist object (nl) and combine the neighbors
        of index "index" up to order "order".

        Params:
        nl : ASE neighborlist
        index: index to return neighbors from
        order: order or how many sites a neighbir is
        bothways: see above...
        """

        Ni = self.get_Ni_from_nl(nl, index)
        numberOfSites = len(nl.positions)
        c = 2
        zero_vector = np.array([0., 0., 0., ])
        current_original_neighbors = [[index, zero_vector]]
        manybody_neighbor_indices = []
        self.combine_to_higher_order(
            nl, manybody_neighbor_indices, Ni, current_original_neighbors, c, bothways, order)
        return manybody_neighbor_indices

    def combine_to_higher_order(self, nl, manybody_neighbor_indices, Ni, current_original_neighbors, c, bothways, order):
        for j in Ni:
            originalNeighborCopy = current_original_neighbors.copy()

            if self.nbrCompare(j, originalNeighborCopy[-1]) and not (bothways and len(originalNeighborCopy) == 1):
                continue

            originalNeighborCopy.append(j)

            N_j_offset = self.translate_all_Ni(
                self.get_Ni_from_nl(nl, j[0]), j[1])
            N_j_filtered = self.filter_Ni_from_smaller(N_j_offset, j)
            intersection_ij = self.get_intersection(Ni, N_j_filtered)

            if len(originalNeighborCopy) + 1 < order:
                self.combine_to_higher_order(
                    nl, manybody_neighbor_indices, intersection_ij, originalNeighborCopy, c + 1, bothways, order)
            if len(intersection_ij) > 0:        
                manybody_neighbor_indices.append(
                    [originalNeighborCopy, intersection_ij])

    def get_intersection(self, Ni, Nj):
        """
        Return intersection of Ni and N_j using the is_j_in_Ni bool method
        """
        N_ij = []
        for j in Ni:
            if self.is_j_in_Ni(j, Nj):
                N_ij.append(j)
        return N_ij

    def is_j_in_Ni(self, j, Ni):
        """
        Returns true if there is a k in Ni that is equal to j
        """
        for k in Ni:
            if k[0] == j[0] and (k[1] == j[1]).all():
                return True

        return False

    def filter_Ni_from_smaller(self, N_i, j):
        """
        Returns all k in N_i that are bigger than j
        """
        N_j_filtered = []
        for k in N_i:
            if not self.nbrCompare(k, k):
                N_j_filtered.append(k)
        return N_j_filtered

    def translate_all_Ni(self, Ni, offset):
        """
        Make a copy of Ni and returns Ni but with all offset "offseted" an addition offset
        """
        N_i_offset = Ni.copy()
        for j in N_i_offset:
            j[1] += offset
        return N_i_offset

    def arrayCompare(self, arr1, arr2):
        """
        Compares two arrays.
        Compares element by element in order.
        Returns true if arr1 < arr2
        """
        assert len(arr1) == len(arr2)
        for i in range(3):
            if arr1[i] < arr2[i]:
                return True
            if arr1[i] > arr2[i]:
                return False
        return False

    def nbrCompare(self, nbr1, nbr2):
        """
        Compare two neighbors as defined in this class
        returns True if nbr1 < nbr2
        returns False otherwise
        """
        if nbr1[0] < nbr2[0]:
            return True
        if nbr1[0] > nbr2[0]:
            return False
        return self.arrayCompare(nbr1[1], nbr2[1])

    def get_Ni_from_nl(self, ase_nl, index):
        """
        Get the neighbors of index in the format used in icet
        """
        indices, offsets = ase_nl.get_neighbors(index)
        Ni = []
        for ind, offs in zip(indices, offsets):
            Ni.append([ind, offs])
        return Ni

        def nbrCompare(self, nbr1, nbr2):
            """
            Compares two neighbors.
            return nbr1 < nbr2 by first looking at index
            and if they are equal look at the offsets.
            """
            if nbr1[0] < nbr2[0]:
                return True
            if nbr1[0] > nbr2[0]:
                return False
            return arrayCompare(nbr1[1], nbr2[1])


"""
def naiveManybodyThirdOrder(nl, index, bothways=True):

    nbr_0 = nl.get_neighbors(index)
    nbr_index = (index, [0, 0, 0])
    nbrs = []
    for j in nbr_0:
        if not bothways and nbrCompare(j, nbr_index):
            continue
        nbr_j = nl.get_neighbors(j[0])
        for k in nbr_j:
            if (nl.is_neighbor(index, k[0], k[1] + j[1])):
                neighbor_k = (k[0], k[1] + j[1])
                manybodyNbr = [index, [0., 0., 0.], j, neighbor_k]
                if not bothways:
                    if nbrCompare(j, neighbor_k):
                        nbrs.append(manybodyNbr)
                else:
                    nbrs.append(manybodyNbr)
    return nbrs
"""