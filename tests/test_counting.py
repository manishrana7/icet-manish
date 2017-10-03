from icetdev import orbitList
from icetdev.structure import structure_from_atoms, Structure
from icetdev.manybodyNeighborlist import *
from ase import Atoms
from ase.build import bulk
from icetdev.neighborlist import get_neighborlists
from icetdev.orbitList import create_orbit_list
from ase.spacegroup import crystal
import time
from icetdev.clusterCounts import *
import random


def test_no_symmetry_local_orbitlist_counting(prim_atoms, cutoffs, repeatInteger):
    """
    Creates a primitive orbitlist and create an orbitlist on supercell and
    then compares the counting over both of them.
    """
    atoms = prim_atoms.copy()
    atoms_repeat = prim_atoms.copy().repeat(repeatInteger)
    for atom in atoms_repeat:
        if random.random() < 0.5:
            atom.symbol = "H"

    """ Set up neighborlists """
    neighborlists_prim = get_neighborlists(atoms=atoms, cutoffs=cutoffs)
    structure = structure_from_atoms(atoms)
    structure_repeat = structure_from_atoms(atoms_repeat)
    neighborlists_supercell = get_neighborlists(
        atoms=atoms_repeat, cutoffs=cutoffs)

    ##########################################################################
    #                                                                        #
    #   set up orbitlists. Important to sort so clusters get same hash       #
    #                                                                        #
    ##########################################################################

    """ Orbitlist primitive """
    prim_orbitlist = orbitList.OrbitList(neighborlists_prim, structure)
    prim_orbitlist.sort()

    """ orbitlist supercell"""
    supercell_orbitlist = orbitList.OrbitList(
        neighborlists_supercell, structure_repeat)
    supercell_orbitlist.sort()

    ################################################################
    #                                                              #
    #             set up clustercounts and count clusters          #
    #                                                              #
    ################################################################

    clusterCount_local = ClusterCounts()

    clusterCount_local.count_each_local_orbitlist(
        structure_repeat, prim_orbitlist)

    clusterCount_supercell = ClusterCounts()

    clusterCount_supercell.count_orbitlist(
        structure_repeat, supercell_orbitlist)

    # Get the clustercount map
    local_cluster_map = clusterCount_local.get_cluster_counts()
    supercell_cluster_map = clusterCount_supercell.get_cluster_counts()

    assert len(local_cluster_map) == len(
        supercell_cluster_map), "lengths of cluster counts in test_no_symmetry_local_orbitlist_counting is not same"

    assert local_cluster_map.keys() == supercell_cluster_map.keys()

    for key in local_cluster_map.keys():
        assert local_cluster_map[key] == supercell_cluster_map[key]


def get_total_count(cluster_count_dict, cluster):
    """
    Returns the total count of a particular cluster
    """
    count = 0
    for counts in cluster_count_dict[cluster]:
        count += cluster_count_dict[cluster][counts]
    return count


def test_no_symmetry_vs_symmetry_count(prim_atoms, cutoffs, repeatInteger):
    """
    Here we use a cutoff so that no extra clusters are found in the symmetry case
    and compare the counts found in both methods (should be equal)
    """
    atoms = prim_atoms.copy()
    atoms_repeat = prim_atoms.copy().repeat(repeatInteger)
    for atom in atoms_repeat:
        if random.random() < 0.3:
            atom.symbol = "H"

    structure_repeat = structure_from_atoms(atoms_repeat)
    structure = structure_from_atoms(atoms)

    ####################################
    #                                  #
    #      no symmetry counting        #
    #                                  #
    ####################################

    # get neighborlist
    neighborlists = get_neighborlists(atoms=atoms, cutoffs=cutoffs)

    # get orbitlist
    orbitlist_no_symmetry = orbitList.OrbitList(neighborlists, structure)
    # setup cluster count and count:
    clusterCount_no_symmetry = ClusterCounts()
    clusterCount_no_symmetry.count_each_local_orbitlist(
        structure_repeat, orbitlist_no_symmetry)

    # get the clustercount map
    clusterCountMap_no_symmetry = clusterCount_no_symmetry.get_cluster_counts()

    ##################################
    #                                #
    #     counting with symmetry     #
    #                                #
    ##################################

    # get orbitlist
    orbitlist_symmetry = create_orbit_list(structure, cutoffs, verbosity=0)

    # setup clustercount and count
    clustercounts_symmetry = ClusterCounts()
    clustercounts_symmetry.count_each_local_orbitlist(
        structure_repeat, orbitlist_symmetry)

    # get the clustercount map
    clusterCountMap_symmetry = clustercounts_symmetry.get_cluster_counts()

    assert orbitlist_symmetry.size() == orbitlist_no_symmetry.size(
    ), "test that orbitlists are of equal size"

    for key in clusterCountMap_no_symmetry.keys():
        assert get_total_count(clusterCountMap_no_symmetry, key) == get_total_count(
            clusterCountMap_symmetry, key), "multiplicity of cluster should be equal"

        for element_key in clusterCountMap_no_symmetry[key]:
            assert clusterCountMap_no_symmetry[key][element_key] == clusterCountMap_symmetry[key][
                element_key], " test that for each cluster that the count for each element combination is the same"


atoms = bulk("Al", "bcc", a=1)
cutoffs = [1.61] * 5
N = 6


test_no_symmetry_local_orbitlist_counting(atoms, cutoffs, N)
test_no_symmetry_vs_symmetry_count(atoms, cutoffs, N)
