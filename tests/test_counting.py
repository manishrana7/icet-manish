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


def setup_test_orbitlist(atoms, cutoffs):
    """
    A helper for getting structure, mbnl and neighborlists used for setting up a orbitlist without symmetry
    """
    structure = structure_from_atoms(atoms)
    mbnl = ManybodyNeighborlist()
    neighborlists = get_neighborlists(atoms=atoms, cutoffs=cutoffs)
    mbnl.build(neighborlists, 0, True)
    return structure, mbnl, neighborlists


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

    # atoms_repeat = atoms.repeat(repeatInteger)

    structure_repeat = structure_from_atoms(atoms_repeat)
    structure, mbnl, neighborlists = setup_test_orbitlist(atoms, cutoffs)

    prim_orbitlist = orbitList.OrbitList(neighborlists, structure)
    clusterCount_local = ClusterCounts()

    clusterCount_local.count_each_local_orbitlist(
        structure_repeat, prim_orbitlist)

    """ Repeat for supercell"""
    structure_repeat, mbnl, neighborlists = setup_test_orbitlist(
        atoms_repeat, cutoffs)

    supercell_orbitlist = orbitList.OrbitList(neighborlists, structure_repeat)
    clusterCount_supercell = ClusterCounts()

    clusterCount_supercell.count_orbitlist(
        structure_repeat, supercell_orbitlist)

    local_cluster_map = clusterCount_local.get_cluster_counts()

    supercell_cluster_map = clusterCount_supercell.get_cluster_counts()

    assert len(local_cluster_map) == len(
        supercell_cluster_map), "lengths of cluster counts in test_no_symmetry_local_orbitlist_counting is not same"

    assert local_cluster_map.keys() == supercell_cluster_map.keys()

    for key in local_cluster_map.keys():
        assert local_cluster_map[key] == supercell_cluster_map[key]
        # print(local_cluster_map[key])
        # print(supercell_cluster_map[key])
        # print(" ")


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

    # no symmetry counting
    structure, mbnl, neighborlists = setup_test_orbitlist(atoms, cutoffs)
    orbitlist_no_symmetry = orbitList.OrbitList(neighborlists, structure)
    clusterCount_no_symmetry = ClusterCounts()
    clusterCount_no_symmetry.count_each_local_orbitlist(
        structure_repeat, orbitlist_no_symmetry)
    clusterCountMap_no_symmetry = clusterCount_no_symmetry.get_cluster_counts()

    # counting with symmetry
    orbitlist_symmetry = create_orbit_list(structure, cutoffs, verbosity=0)
    clustercounts_symmetry = ClusterCounts()
    clustercounts_symmetry.count_each_local_orbitlist(
        structure_repeat, orbitlist_symmetry)
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
