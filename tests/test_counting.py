from icetdev import orbitList
from icetdev.structure import structure_from_atoms, Structure
from icetdev.manybodyNeighborlist import *
from ase import Atoms
from ase.db import connect
from icetdev.neighborlist import get_neighborlists
from icetdev.orbitList import create_orbit_list
from ase.spacegroup import crystal
import time
from icetdev.clusterCounts import *
import random
import numpy as np


def get_equivalent_clustermap_key(key1, clustermap_keys, tol=1e-3):
    """
    Will search for key1 in clustermap keys such that the distances are all roughly the same and
     the geometrical size is rougly the same.
     It will return the key in clustermap_keys that are within this threshhold if none is find
     it will throw an exception
    """

    for key2 in clustermap_keys:
        if len(key1.get_sites()) == len(key2.get_sites()):
            if np.linalg.norm(np.array(key1.get_distances()) - np.array(key2.get_distances())) < tol:
                if np.abs(key1.get_geometrical_size() - key2.get_geometrical_size()) < tol:
                    if np.linalg.norm(np.array(key1.get_sites()) - np.array(key2.get_sites())) < tol:
                        return key2

    raise Exception("Did not find a matching key")


def test_no_symmetry_local_orbitlist_counting(atoms_primitive, atoms_tag, cutoffs, repeatInteger):
    """
    This function creates a primitive orbitlist (no-symmetry) and an orbitlist on supercell
    (symmetry) and then compares the cluster counting over both of them.


    Args:
        atoms_primitive : ASE atom object
        atoms_tag : Human readable string  describing ASE atom object
        cutoffs (list of float): List of cutoff radii, one for each order of cluster
                                 (pairs, triplet, quadruplets, etc)
        repeatInteger (int): Create a number of repeated ASE atom object

    Raises:
        AssertionError : If total count or keys or count of clusters obtained from
                         symmetry and no-symmetry orbitlist are not equal.
        RuntimeError: Testing for rocksalt, hcp, perovskite, non-pbc structures

    """
    atoms = atoms_primitive.copy()
    atoms_repeat = atoms_primitive.copy().repeat(repeatInteger)
    for atom in atoms_repeat:
        if random.random() < 0.5:
            atom.symbol = "H"

    """ Set up neighborlists """
    neighborlists_primitive = get_neighborlists(atoms=atoms, cutoffs=cutoffs)
    structure = structure_from_atoms(atoms)
    structure_repeat = structure_from_atoms(atoms_repeat)
    neighborlists_supercell = get_neighborlists(
        atoms=atoms_repeat, cutoffs=cutoffs)

    """ Orbitlist primitive """
    orbitlist_primitive = orbitList.OrbitList(
        neighborlists_primitive, structure)
    orbitlist_primitive.sort()

    """ Orbitlist supercell """
    orbitlist_supercell = orbitList.OrbitList(
        neighborlists_supercell, structure_repeat)
    orbitlist_supercell.sort()

    """ Set up clustercounts and count cluster """
    clusterCount_local = ClusterCounts()

    clusterCount_local.count_each_local_orbitlist(
        structure_repeat, orbitlist_primitive, False)

    clusterCount_supercell = ClusterCounts()

    clusterCount_supercell.count_orbitlist(
        structure_repeat, orbitlist_supercell, False)

    """ Get the clustercount map """
    cluster_map_local = clusterCount_local.get_cluster_counts()

    cluster_map_supercell = clusterCount_supercell.get_cluster_counts()

    assert len(cluster_map_local) == len(
        cluster_map_supercell), "Testing total count of clusters using no symmetry "\
        "local orbitlist failed for structure {}".format(atoms_tag)

    assert cluster_map_local.keys() == cluster_map_supercell.keys(
        ), "Testing keys of cluster count using no symmetry local orbitlist "\
        "failed for structure {}".format(atoms_tag)
    
    for key1, keyt in zip(sorted(cluster_map_local.keys()), sorted(cluster_map_supercell.keys())):
       key2 = get_equivalent_clustermap_key(key1, cluster_map_supercell.keys())
       assert cluster_map_local[key1] == cluster_map_supercell[key2], "Testing cluster "\
            "count using no symmetry local orbitlist failed for {} when "\
            "counts {} != {}".format(
                atoms_tag, cluster_map_local[key1], cluster_map_supercell[key2])


def get_total_count(cluster_count_dict, cluster):
    """
    Returns the total count of a particular cluster

    Args:
        cluster_count_dict : cluster counting dictionary
        cluster : key of a particular cluster

    """
    count=0
    for counts in cluster_count_dict[cluster]:
        count += cluster_count_dict[cluster][counts]
    return count


def test_no_symmetry_vs_symmetry_count(atoms_primitive, atoms_tag, cutoffs, repeatInteger):
    """
    This function creates a no-symmetry local orbitlist from the neighborlist of a primitive
    cell and compares the cluster counting with the one obtained from symmetry orbitlist created
    from the primitive structure and the cutoff

    Args:
        atoms_primitive : ASE atom object
        atoms_tag : Human readable string  describing ASE atom object
        cutoffs (list of float): List of cutoff radii, one for each order of cluster
                                 (pairs, triplet, quadruplets, etc)
        repeatInteger (int): Create a number of repeated ASE atom object

    Raises:
        AssertionError : If total count or keys or count of clusters obtained from
                         symmetry and no-symmetry orbitlist are not equal.
        RuntimeError: Testing for rocksalt, hcp, perovskite, non-pbc structures

    """
    atoms=atoms_primitive.copy()
    atoms_repeat=atoms_primitive.copy().repeat(repeatInteger)
    for atom in atoms_repeat:
        if random.random() < 0.3:
            atom.symbol="H"

    structure_repeat=structure_from_atoms(atoms_repeat)
    structure=structure_from_atoms(atoms)

    ####################################
    #                                  #
    #      no symmetry counting        #
    #                                  #
    ####################################

    # get neighborlist
    neighborlists=get_neighborlists(atoms=atoms, cutoffs=cutoffs)

    structure=structure_from_atoms(atoms)

    neighborlists=get_neighborlists(atoms=atoms, cutoffs=cutoffs)

    """ Get orbitlist no symmetry case """
    orbitlist_no_symmetry=orbitList.OrbitList(neighborlists, structure)
    orbitlist_no_symmetry.sort()
    """ Set up cluster count and count """
    clusterCount_no_symmetry=ClusterCounts()
    clusterCount_no_symmetry.count_each_local_orbitlist(
        structure_repeat, orbitlist_no_symmetry, False), "Here we use a cutoff so that no extra clusters are found in the symmetry case" \
                                        "and compare the counts found in both methods"

    """ Get the clustercount map """
    cluster_map_no_symmetry=clusterCount_no_symmetry.get_cluster_counts()

    """ Get orbitlist symmetry case """
    orbitlist_symmetry=create_orbit_list(structure, cutoffs, verbosity=0)
    orbitlist_symmetry.sort()

    """ Set up clustercount and count """
    clusterCount_symmetry=ClusterCounts()
    clusterCount_symmetry.count_each_local_orbitlist(
        structure_repeat, orbitlist_symmetry, False)

    """ Get the clustercount map """
    cluster_map_symmetry=clusterCount_symmetry.get_cluster_counts()

    assert orbitlist_symmetry.size() == orbitlist_no_symmetry.size(
    ), "Testing count orbitlist of symmetry case failed for {} when counts "\
    "{}!={}".format(atoms_tag, orbitlist_symmetry.size(),
                    orbitlist_no_symmetry.size())


    for key1, key2 in zip(sorted(cluster_map_no_symmetry.keys()), sorted(cluster_map_symmetry.keys())):
        assert (key1.get_geometrical_size() - key2.get_geometrical_size()
                ) < 1e-3, "sizes for keys are not equal {} != {}".format(key1.get_geometrical_size(), key2.get_geometrical_size())
        assert len(key1.get_distances()) == len(key2.get_distances(
            )), "clusters are not of the same order. {} != {}".format(key1.get_distances(), key2.get_distances())
        assert np.linalg.norm(np.array(key1.get_distances()) - np.array(key2.get_distances())
                              ) < 1e-3, "sizes for keys are not equal {} ! {} ".format(key1.get_distances(), key2.get_distances())

        assert get_total_count(cluster_map_no_symmetry, key1) == get_total_count(
            cluster_map_symmetry, key2), "Testing cluster multiplicity of symmetry case failed for "\
            "{} when total count {}!={} ".format(atoms_tag, get_total_count(cluster_map_symmetry, key1),
                                                     get_total_count(cluster_map_no_symmetry, key1))
        for element_key in cluster_map_no_symmetry[key1]:
            assert cluster_map_no_symmetry[key1][element_key] == cluster_map_symmetry[key2][
                element_key], "Testing element combination in cluster count of symmetry case "\
                "failed for {}. {} != {}".format(
                    atoms_tag, cluster_map_no_symmetry[key1], cluster_map_symmetry[key2])

db=connect("structures_for_testing.db")

for row in db.select():
    atoms_row=row.toatoms()
    atoms_tag=row.tag
    N=3
    cutoffs=[1.4] * 3
    if "NaCl" in atoms_tag:
        cutoffs=[1.1] * 4            
    if "Ni-hcp-hexagonal" in atoms_tag:
        cutoffs = [1.4, 0.5, 0.8]
    if "BaZrO3-perovskite" in atoms_tag or "distorted" in atoms_tag: #these system has wyckoff sites so even with cutoffs at zero the singlets give different size of the orbitlists
        continue
    if atoms_row.get_pbc().all() == True:
        print("Testing structure: {} with cutoffs {}".format(atoms_tag, cutoffs))
        # start = time.time()
        atoms_row.wrap()
        test_no_symmetry_local_orbitlist_counting(
            atoms_row, atoms_tag, cutoffs, N)
        test_no_symmetry_vs_symmetry_count(atoms_row, atoms_tag, cutoffs, N)
        # elapsed_time = time.time() - start
        # print("Test succeed for {} after {}
        # secs".format(atoms_tag,elapsed_time))
