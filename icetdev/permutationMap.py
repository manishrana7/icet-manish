import numpy as np
import spglib as spglib

from ase import Atoms
from ase.visualize import view
from _icetdev import PermutationMap
from icetdev.neighborlist import *
from icetdev.latticeNeighbor import LatticeNeighbor
from icetdev.structure import structure_from_atoms
from icetdev.tools.geometry import get_scaled_positions


def __get_primitive_structure(atoms):
    """
    Returns primitive atoms object.
    """
    lattice, scaled_positions, numbers = spglib.standardize_cell(
        atoms, to_primitive=True)

    # create the primitive atoms object
    atoms_prim = Atoms(scaled_positions=scaled_positions,
                       numbers=numbers, cell=lattice, pbc=atoms.pbc)
    return atoms_prim


def __get_neigbhorlists(atoms, cutoffs):
    """
    Get a list of neigbhorlist objects (one for each cutoff)
    Also returns the created structure object
    """
    structure = structure_from_atoms(atoms)
    neighborlists = []
    for co in cutoffs:
        nl = Neighborlist(co)
        nl.build(structure)
        neighborlists.append(nl)
    return structure, neighborlists


def __get_fractional_positions_from_nl(structure, neighborlist):
    """
    Returns the fractional positions in structure from the neighbors in the nl

    """
    position_of_neighbors = []
    fractional_positions = []
    latnbr_i = LatticeNeighbor(0, [0, 0, 0])
    for i in range(neighborlist.size()):
        latnbr_i.index = i
        position = structure.get_position(latnbr_i)
        position_of_neighbors.append(position)
        for latNbr in neighborlist.get_neighbors(i):
            position = structure.get_position(latNbr)
            position_of_neighbors.append(position)
    if len(position_of_neighbors) > 0:
        fractional_positions = get_scaled_positions(np.array(position_of_neighbors),
                                                    structure.cell, wrap=False,

                                                    pbc=structure.pbc)
    return fractional_positions


def permutation_maps_from_atoms(atoms, cutoffs=None, find_prim=True, verbosity=0):
    """
    Setup a list of permutation maps from an atoms object.

    Keyword arguments:
        cutoffs -- list of cutoffs for the clusters wanted.
        find_primitive -- if true it will take symmetries from the primitive
        atoms object (default True)
        verbosity -- 0 or lower is no verbosity, 1 is some output 3 is debug mode (default 0)
    """
    atoms = atoms.copy()
    # set each element to the same since we only care about geometry when
    # taking primitive
    atoms.set_chemical_symbols(len(atoms) * [atoms[0].symbol])

    atoms_prim = atoms
    if find_prim:
        atoms_prim = __get_primitive_structure(atoms)

    if verbosity >= 3:
        print("size of atoms_prim {}".format(len(atoms_prim)))
    # Get symmetry information and load into a permutation map object
    symmetry = spglib.get_symmetry(atoms_prim)
    translations = symmetry['translations']
    rotations = symmetry['rotations']
    permutation_maps = [PermutationMap(
        translations, rotations) for i in range(len(cutoffs))]

    # Create neighborlists from the different cutoffs
    prim_structure, neighborlists = __get_neigbhorlists(atoms_prim, cutoffs)

    # get fractional positions for each neighborlist
    for i, neighborlist in enumerate(neighborlists):
        if verbosity >= 3:
            print("building permutation map {}/{}".format(i, len(neighborlists)))
        frac_positions = __get_fractional_positions_from_nl(
            prim_structure, neighborlist)
        if verbosity >= 3:
            print("number of fractional positions: {} ".format(len(frac_positions)))
        if len(frac_positions) > 0:
            permutation_maps[i].build(frac_positions)

    return permutation_maps, prim_structure, neighborlists
