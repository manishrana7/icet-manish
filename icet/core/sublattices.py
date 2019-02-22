
from copy import deepcopy
from icet.core.structure import Structure
from typing import List
from ase import Atoms


class Sublattices:
    """
    This class stores and provides information about the sublattices
    of a structure


    Parameters
    ----------
    allowed_species
        list of the allowed species on each site of the primitve
        structure. For example this can be the chemical_symols from
        a cluster space
    primitive_structure
        the primitive structure the allowed species reference to
    structure
        the structure that the sublattices will be based on
    """

    def __init__(self, allowed_species: List[List[str]],
                 primitive_structure: Atoms, structure: Atoms):

        unique_sites = list(set(tuple(sorted(symbols)) for symbols in allowed_species))

        # sorted unique sites, this basically decides A, B, C... sublattices
        self._allowed_species = sorted(unique_sites)
        n_sublattices = len(unique_sites)
        self._species_to_sublattice = {}
        self._index_to_sublattice = {}

        for i, symbols in enumerate(self._allowed_species):
            for symbol in symbols:
                self._species_to_sublattice[symbol] = i

        cpp_prim_structure = Structure.from_atoms(primitive_structure)

        self._sublattice_to_indices = [[] for _ in range(n_sublattices)]
        for index, position in enumerate(structure.get_positions()):
            lattice_site = cpp_prim_structure.find_lattice_site_by_position(
                position)
            species = allowed_species[lattice_site.index]

            sublattice = self._allowed_species.index(tuple(sorted(species)))
            self._index_to_sublattice[index] = sublattice
            self._sublattice_to_indices[sublattice].append(index)

    def get_sublattice_index(self, symbol: str = None,
                             index: int = None) -> int:
        """ Returns the index of the sublattice the symbol
            or index in the structure belongs to

            Parameters
            -----------
            symbol
                species symbol
            index
                index of site in the structure
        """

        if symbol is not None:
            return self._species_to_sublattice[symbol]
        elif index is not None:
            return self._index_to_sublattice[index]
        else:
            raise ValueError("either symbol or index must be supplied")

    @property
    def allowed_species(self)->List[List[str]]:
        """Lists of the allowed species on each sublattice, in order"""
        return deepcopy(self._allowed_species)

    def get_sublattice_sites(self, index) ->List[int]:
        """Returns the sites that belong to the sublattice with the
            corresponding index

            Parameters
            ----------
            index
                index of the sublattice
        """
        return self._sublattice_to_indices[index]
