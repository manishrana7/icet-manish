from _icetdev import Structure


def structure_from_atoms(conf):
    """ 
    Returns an icet structure object from ASE atoms object.

    Parameters
    ----------
    conf:
        ASE atoms object.

    Returns
    -------
    icet structure object
    """

    return Structure(conf.positions,
                     conf.get_chemical_symbols(),
                     conf.cell,
                     conf.pbc.tolist())


def __Structure_to_atoms(self):
    """
    Returns the structure as an ASE atoms object.


    Returns
    -------
    ASE atoms object
        atomic structure
    """
    import ase
    conf = ase.Atoms(pbc=self.pbc)
    conf.set_cell(self.cell)
    for element, position in zip(self.elements, self.positions):
        conf.append(ase.Atom(element, position))
    conf.set_positions(self.get_positions())
    conf.set_chemical_symbols(self.get_elements())
    return conf


Structure.to_atoms = __Structure_to_atoms



def __structure_pbc(self):
    """
    list of three booleans: Periodic boundary conditions.
    """
    return self.get_pbc()


Structure.pbc = property(__structure_pbc)


def __structure_cell(self):
    """
    numpy matrix : unit cell
    """
    return self.get_cell()


Structure.cell = property(__structure_cell)


def __structure_positions(self):
    """
    N x 3 numpy matrix : positions
    """
    return self.get_positions()


Structure.positions = property(__structure_positions)


def __structure_elements(self):
    """
    Numpy array  of type str(check this) : elements
    """
    return self.get_elements()


Structure.elements = property(__structure_elements)


def __get_size(self):
    """
    Get size of structure (number of elements)
    """
    return len(self.elements)


Structure.size = __get_size


def __find_index_from_position(self, position, tolerance=1e-6):

    index = self.find_index_of_position_pybind(position, tolerance)
    if index == -1:
        raise ValueError(
            "Did not find index of position in function find_latticeNeighbor_from_position")
    else:
        return index


Structure.find_index_of_position = __find_index_from_position


def __get_len(self):
    """
    returns size of structure
    """
    return self.size()


Structure.__len__ = __get_len


def __str_function(self):
    """
    Define the str function for structure
    """

    str = "Cell: \n"
    str += """ {} \n \n""".format(self.cell)
    str += """Element and positions: \n"""
    for pos, element in zip(self.positions, self.elements):
        str += " {}  {} \n".format(element, pos)

    return str
Structure.__str__ = __str_function
