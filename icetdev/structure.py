from example import Structure



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
