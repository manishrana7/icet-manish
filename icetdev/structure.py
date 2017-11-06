from _icetdev import Structure


def structure_from_atoms(conf):
    '''
    Creates an icet structure object from an ASE atoms object.

    Parameters
    ----------
    conf : ASE atoms object
        input configuration

    Returns
    -------
    icet structure object
        output configuration
    '''
    return Structure(conf.positions,
                     conf.get_chemical_symbols(),
                     conf.cell,
                     conf.pbc.tolist())


def __Structure_to_atoms(self):
    '''
    Returns the structure as an ASE atoms object.

    Returns
    -------
    ASE atoms object
        atomic configuration
    '''
    import ase
    conf = ase.Atoms(pbc=self.pbc)
    conf.set_cell(self.cell)
    for symbol, position in zip(self.chemical_symbols, self.positions):
        conf.append(ase.Atom(symbol, position))
    conf.set_positions(self.get_positions())
    conf.set_chemical_symbols(self.get_chemical_symbols())
    return conf


Structure.to_atoms = __Structure_to_atoms


def __repr_function(self):
    s = ['Cell:']
    s += ['{}\n'.format(self.cell)]
    s += ['Element and positions:']
    for symbol, position in zip(self.chemical_symbols, self.positions):
        s += [' {}  {}'.format(symbol, position)]
    return '\n'.join(s)


Structure.__repr__ = __repr_function
