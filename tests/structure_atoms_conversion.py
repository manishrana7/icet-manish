from ase.db import connect
import icetdev as icet

'''
Testing icet structure against ASE
'''

''' Fetch structures from database '''
db = connect('structures_for_testing.db')

for row in db.select():

    atoms_row = row.toatoms()

    ''' Convert ASE atoms to icet structures '''
    structure = icet.structure_from_atoms(atoms_row)

    ''' Test that each position is equal '''
    msg = 'Test for len of positions failed for structure {}'.format(row.tag)
    assert len(atoms_row) == len(structure.get_positions()), msg

    for ase_pos, struct_pos in zip(atoms_row.positions,
                                   structure.get_positions()):
        msg = 'Test for positions failed for structure {}'.format(row.tag)
        assert (ase_pos == struct_pos).all(), msg

    ''' Test that each element is equal '''
    msg = 'Test for len of elements failed for structure {}'.format(row.tag)
    assert len(atoms_row) == len(structure.get_elements()), msg

    for ase_symbol, struct_symbol in zip(atoms_row.get_atomic_numbers(),
                                         structure.get_elements()):
        msg = 'Test for elements failed for structure {}'.format(row.tag)
        msg += '; {} != {}'.format(ase_symbol, struct_symbol)
        assert ase_symbol == struct_symbol, msg

    ''' Test periodic boundary conditions '''
    msg = 'Test for periodic boundary conditions failed'
    msg += ' for structure {}'.format(row.tag)
    assert (structure.pbc == atoms_row.pbc).all(), msg

    ''' Test unit cell equality '''
    msg = 'Test for cell failed for structure {}'.format(row.tag)
    assert (structure.cell == atoms_row.cell).all(), msg

    ''' Assert that structure return as an equal ASE atoms object '''
    atoms_icet = structure.to_atoms()
    msg = 'Test for conversion back to atoms failed'
    msg += ' for structure {}'.format(row.tag)
    assert atoms_icet == atoms_row, msg
