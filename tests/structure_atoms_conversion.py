from ase import Atoms
from ase.build import bulk
import numpy as np
from icetdev import *
from icetdev.structure import *
#Setup bcc WTi
a = bulk('Na').repeat(2)
for at in a:
    if at.index % 2 == 0:
        at.symbol ="Cl"
a.set_pbc((True, True, True))

structure = structure_from_atoms(a)
# Test that each position is equal
assert len(a.positions) == len(structure.get_positions())
for ase_pos, struct_pos in zip(a.positions, structure.get_positions()):
    assert (ase_pos == struct_pos).all()

# Test that each element is equal
assert len(a.get_chemical_symbols()) == len(structure.get_elements())
for ase_symbol, struct_symbol in zip(a.get_chemical_symbols(),
                                     structure.get_elements()):
    assert ase_symbol == struct_symbol

#Test pbc equality
assert (structure.pbc == a.pbc).all()

#Test unit cell equality
assert (structure.cell == a.cell).all()

#assert that structure return as an equal ASE atoms object
conf = structure.to_atoms()
assert conf == a

