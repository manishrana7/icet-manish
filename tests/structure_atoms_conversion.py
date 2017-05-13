from ase import Atoms
from ase.build import bulk
import numpy as np
from icetdev import *
from icetdev.structure import *
#Setup bcc WTi
a = bulk('TiW',"bcc",a=3.321).repeat(2)
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



#True to return structure as ASE atoms object
conf = structure.to_atoms()