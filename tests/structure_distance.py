from example import Structure
from ase import Atoms
# Setup a chain of H,O,C
# H-O Dist = 2
# O-C Dist = 3
# C-H Dist = 5 with mic=False
# C-H Dist = 4 with mic=True
a = Atoms('HOC', positions=[(1, 1, 1), (3, 1, 1), (6, 1, 1)])
a.set_cell((9, 2, 2))
a.set_pbc((True, False, False))

ex_structure = Structure(a.positions, a.get_chemical_symbols())