from ase.build import bulk
from ase.db import connect
from ase.calculators.emt import EMT
from ase.optimize import BFGS
from icetdev.enumeration import enumerate_structures

# step 1: Prepare database and set up basic structure
db = connect('structures.db')
prim = bulk('Au')
subelements = ['Ag', 'Au']

# step 2: Enumerate structures, then relax and add them to the database
sizes = range(1, 9)
for atoms in enumerate_structures(prim, sizes, subelements):
    # remember the original (unrelaxed) positions
    original_positions = atoms.get_positions()
    # relax structure
    atoms.calc = EMT()
    dyn = BFGS(atoms)
    dyn.run(fmax=0.01)
    # add structure to database
    db.write(atoms, data={'original_positions': original_positions})
