from icetdev import ClusterSpace, StructureContainer
from ase.db import connect
from ase.build import bulk

# step 1: Setting up the basic structure and a cluster space
prim = bulk('Au')
cutoffs = [6.0, 5.0, 4.0]
subelements = ['Ag', 'Au']
cs = ClusterSpace(prim, cutoffs, subelements)
print(cs)

# step 2: Parsing input structures and setting up a structure container
db = connect('structures.db')
atoms_list = []
properties = []
for row in db.select():
    atoms = row.toatoms()
    properties.append({'energy': row.energy})
    atoms.set_positions(row.data['original_positions'])
    atoms_list.append(atoms)
sc = StructureContainer(cs, atoms_list, properties)
print(sc)

# step 3: Training cluster expansions
