from ase.db import connect
from ase.build import bulk
from icet import (ClusterSpace,
                  StructureContainer,
                  Optimizer,
                  ClusterExpansion)

# step 1: Set up the basic structure and a cluster space
prim = bulk('Au')
cutoffs = [6.0, 5.0, 4.0]
subelements = ['Ag', 'Au']
cs = ClusterSpace(prim, cutoffs, subelements)
print(cs)

# step 2: Parse the input structures and set up a structure container
db = connect('structures.db')
sc = StructureContainer(cs)
for row in db.select():
    sc.add_structure(row.toatoms(), user_tag=str(row.structure_id),
                     properties={'energy': row.emix})
print(sc)

# step 3: Train parameters
opt = Optimizer(sc.get_fit_data())
opt.train()
print(opt)

# step 4: Construct cluster expansion and write it to file
ce = ClusterExpansion(cs, opt.parameters)
ce.write('cluster_expansion.icet')
