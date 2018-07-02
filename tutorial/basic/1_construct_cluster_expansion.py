# This scripts runs in about 2 seconds on an i7-6700K CPU.

from ase.db import connect
from ase.build import bulk
from icet import (ClusterSpace,
                  StructureContainer,
                  Optimizer,
                  ClusterExpansion)

# step 0: Setup
subelements = ['Ag', 'Pd']
prim = bulk(subelements[0])

# step 1: Set up the basic structure and a cluster space
cutoffs = [8.0, 7.0, 6.0]
cs = ClusterSpace(atoms=prim, cutoffs=cutoffs, chemical_symbols=subelements)
print(cs)

# step 2: Parse the input structures and set up a structure container
db = connect('reference-data.db')
sc = StructureContainer(cluster_space=cs)
for row in db.select('natoms<=6'):
    sc.add_structure(atoms=row.toatoms(),
                     properties={'mixing-energy': row.mixing_energy})
print(sc)

# step 3: Train parameters
opt = Optimizer(fit_data=sc.get_fit_data(key='mixing-energy'),
                fit_method='lasso')
opt.train()
print(opt)

# step 4: Construct cluster expansion and write it to file
ce = ClusterExpansion(cluster_space=cs, parameters=opt.parameters)
ce.write('mixing-energy.ce')
