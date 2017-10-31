#Step 1
from icetdev.clusterspace import create_clusterspace
from ase.build import bulk
from ase import Atoms


conf = bulk("Si")
cutoffs = [5.0, 5.0, 5.0]
subelements = ["Si", "Ge"]

clusterspace = create_clusterspace(conf, cutoffs, subelements)
#Step 1.1
print(clusterspace)

#Step 2

supercell = bulk("Si").repeat(2)

cv = clusterspace.get_clustervector(supercell)
print(cv)