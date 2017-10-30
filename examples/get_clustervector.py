from icetdev.clusterspace import create_clusterspace
from ase.build import bulk
from ase import Atoms

conf = bulk("Si")
cutoffs = [5.0, 5.0, 5.0]
subelements = ["Si", "Ge"]
print(subelements)
print(type(conf))
exit(1)
clusterspace = create_clusterspace(conf, cutoffs, subelements)
