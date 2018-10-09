from ase.build import make_supercell
from numpy import arange, array
from time import time

from icet import ClusterExpansion
from mchammer.calculators import ClusterExpansionCalculator
from mchammer.ensembles import SemiGrandCanonicalEnsemble

# step 1: set up the structure to simulate
ce = ClusterExpansion.read('mixing_energy.ce')
chemical_symbols = ce.cluster_space.chemical_symbols
atoms = ce.cluster_space.primitive_structure.repeat(3)
atoms = make_supercell(ce.cluster_space.primitive_structure,
                       3*array([[-1, 1, 1],
                                [1, -1, 1],
                                [1, 1, -1]]))

# step 2: set up calculator
print('setting up calculator')
start = time()
calculator = ClusterExpansionCalculator(atoms, ce)
print(' done in {:.2f} sec'.format(time() - start))

# step 3: set up ensemble
# TODO: remove temperature and chemical_potentials once possible
mc = SemiGrandCanonicalEnsemble(
    calculator=calculator,
    atoms=atoms,
    ensemble_data_write_interval=10,
    temperature=1,
    chemical_potentials={chemical_symbols[0]: 0,
                         chemical_symbols[1]: 0})

# step 4: execute Monte Carlo runs
n_production = 1000
n_equilibration = 400
for temperature in [600, 900, 300]:
    for dmu in arange(-0.6, 0.53, 0.04):
        print('temperature: {}   dmu: {:.3f}'.format(temperature, dmu))
        start = time()
        mc.temperature = temperature
        mc.chemical_potentials = {chemical_symbols[0]: 0,
                                  chemical_symbols[1]: dmu}
        mc.run(n_equilibration)
        mc.reset_data_container()
        mc.run(n_production)
        # TODO: change the next line once mc.data_container is writable
        mc.data_container.write('sgc-T{}-dmu{:.3f}.dc'
                                .format(temperature, dmu))
        print(' done in {:.2f} sec'.format(time() - start))
