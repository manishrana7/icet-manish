from icet import ClusterExpansion
from mchammer.calculators import ClusterExpansionCalculator
from mchammer.ensembles import SemiGrandCanonicalEnsemble
from numpy import arange

# step 1: set up the structure to simulate
ce = ClusterExpansion.read('mixing-energy.ce')
chemical_symbols = ce.cluster_space.chemical_symbols
atoms = ce.cluster_space.primitive_structure.repeat(3)
calculator = ClusterExpansionCalculator(atoms, ce)

# step 2: set up calculator and ensemble
ntrials = 500
equil = 200
for temperature in [600, 900, 300]:
    for dmu in arange(-2.0, 2.0, 0.1):
        mc = SemiGrandCanonicalEnsemble(
            calculator=calculator, atoms=atoms,
            data_container='sgc.dc',
            random_seed=42, temperature=temperature,
            chemical_potentials={chemical_symbols[0]: 0,
                                 chemical_symbols[1]: dmu},
            ensemble_data_write_interval=10)
        mc.run(ntrials)
