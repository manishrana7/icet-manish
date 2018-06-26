from icet import ClusterExpansion
from mchammer.calculators import ClusterExpansionCalculator
from mchammer.ensembles import SemiGrandCanonicalEnsemble
from numpy import arange, array

# step 1: set up the structure to simulate
ce = ClusterExpansion.read('cluster_expansion.icet')
atoms = ce.cluster_space.primitive_structure.repeat(3)
calculator = ClusterExpansionCalculator(atoms, ce)

# step 2: set up calculator and ensemble
ntrials = 500
equil = 200
for temperature in [1200, 600, 900, 300]:
    for dmu in arange(-1.78, 1.91, 0.02):
        mc = SemiGrandCanonicalEnsemble(
            calculator=calculator, atoms=atoms,
            data_container='sgc.dc',
            random_seed=42, temperature=temperature,
            chemical_potentials={'Au': 0, 'Ag': dmu},
            ensemble_data_write_interval=2)
        mc.run(ntrials)
