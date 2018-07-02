# This scripts runs in about 4.5 minutes on an i7-6700K CPU.

import matplotlib.pyplot as plt
from numpy import array
from icet import ClusterExpansion
from icet.tools import ConvexHull, enumerate_structures

# step 1: Predict energies for enumerated structures
ce = ClusterExpansion.read('mixing-energy.ce')
subelements = ce.cluster_space.chemical_symbols
data = {'concentration': [], 'mixing-energy': []}
structures = []
for atoms in enumerate_structures(
        atoms=ce.cluster_space.primitive_structure,
        sizes=range(1, 13),
        chemical_symbols=ce.cluster_space.chemical_symbols):
    conc = atoms.get_chemical_symbols().count(subelements[1]) / len(atoms)
    data['concentration'].append(conc)
    data['mixing-energy'].append(ce.predict(atoms))
    structures.append(atoms)
print('Predicted energies for {} structures'.format(len(structures)))

# step 2: Construct convex hull
hull = ConvexHull(data['concentration'], data['mixing-energy'])

# step 3: Plot the results
fig, ax = plt.subplots(figsize=(4, 3))
ax.set_xlabel(r'Concentration')
ax.set_ylabel(r'Mixing energy (meV/atom)')
ax.set_xlim([0, 1])
ax.set_ylim([-69, 15])
ax.scatter(data['concentration'], 1e3 * array(data['mixing-energy']),
           marker='x')
ax.plot(hull.concentrations, 1e3 * hull.energies, '-o', color='green')
plt.savefig('mixing-energy-predicted.png', bbox_inches='tight')

# step 4: Extract candidate ground state structures
tol = 0.0005
low_energy_structures = hull.extract_low_energy_structures(
    data['concentration'], data['mixing-energy'], tol, structures)
print('Found {} structures within {} meV/atom of the convex hull'.
      format(len(low_energy_structures), 1e3*tol))
