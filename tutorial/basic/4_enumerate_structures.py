import matplotlib.pyplot as plt
import numpy as np
from icet import ClusterExpansion
from icet.tools import ConvexHull, enumerate_structures

# step 1: Predict energies for enumerated structures
ce = ClusterExpansion.read('cluster_expansion.icet')
data = []
structures = []
for atoms in enumerate_structures(ce.cluster_space.primitive_structure,
                                  range(1, 9),
                                  ce.cluster_space.chemical_symbols):
    conc = float(atoms.get_chemical_symbols().count('Ag')) / len(atoms)
    emix = ce.predict(atoms)
    data.append([conc, emix])
    structures.append(atoms)
print('Predicted energies for {} structures'.format(len(data)))
data = np.array(data).T

# step 2: Construct convex hull
hull = ConvexHull(data[0], data[1])

# step 3: Plot the results
fig, ax = plt.subplots(figsize=(4, 3))
ax.set_xlabel(r'Ag concentration')
ax.set_ylabel(r'Mixing energy (meV/atom)')
ax.set_xlim([0, 1])
ax.scatter(data[0], 1e3 * data[1], marker='x')
ax.plot(hull.concentrations, 1e3 * hull.energies, '-o', color='green')
plt.savefig('mixing-energy-predicted.png', bbox_inches='tight')

# step 4: Extract candidate ground state structures
tol = 0.0005
low_energy_structures = hull.extract_low_energy_structures(data[0], data[1],
                                                           tol, structures)
print('Found {} structures within {} meV/atom of the convex hull'.
      format(len(low_energy_structures), 1e3*tol))
