# This scripts runs in about 2 seconds on an i7-6700K CPU.

import matplotlib.pyplot as plt
from ase.db import connect
from icet import ClusterExpansion
from pandas import DataFrame

# step 1: Load cluster expansion
ce = ClusterExpansion.read('mixing-energy.ce')

# step 2: Compile predicted and reference data for plotting
data = {'concentration': [], 'reference-energy': [], 'predicted-energy': []}
db = connect('reference-data.db')
for row in db.select('natoms<=6'):
    data['concentration'].append(row.concentration)
    data['reference-energy'].append(row.mixing_energy)
    data['predicted-energy'].append(ce.predict(row.toatoms()))
data = DataFrame.from_dict(data)

# step 3: Plot results
fig, ax = plt.subplots(figsize=(4, 3))
ax.set_xlabel(r'Concentration')
ax.set_ylabel(r'Mixing energy (meV/atom)')
ax.set_xlim([0, 1])
ax.set_ylim([-69, 15])
ax.scatter(data['concentration'], 1e3 * data['reference-energy'],
           marker='o', label='reference')
ax.scatter(data['concentration'], 1e3 * data['predicted-energy'],
           marker='x', label='CE prediction')
plt.savefig('mixing-energy-comparison.png', bbox_inches='tight')
