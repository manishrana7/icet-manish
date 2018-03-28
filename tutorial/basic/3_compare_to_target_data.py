import matplotlib.pyplot as plt
import numpy as np
from ase.db import connect
from icet import ClusterExpansion

# step 1: Compare predicted and target data
ce = ClusterExpansion.read('cluster_expansion.icet')
data = []
db = connect('structures.db')
for row in db.select():
    emix_ce = ce.predict(row.toatoms())
    data.append([row.conc, row.emix, emix_ce])
data = np.array(data).T

# step 2: Plot results
fig, ax = plt.subplots(figsize=(4, 3))
ax.set_xlabel(r'Ag concentration')
ax.set_ylabel(r'Mixing energy (meV/atom)')
ax.set_xlim([0, 1])
ax.scatter(data[0], 1e3 * data[1], marker='o')
ax.scatter(data[0], 1e3 * data[2], marker='x')
plt.savefig('mixing-energy-comparison.png', bbox_inches='tight')
