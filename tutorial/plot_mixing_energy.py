import matplotlib.pyplot as plt
from ase.db import connect
import numpy as np

db = connect('structures.db')
data = []
enAu = 0.0
enAg = 0.0
for row in db.select():
    nAg = row.count_atoms().get('Ag', 0)
    nAu = row.count_atoms().get('Au', 0)
    if row.natoms == 1:

        enAg = row.energy
    if row.natoms == 1 and nAu == 1:
        enAu = row.energy
    conc = float(nAg) / row.natoms
    emix = row.energy / row.natoms - conc * enAg - (1.0 - conc) * enAu
    data.append([conc, emix])
data = np.array(data).T

fig, ax = plt.subplots()
ax.set_xlabel(r'Ag concentration')
ax.set_ylabel(r'Mixing energy (meV/atom)')
ax.set_xlim([0, 1])
ax.scatter(data[0], 1e3 * data[1], marker='x')
plt.savefig('mixing-energy.pdf', bbox_inches='tight')
