import matplotlib.pyplot as plt
import numpy as np
from icet import ClusterExpansion
from collections import OrderedDict

# step 1: Collect ECIs in dictionary
ce = ClusterExpansion.read('cluster_expansion.icet')
ecis = OrderedDict()
for order in range(len(ce.cluster_space.cutoffs)+2):
    for orbit in ce.cluster_space.orbit_data:
        if orbit['order'] != order:
            continue
        if order not in ecis:
            ecis[order] = []
        ecis[order].append([orbit['size'], ce.parameters[orbit['index']]])

# step 2: Plot ECIs
fig, axs = plt.subplots(1, 2, sharey=True, figsize=(5, 3))
for k, (order, data) in enumerate(ecis.items()):
    if k < 2 or k >= 4:
        continue
    ax = axs[k-2]
    ax.set_xlim((1.2, 3.2))
    ax.set_ylim((-0.5, 7.5))
    ax.set_xlabel(r'Cluster radius (A)')
    if order == 2:
        ax.set_ylabel(r'Effective cluster interaction (meV)')
        ax.text(1.65, 4.8, 'zerolet: {:.1f} meV'.format(1e3*ecis[0][0][1]))
        ax.text(1.65, 3.8, 'singlet: {:.1f} meV'.format(1e3*ecis[1][0][1]))
    data = np.transpose(data)
    ax.plot([0, 5], [0, 0], color='black')
    ax.bar(data[0], 1e3*data[1], width=0.1)
    ax.text(0.05, 0.91, 'order: {}'.format(order),
            transform=ax.transAxes,)
plt.savefig('ecis.png', bbox_inches='tight')
