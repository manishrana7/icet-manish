# This scripts runs in about 1 second on an i7-6700K CPU.

import matplotlib.pyplot as plt
from collections import OrderedDict
from icet import ClusterExpansion
from numpy import array, count_nonzero

# step 1: Collect ECIs in dictionary
ce = ClusterExpansion.read('mixing-energy.ce')
ecis = OrderedDict()
for order in range(len(ce.cluster_space.cutoffs)+2):
    for orbit in ce.cluster_space.orbit_data:
        if orbit['order'] != order:
            continue
        if order not in ecis:
            ecis[order] = {'radius': [], 'parameters': []}
        ecis[order]['radius'].append(orbit['radius'])
        ecis[order]['parameters'].append(ce.parameters[orbit['index']])

# step 2: Plot ECIs
fig, axs = plt.subplots(1, 3, sharey=True, figsize=(7.5, 3))
for k, (order, data) in enumerate(ecis.items()):
    if k < 2 or k > 4:
        continue
    ax = axs[k-2]
    ax.set_ylim((-6, 34))
    ax.set_xlabel(r'Cluster radius (A)')
    if order == 2:
        ax.set_xlim((1.2, 4.2))
        ax.set_ylabel(r'Effective cluster interaction (meV)')
    if order == 3:
        ax.set_xlim((1.5, 3.9))
    if order == 4:
        ax.set_xlim((1.5, 3.9))
        ax.text(0.05, 0.75, 'zerolet: {:.1f} meV'
                .format(1e3*ecis[0]['parameters'][0]),
                transform=ax.transAxes)
        ax.text(0.05, 0.65, 'singlet: {:.1f} meV'
                .format(1e3*ecis[1]['parameters'][0]),
                transform=ax.transAxes)
    ax.plot([0, 5], [0, 0], color='black')
    ax.bar(data['radius'], 1e3*array(data['parameters']), width=0.05)
    ax.scatter(data['radius'], len(data['radius']) * [-5],
               marker='o', s=2.0)
    ax.text(0.05, 0.91, 'order: {} (#pars {}/{})'
            .format(order, len(data['parameters']),
                    count_nonzero(data['parameters'])),
            transform=ax.transAxes,)
plt.savefig('ecis.png', bbox_inches='tight')
