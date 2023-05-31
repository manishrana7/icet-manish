"""
This examples demonstrates how one can generate trainings structures
"""

# Import modules
import numpy as np
from ase.build import bulk
from icet import ClusterSpace
from icet.tools.structure_generation import occupy_structure_randomly
from icet.tools.training_set_generation import structure_selection_annealing

# For plotting
import matplotlib.pyplot as plt


# Convenience function for supercell size generation
def get_random_supercell_size(max_repeat, max_atoms, n_atoms_in_prim):
    while True:
        nx, ny, nz = np.random.randint(1, max_repeat + 1, size=3)
        if nx * ny * nz * n_atoms_in_prim < max_atoms:
            break
    return nx, ny, nz


# Create the primitive structure and cluster space.
# The possible occupations are Au and Pd
primitive_structure = bulk('Au', 'fcc', 4.0)
cutoffs = [10.0, 6.0, 4.0]
subelements = ['Au', 'Pd']
cluster_space = ClusterSpace(primitive_structure, cutoffs, subelements)

# Create a random structure pool
n_random_structures = 10000
max_repeat = 8
max_atoms = 50

structures = []
for _ in range(n_random_structures):
    # Create random supercell.
    supercell = get_random_supercell_size(max_repeat, max_atoms, len(primitive_structure))
    structure = primitive_structure.repeat(supercell)

    # Randomize concentrations in the supercell
    n_atoms = len(structure)
    n_Au = np.random.randint(0, n_atoms)
    n_Pd = n_atoms - n_Au
    concentration = {'Au': n_Au / n_atoms, 'Pd': n_Pd / n_atoms}

    # Occupy the structure randomly and store it.
    occupy_structure_randomly(structure, cluster_space, concentration)
    structures.append(structure)

# We want to add 2 times the number of parameters structures and we want the annealing to run
# for 1e4 steps.
n_structures_to_add = 2 * len(cluster_space)
n_steps = 10000

# start the annealing procedure to minimize the condition number of
# the fit matrix.
indices, traj = structure_selection_annealing(cluster_space, structures, n_structures_to_add,
                                              n_steps)
condition_number_annealing = traj[-1]

# Since we start with random structures from the pool
# this represents choosing structures on random
condition_number_random_structures = traj[0]

# Collect the structures that were found to be good training structure candidates.
training_structures = [structures[ind] for ind in indices]

# Plot the metric vs accepted trials.
fig, ax = plt.subplots()
ax.plot(traj)
ax.set_xlabel('Accepted trials')
ax.set_ylabel('Condition number')
fig.tight_layout()
fig.savefig('training_set_generation_cond_traj.svg')

print('='*20 + ' Condition number for annealing and random structures ' + '='*20 + '\n'
      f'annealing structures: {condition_number_annealing}\n'
      f'random structures: {condition_number_random_structures}')
