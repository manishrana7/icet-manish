import numpy as np
import os

from ase import Atoms
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.io import write

from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.file_IO import write_force_constants_to_hdf5


# Start reference_forces
size = 3
number_of_structures = 16
calc = EMT()

structures = []
for k in range(number_of_structures):
    # prepare configuration with displacements
    atoms = bulk('Al', cubic=True).repeat(size)
    ideal_positions = atoms.get_positions()
    if k == 0:
        # The first configuration in the training set mimics exactly the
        # displacement pattern requested by phonopy (see below).
        atoms[0].position += [0.01, 0, 0]
    elif k <= 5:
        atoms.rattle(0.05, seed=np.random.randint(1, 1e8))
    elif k <= 10:
        atoms.rattle(0.10, seed=np.random.randint(1, 1e8))
    else:
        atoms.rattle(0.15, seed=np.random.randint(1, 1e8))
    displacements = atoms.positions - ideal_positions

    # compute forces
    atoms.set_calculator(calc)
    forces = atoms.get_forces()

    # attach forces and displacements to Atoms object
    atoms.new_array('forces', forces)
    atoms.new_array('displacements', displacements)
    atoms.positions = ideal_positions
    structures.append(atoms.copy())

# write reference data to file
try:
    os.makedirs('reference-data')
except:
    pass
write('reference-data/training-configurations.xyz', structures)
# End reference_forces

# Start reference_force_constants
# prepare primitive unit cell
unitcell_ASE = bulk('Al')
unitcell = PhonopyAtoms(symbols=unitcell_ASE.get_chemical_symbols(),
                        cell=unitcell_ASE.cell,
                        scaled_positions=unitcell_ASE.get_scaled_positions())

# prepare supercell
phonon = Phonopy(unitcell,
                 3 * np.array([[-1, 1, 1],
                               [1, -1, 1],
                               [1, 1, -1]]))
phonon.generate_displacements(distance=0.01)
supercells = phonon.get_supercells_with_displacements()

# compute force constant matrix
forces = []
for atoms in supercells:
    atoms_ASE = Atoms(symbols=atoms.get_chemical_symbols(),
                      cell=atoms.get_cell(),
                      scaled_positions=atoms.get_scaled_positions(),
                      pbc=True)
    atoms_ASE.set_calculator(EMT())
    forces.append(atoms_ASE.get_forces().copy())
phonon.set_forces(forces)
phonon.produce_force_constants()

# write force constant matrix and ideal supercell to file
fname = 'reference-data/fc2_phonopy.hdf5'.format(len(atoms_ASE))
write_force_constants_to_hdf5(phonon.force_constants, filename=fname)
fname = 'reference-data/SPOSCAR.phonopy'.format(len(atoms_ASE))
atoms = phonon.get_supercell()
atoms_ASE = Atoms(symbols=atoms.get_chemical_symbols(),
                  cell=atoms.get_cell(),
                  scaled_positions=atoms.get_scaled_positions(),
                  pbc=True)
write(fname, atoms_ASE)
# End reference_force_constants
