from ase import Atoms
from ase.build import bulk, fcc111
from ase.cluster.icosahedron import Icosahedron
from ase.db import connect

"""
Database containing different type of structures used for testing.
Some tests may use partially the following database.

In case you re-edit the id for the structures please be sure to delete any exisiting 
database otherwise the structure will be duplicated.
"""

structures = []

""" FCC (single element, primitive cell, pbc=True) """
atoms=bulk('Al','fcc', a=1.0)
atoms.tag = 'Al-fcc-prim'
structures.append(atoms)

""" FCC (single element, supercell, pbc=True) """
atoms = bulk('Al', 'fcc', a=1.0).repeat(2)
atoms.tag = 'Al-fcc-supercell'
structures.append(atoms)

""" FCC (single element, distorted structure, pbc=True) """
atoms = bulk('Al', 'fcc', a=1.0).repeat(2)
atoms.rattle(stdev=0.001, seed=42)
atoms.tag = 'Al-fcc-distorted'
structures.append(atoms)

""" BCC (two elements, cubic structure, pbc=True) """
atoms = bulk('Ti','bcc', a=1.0, cubic=True).repeat(2)
for atom in atoms:
    if atom.index % 2 == 0:
        atom.symbol='W'
atoms.tag = 'TiW-bcc-supercell'
structures.append(atoms)

""" rocksalt (two elements, cubic structure) """
atoms = bulk('NaCl', 'rocksalt', a=1.0)
atoms.tag = 'NaCl-rocksalt-prim'
structures.append(atoms)

""" HCP (single element, hexagonal structure) """
atoms = bulk('Ni', 'hcp', a=0.625, c=1.0)
atoms.tag = 'Ni-hexagonal-prim'
structures.append(atoms)

""" perovskite (three elements, cubic structure) """
a=1.0
atoms = Atoms('BaZrO3', positions=[(0, 0, 0), (a/2,a/2,a/2), 
                (a/2,a/2,0), (a/2,0,a/2), (0,a/2,a/2)], cell=[a, a, a], pbc=True)
atoms.tag = 'BaZrO3-perovskite'
structures.append(atoms)


""" surface slab (two elements, pbc=[True, True, False]) """
atoms = fcc111('Pd', a=1.0, size=(2,2,3))
atoms.center(vacuum=4.0, axis=2)
atoms.tag = 'Pd-slab-surface'
structures.append(atoms)

""" Nanoparticle (single element, pbc=False) """
atoms = Icosahedron('Au', noshells=3, latticeconstant=1.0)
atoms.center(vacuum=6.0, axis=(0, 1, 2))
atoms.tag = 'Au-nanoparticle'
structures.append(atoms)



with connect('structures_for_testing.db') as db:
    tags = [row.tag for row in db.select()]
    for atoms in structures:
        if atoms.tag in tags:
            continue
        else: 
            db.write(atoms, tag=atoms.tag)
 