.. _example_enumerate_structures:
.. highlight:: python
.. index::
   single: Examples; Structure enumeration

Structure enumeration
=====================

The purpose of this example is to demonstrate how to enumerate structures. In
the present context, this structural enumeration means the generation of all
inequivalent structures derived from a primitive structure up to a certain
size.

Import modules
--------------

The :func:`enumerate_structures
<icet.tools.structure_enumeration.enumerate_structures>` function needs to
be imported together with some additional functions from `ASE
<https://wiki.fysik.dtu.dk/ase>`_.

.. literalinclude:: ../../../../examples/enumerate_structures.py
   :start-after: # Import modules
   :end-before: # Generate all binary

Generate binary Au/Pd structures
--------------------------------

Before being able to perform the structural enumeration, it is first necessary
to generate a primitive structure. In this case, an Au fcc :class:`ASE Atoms`
object is created using the :func:`ase.build.bulk` function. Then a database
``AuPd-fcc.db`` is initialized, in which the enumerated structures will be
stored. All possible binary Au/Pd structures with up to 6 atoms per unit cell
are subsequently generated and stored in this database.

.. literalinclude:: ../../../../examples/enumerate_structures.py
   :start-after: # and save them
   :end-before: # Enumerate all palladium

Generate PdH structures with vacancies
--------------------------------------

The steps above are now repeated to enumerate all palladium hydride structures
based on up to four primitive cells, which contain up to 4 Pd atoms and between
0 and 4 H atoms. Vacancies, represented by vanadium, are included, which
results in a ternary system. The structures thus obtained are stored in a
database named ``PdHVac-fcc.db``.

.. literalinclude:: ../../../../examples/enumerate_structures.py
   :start-after: # either a hydrogen
   :end-before: # Enumerate a copper surface

Generate surface slabs with adsorbates
--------------------------------------

Lower dimensional systems can also be enumerated. Here, this is demonstrated
with a copper surface with oxygen atoms adsorbed in hollow sites on a {111}
surface. The key to trigger a two- or one-dimensional enumeration is to make
sure that the periodic boundary conditions of the input structure reflect the
desired behavior. For the surface system, this means that the the boundary
conditions are *not* periodic in the direction of the normal to the surface.
This is the default behavior with ASE:s surface building functions, but is in
the below example enforced for clarity.

.. literalinclude:: ../../../../examples/enumerate_structures.py
   :start-after: # fcc and hcp hollow sites

Source code
-----------

.. container:: toggle

    .. container:: header

       The complete source code is available in
       ``examples/enumerate_structures.py``

    .. literalinclude:: ../../../../examples/enumerate_structures.py
