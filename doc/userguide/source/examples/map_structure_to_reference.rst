.. _example_map_structure_to_reference:
.. highlight:: python
.. index::
   single: Tutorial; Structure mapping

Structure mapping
=================

A cluster vector calculation requires that all atoms reside on a fixed
lattice. The interesting properties, on the other hand, are typically
calculated for a structure in which the cell metric and the atoms and have
been allowed to relax. Unless the ideal structures have been saved prior to
relaxation, one is therefore faced with the task of mapping back the relaxed
structure onto the ideal one. This is the purpose of the function
:func:`map_structure_to_reference
<icetdev.tools.map_sites.map_structures_to_reference>`. The function is also
useful to analyze whether the relaxation has gone too far for the cluster
expansion to viable, i.e., whether the ideal structure from which the
relaxation started is not a valid representation of the structure for which
the property has been obtained.

Import modules
--------------

The :func:`map_structure_to_reference
<icetdev.tools.map_sites.map_structures_to_reference>`
function needs to be imported together with some additional functions from `ASE
<https://wiki.fysik.dtu.dk/ase>`_.

.. literalinclude:: ../../../../examples/map_structure_to_reference.py
   :start-after: # Import modules
   :end-before: # End import

Create structures
-----------------

First, a reference structure defining the ideal lattice is created, and a
supercell thereof is scaled and rattled to simulate relaxation in an energy
minimixation.

.. literalinclude:: ../../../../examples/map_structure_to_reference.py
   :start-after: # simulate a relaxed structure.
   :end-before: # Map the

Map the relaxed structure onto an ideal structure
-------------------------------------------------

The structure can now be mapped onto a structure in which all atoms are in
ideal positions. The function returns the ideal structure, as well as the
maximum and average displacement.

.. literalinclude:: ../../../../examples/map_structure_to_reference.py
   :start-after: # Map the "relaxed"
   :end-before: # Map a structure

Structures with vacancies
-------------------------

The function has to be informed if the structure to be mapped contains
vacancies. Two additional keywords should be given, (1) `vacancy_type`, the
chemical symbol that signifies vacancies, and (2) `inert_species`, a list of
elements that are never substituted for a vacancy. The latter can be omitted,
but the mapping is then more likely to fail, at least if the cell has been
relaxed a lot.

In the below example, a Au-Pd-H-vacancy system is created. Vanadium (`'V'`)
will signify vacancies. Further, it is only the hydrogen that may be changed
for a vacancy, so we can set `inert_species = ['Au', 'Pd']`.

.. literalinclude:: ../../../../examples/map_structure_to_reference.py
   :start-after: # Pd and Au share

Source code
-----------

.. container:: toggle

    .. container:: header

       The complete source code is available in
       ``examples/map_structure_to_reference.py``

    .. literalinclude:: ../../../../examples/map_structure_to_reference.py
