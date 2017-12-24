.. _example_get_neighborlist:
.. highlight:: python
.. index::
   single: Tutorial; Neighbor list

Neighbor list
=============

The purpose of this example is to demonstrate how to use functionalities, which are built into :program:`iceT`, to extract a list of neighbors.

Import modules
--------------

It is necessary to import the :class:`Structure <icetdev.structure.Structure>` class, which  used to store information regarding a specific structure as well as the :func:`get_neighbor_lists <icetdev.neighbor_list.get_neighbor_lists>` function. Additionally, the `ASE <https://wiki.fysik.dtu.dk/ase>`_ function :func:`ase.build.bulk` will be needed to generate the structures.

.. literalinclude:: ../../../../examples/get_neighborlist.py
   :start-after: # Import modules
   :end-before: # Generate an iceT

.. _generate-prototype-al-supercell
Generate prototype Al supercell
-------------------------------

The next step is to build a prototype structure, in the form of a alumininum :math:`2\times2\times2` supercell.

.. literalinclude:: ../../../../examples/get_neighborlist.py
   :start-after: # Generate an iceT
   :end-before: # Construct a list

.. _obtain-neighbor-list
Obtain neighbor list
--------------------

A list of all neighbors within a cutoff distance of 1.5 Å is obtained by calling the :func:`get_neighbor_lists <icetdev.neighbor_list.get_neighbor_lists>` function and providing a :class:`Structure <icetdev.structure.Structure>` object, which was defined :ref:`above <_generate-prototype-al-supercell>`, and a :class:`list` of cutoffs as input arguments. One should note that the output from this functions is a :class:`list` that containst one :class:`NeighborList <icetdev.neighbor_list.NeighborList>` object per cutoff.

.. literalinclude:: ../../../../examples/get_neighborlist.py
   :start-after: # Construct a list
   :end-before: # Loop over all

Print all neighbors
-------------------

All relevant information can be extracted from the :class:`NeighborList <icetdev.neighbor_list.NeighborList>` object, obtained :ref:`earlier <_obtain-neighbor-list> using varius the built-in methods for this class. Here, the goal is to loop over the atoms in the structure and print the indices, offsets and distances for all the neighbors. The first step is to use :meth:`NeighborList.get_neighbors` to obtain the list of neighbors, for the atom with the specified index. While the indeces and offsets are attributes of the objects in the :class:`NeighborList <icetdev.neighbor_list.NeighborList>`, the distances are not. For this reason, it is calculated using the :meth:`NeighborList.get_distance` method, which takes the indices and offsets of the atom and the neighbor as input arguments.

.. literalinclude:: ../../../../examples/get_neighborlist.py
   :start-after: # Loop over all

These lines should give the following (partial) output: ::
   Neighbors of atom with index 0
   1 [ 0.  0. -1.] 1.41421
   1 [ 0.  0.  0.] 1.41421
   2 [ 0. -1.  0.] 1.41421
   2 [ 0.  0.  0.] 1.41421
   3 [ 0. -1.  0.] 1.41421
   3 [ 0.  0. -1.] 1.41421
   4 [-1.  0.  0.] 1.41421
   4 [ 0.  0.  0.] 1.41421
   5 [-1.  0.  0.] 1.41421
   5 [ 0.  0. -1.] 1.41421
   6 [-1.  0.  0.] 1.41421
   6 [ 0. -1.  0.] 1.41421
   ...
   Neighbors of atom with index 7
   1 [ 0.  1.  0.] 1.41421
   1 [ 1.  0.  0.] 1.41421
   2 [ 0.  0.  1.] 1.41421
   2 [ 1.  0.  0.] 1.41421
   3 [ 0.  0.  0.] 1.41421
   3 [ 1.  0.  0.] 1.41421
   4 [ 0.  0.  1.] 1.41421
   4 [ 0.  1.  0.] 1.41421
   5 [ 0.  0.  0.] 1.41421
   5 [ 0.  1.  0.] 1.41421
   6 [ 0.  0.  0.] 1.41421
   6 [ 0.  0.  1.] 1.41421

   fcc has 12 nearest neighbors

Source code
-----------

.. container:: toggle

    .. container:: header

       The complete source code is available in
       ``examples/get_neighborlist.py``

    .. literalinclude:: ../../../../examples/get_neighborlist.py
