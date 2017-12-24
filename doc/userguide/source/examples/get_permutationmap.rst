.. _example_get_permutationmap:
.. highlight:: python
.. index::
   single: Tutorial; Permutation map

Permutation map
===============

The purpose of this example is to demonstrate how to use functionalities, which are built into :program:`iceT`, to extract a permutation map.

Import modules
--------------

To extract permutation maps, only a single :program:`iceT` function is required, namely the :func:`permutation_matrix_from_atoms <icetdev.permutation_map.permutation_matrix_from_atoms>`. Additionally, the `ASE <https://wiki.fysik.dtu.dk/ase>`_ function :func:`ase.build.bulk` will be needed to generate the :class:`Atoms` object.

.. literalinclude:: ../../../../examples/get_permutationmap.py
   :start-after: # Import modules
   :end-before: # Create a prototype

.. _generate-prototype-al-unitcell
Generate prototype Al unit cell
-------------------------------

The next step is to build a prototype :class:`Atoms` object, in the form of a alumininum fcc unit cell.

.. literalinclude:: ../../../../examples/get_permutationmap.py
   :start-after: # Create a prototype
   :end-before: # Generate a permutation

.. _obtain-a-permutation-map
Obtain a permutation map
------------------------

Next, the :func:`permutation_matrix_from_atoms <icetdev.permutation_map.permutation_matrix_from_atoms>` function is called with the :class:`Atoms>` object, which was defined :ref:`above <_generate-prototype-al-unitcell>`, and the neighbor cutoff, in this case 2.0 Å, as arguments. The resulting output is a :class:`list` that containst three elements, namely a :class:`PermutationMap <icetdev.permutation_map.PermutationMap>`, a primitive :class:`Structure <icetdev.structure.Structure>` and a :class:`NeighborList <icetdev.neighbor_list.NeighborList>`.

.. literalinclude:: ../../../../examples/get_permutationmap.py
   :start-after: # (2.0 A).
   :end-before: # Extract the permutated,

Since the verbosity was set :math:`\geq 3`, the following extra information should be printed by :func:`permutation_matrix_from_atoms <icetdev.permutation_map.permutation_matrix_from_atoms>`: ::
   size of atoms_prim 1
   number of positions: 19

.. _extract-permuted-positions
Extract permuted positions
--------------------------

The permutated as well as the indexed and unique positions can be extracted from the :class:`PermutationMap <icetdev.permutation_map.PermutationMap>` object, generated :ref:`earlier <_obtain-a-permutation-map>`, using the methods :meth:`PermutationMap.get_permutated_positions` and
:meth:`PermutationMap.get_indiced_positions`, respectively.

.. literalinclude:: ../../../../examples/get_permutationmap.py
   :start-after: # Extract the permutated,
   :end-before: # Print the permutated,

Print the positions
-------------------

Finally the permutated as well as the indexed and unique positions, obtained in the :ref:`previous section <_extract-permuted-positions>`, are printed using the following syntax:

.. literalinclude:: ../../../../examples/get_permutationmap.py
   :start-after: # Print the permutated,

These lines should give the following (partial) output: ::
   Permutated fractional coordinates
   [ 0.  0.  0.]
   [ 1.  1. -1.] [-1.  1.  1.] [-1. -1.  1.] [ 1. -1. -1.] [ 1. -1.  1.] [-1.  1. -1.]
   ...
   [ 1.  1. -1.] [-1.  1.  1.] [-1. -1.  1.] [ 1. -1. -1.] [ 1. -1.  1.] [-1.  1. -1.]
   Permutated indices and positions
   0 1 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
   1 6 [1, 2, 3, 4, 2, 1, 4, 3, 2, 1, 3, 4, 1, 2, 4, 3, 5, 6, 2, 1, 6, 5, 1, 2, 6, 5, 2, 1, 5, 6, 1, 2, 4, 3, 6, 5, 3, 4, 5, 6, 3, 4, 6, 5, 4, 3, 5, 6]
   ...
   18 6 [2, 1, 4, 3, 1, 2, 3, 4, 1, 2, 4, 3, 2, 1, 3, 4, 6, 5, 1, 2, 5, 6, 2, 1, 5, 6, 1, 2, 6, 5, 2, 1, 3, 4, 5, 6, 4, 3, 6, 5, 4, 3, 5, 6, 3, 4, 6, 5]
   Unique permutated indices and positions
   0 [ 0.  0.  0.]
   1 [-1. -1.  1.]
   ...
   18 [ 0.  1.  0.]

Source code
-----------

.. container:: toggle

    .. container:: header

       The complete source code is available in
       ``examples/get_permutationmap.py``

    .. literalinclude:: ../../../../examples/get_permutationmap.py
