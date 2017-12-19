.. _example_get_cluster_vector:
.. highlight:: python
.. index::
   single: Tutorial; Cluster vectors

Cluster vectors
===============

The purpose of this example is to demonstrate how to use functionalities, which are built into :program:`iceT`, to construct cluster vectors.

Import modules
--------------

Firstly, one needs to import the class :class:`ClusterSpace <icetdev.ClusterSpace>`, which is used to store information regarding a given cluster space. Additionally, the `ASE <https://wiki.fysik.dtu.dk/ase>`_ function :func:`ase.build.bulk` will be needed to generate the structures.

.. literalinclude:: ../../../../examples/get_cluster_vector.py
   :start-after: # Start import
   :end-before: # End import

.. _`discussed above`
Generate prototype structure
----------------------------

The next step is to build a prototype structure, in the form of a silicon, bulk, unit cell. It is, in addition, decided that the cluster vectors will be created by populating the sites with either silicon or germanium. Also, the cutoffs for pairs, triplets and quadruplets are all set to 5.0 Å.

.. literalinclude:: ../../../../examples/get_cluster_vector.py
   :start-after: # Start setup
   :end-before: # End setup

Create the cluster space
------------------------

The cluster space is created by simpling initiating a :class:`ClusterSpace <icetdev.ClusterSpace>` object and providing the prototype structure, cutoffs and subelements, `discussed above`_, as arguments. Next, the built in :func:`print` function is used to print all relevant information regarding the cluster space in a table format.

.. literalinclude:: ../../../../examples/get_cluster_vector.py
   :start-after: # Start clusterspace
   :end-before: # End clusterspace

Specifically, the final call should produce the following (partial) output: ::
   ------------------------- Cluster Space -------------------------
   subelements: Si Ge
   cutoffs: 5.0 5.0 5.0
   number of orbits: 21
   -----------------------------------------------------------------
   order |  radius  | multiplicity | index | orbit |    MC vector
   -----------------------------------------------------------------
     1   |   0.0000 |        2     |    0  |    0  |    [0]
     2   |   1.1756 |        4     |    1  |    1  |  [0, 0]
     2   |   1.9198 |       12     |    2  |    2  |  [0, 0]
   ...
     4   |   2.5525 |        8     |   20  |   20  | [0, 0, 0, 0]
   -----------------------------------------------------------------


.. _`previous section`
Cluster vectors for Si supercells
---------------------------------

After building a new structure in the form of a :math:`2\times2\times2` pure Si supercell. Using this :class:`ase.Atoms` object as input, the cluster vectors are constructed using the :func:`get_cluster_vector` which is inherent to the :class:`ClusterSpace <icetdev.ClusterSpace>` class. The cluster vectors are printed, as a sequence of tables, with help of the :func:`print` function.

.. literalinclude:: ../../../../examples/get_cluster_vector.py
   :start-after: # Start cluster_vector1
   :end-before: # End cluster_vector1

These lines ought to yield the following result: ::
   [1.0, -1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

Cluster vectors for Si supercells
---------------------------------

Finally, the steps described in the `previous section` are repeated after substituting one of the Si atoms in the supercell with Ge.

.. literalinclude:: ../../../../examples/get_cluster_vector.py
   :start-after: # Start cluster_vector2
   :end-before: # End cluster_vector2

In this case resulting output should be: ::
   [1.0, -0.875, 0.75, 0.75, 0.75, -0.625, -0.625, -0.625, -0.625, -0.625, -0.625, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]

Source code
-----------

.. container:: toggle

    .. container:: header

       The complete source code is available in
       ``examples/get_cluster_vector.py``

    .. literalinclude:: ../../../../examples/get_cluster_vector.py
