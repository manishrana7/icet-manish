.. _example_get_native_clusters:
.. highlight:: python
.. index::
   single: Tutorial; Native clusters

Native clusters
===============

The purpose of this example is to demonstrate how to use functionalities, which are built into :program:`iceT`, to extract native clusters.

Import modules
--------------

In the present case it is necessary to import two :program:`iceT` classes, namely :class:`ClusterSpace <icetdev.cluster_space.ClusterSpace>` and :class:`Structure <icetdev.structure.Structure>`. In particular, these objects are used to store information regarding a specific cluster space and structure, respectively. Additionally, the `ASE <https://wiki.fysik.dtu.dk/ase>`_ function :func:`ase.build.bulk` will be needed to generate the structures.

.. literalinclude:: ../../../../examples/get_native_clusters.py
   :start-after: # Start import
   :end-before: # End import

.. _`discussed above`
Generate prototype structure
----------------------------

The next step is to build a prototype structure, in the form of a silicon, bulk, unit cell. It is, in addition, decided that the cluster vectors will be created by populating the sites with either silicon or germanium. Also, the cutoffs for pairs set to 10.0 Å.

.. literalinclude:: ../../../../examples/get_native_clusters.py
   :start-after: # Start setup
   :end-before: # End setup

Create the cluster space
------------------------

The cluster space is created by simpling initiating a :class:`ClusterSpace <icetdev.ClusterSpace>` object and providing the prototype structure, cutoffs and subelements, `discussed above`_, as arguments.

.. literalinclude:: ../../../../examples/get_native_clusters.py
   :start-after: # Start clusterspace
   :end-before: # End clusterspace

.. _`previous section`
Cluster vectors for Si supercells
---------------------------------

After building a new structure in the form of a :math:`2\times2\times2` pure Si supercell. Using this :class:`ase.Atoms` object as input, the cluster vectors are constructed using the :func:`get_native_clusters` which is inherent to the :class:`ClusterSpace <icetdev.ClusterSpace>` class. The cluster vectors are printed, as a sequence of tables, with help of the :func:`print` function.

.. literalinclude:: ../../../../examples/get_native_clusters.py
   :start-after: # Start cluster_vector1
   :end-before: # End cluster_vector1

These lines ought to yield the following result: ::
   [1.0, -1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

Cluster vectors for Si supercells
---------------------------------

Finally, the steps described in the `previous section` are repeated after substituting one of the Si atoms in the supercell with Ge.

.. literalinclude:: ../../../../examples/get_native_clusters.py
   :start-after: # Start cluster_vector2
   :end-before: # End cluster_vector2

In this case resulting output should be: ::
   [1.0, -0.875, 0.75, 0.75, 0.75, -0.625, -0.625, -0.625, -0.625, -0.625, -0.625, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]

Source code
-----------

.. container:: toggle

    .. container:: header

       The complete source code is available in
       ``examples/get_native_clusters.py``

    .. literalinclude:: ../../../../examples/get_native_clusters.py
