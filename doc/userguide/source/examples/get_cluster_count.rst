.. _example_get_cluster_count:
.. highlight:: python
.. index::
   single: Tutorial; Cluster counts

Cluster counts
===============

The purpose of this example is to demonstrate how to use functionalities, which are built into :program:`iceT`, to construct cluster counts.

Import modules
--------------

Firstly, one needs to import the classes :class:`ClusterCounts <icetdev.cluster_counts.ClusterCounts>` and :class:`Structure <icetdev.structure.Structure>`, which are used to store information regarding the cluster counts and :program:`iceT` structures, respectively, in addition to the functions :func:`create_orbit_list <icetdev.orbit_list.create_orbit_list>` and :func:`__get_primitive_structure <icetdev.permutation_map.__get_primitive_structure>`. In particular, the latter will be used to generate a primitive structure, for which an orbitlist will be created, using the former. Also, the `ASE <https://wiki.fysik.dtu.dk/ase>`_ function :func:`ase.build.bulk` is needed to generate the prototype structure.

.. literalinclude:: ../../../../examples/get_cluster_count.py
   :start-after: # Import modules
   :end-before: # Create a titanium,

.. _generate-prototype-ti-supercell
Generate a prototype Ti supercell
---------------------------------

The next step is to build a prototype supercell, in the form of a, :math:`2\times2\times1`, single layer sheet that is randomly populated with titanium and tungsten atoms. Moreover, the cutoffs for pairs is set to 4 Å.

.. literalinclude:: ../../../../examples/get_cluster_count.py
   :start-after: # sites with W atoms.
   :end-before: # Create the orbit

.. _create-a-primitive-orbit-list
Create a primitive orbit list
-----------------------------

As will be seen :ref:`later <_calculate-cluster-counts>`, the orbit list corresponding to the primitive structure of the :class:`ase.Atoms` object, defined :ref:`earlier < _generate-prototype-ti-supercell>`, is required to calculate the cluster counts. More precisely, this is achieved by first calling the :func:`__get_primitive_structure <icetdev.permutation_map.__get_primitive_structure>` followed by the :meth:`Structure.from_atoms` method. Thereafter, the primitive structure thus obtained is provided as input to :func:`create_orbit_list <icetdev.orbit_list.create_orbit_list>`, together with the :class:`list` of cutoff distances for the clusters.

.. literalinclude:: ../../../../examples/get_cluster_count.py
   :start-after: # Create the orbit
   :end-before: # Use the primitive


.. _count-the-clusters
Count the clusters
------------------

To count the number of clusters a :class:`ClusterSpace <icetdev.ClusterSpace>` object is first initiated. Thereafter, a :class:`Structure <icetdev.structure.Structure>` is generated from the supercell that was definied :ref:`earlier <_generate-prototype-ti-supercell>` with help of the :meth:`Structure.from_atoms` method. To calculate the number of clusters, this structure together with the primitive orbit list that was created in the :ref:`previous section <_create-a-primitive-orbit-list>` are given as input arguments to the :meth:`ClusterSpace.count_clusters` method.

.. literalinclude:: ../../../../examples/get_cluster_count.py
   :start-after: # Use the primitive
   :end-before: # Print all of

Print the cluster counts
------------------------

Finally, the cluster counts, that were calculated by the :ref:`previous code block <_count-the-clusters>`, are printed using the :meth:`ClusterSpace.print` method.

.. literalinclude:: ../../../../examples/get_cluster_count.py
   :start-after: # Print all of

The resulting output should be similar to the following: ::
   number of atoms 4
   Found 3 clusters
   3.43  :: 0 0 1.715
   ==============
   Ti Ti 2
   W W 2
   Total: 4

   2.97047  :: 0 0 1.48523
   ==============
   Ti W 8
   Total: 8

    :: 0 0
   ==============
   Ti 2
   W 2
   Total: 4

Source code
-----------

.. container:: toggle

    .. container:: header

       The complete source code is available in
       ``examples/get_cluster_count.py``

    .. literalinclude:: ../../../../examples/get_cluster_count.py
