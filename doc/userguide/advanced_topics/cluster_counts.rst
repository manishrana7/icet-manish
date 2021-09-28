.. _advanced_topics_get_cluster_count:
.. highlight:: python
.. index::
   single: Advanced topics; Cluster counts

Cluster counts
===============

A cluster vector is essentially a count of clusters averaged and
aggregated in a compact form. Sometimes it may be of interest to study
the underlying cluster counts rather than the cluster vector. In a
binary alloy with species A and B, for example, one may wonder how
many nearest neighbor pairs are A-A, A-B, and B-B, respectively. This
tutorial demonstrates how such information can be extracted with
:program:`icet`.

Import modules
--------------

In addition to the :class:`~icet.core.orbit_list.OrbitList` class, we also
use the `ASE <https://wiki.fysik.dtu.dk/ase>`_ function
:func:`~ase.build.bulk` to generate a prototype structure.

.. literalinclude:: ../../../examples/advanced_topics/get_cluster_counts.py
   :start-after: # Start import
   :end-before: # End import

Generate prototype cell
-----------------------

The next step is to build a prototype supercell. Here, we create a toy
structure based on Ti and W occupying a simple cubic :math:`2\times1\times1`
cell. Moreover, we define a cutoff for pairs at 5 Å. We also keep the
primitive cell, i.e., the :math:`1\times1\times1` simple cubic cell.

.. literalinclude:: ../../../examples/advanced_topics/get_cluster_counts.py
   :start-after: # Start setup
   :end-before: # End setup

Create primitive orbit list
---------------------------

Next, we create an orbit list based on the above primitive structure
and cutoff. To this end, we use the
:class:`~icet.core.orbit_list.OrbitList` class. The orbit list defines
the orbits that are included in the cluster count (every cluster
corresponds to an orbit).

.. literalinclude:: ../../../examples/advanced_topics/get_cluster_counts.py
   :start-after: # pair clusters within the
   :end-before: # Use the primitive


Count clusters
--------------

We are now ready to count the clusters. We do this with the
:func:`~icet.core.orbit_list.get_cluster_counts` function, and
we provide the supercell as an argument.

.. literalinclude:: ../../../examples/advanced_topics/get_cluster_counts.py
   :start-after: # Use the primitive
   :end-before: # Print all of

Print cluster counts
--------------------

Finally, we print information about the orbit and the cluster counts.

.. literalinclude:: ../../../examples/advanced_topics/get_cluster_counts.py
   :start-after: # Print all of

The resulting output should be similar to the following::

  Number of atoms: 2
  Found 3 orbits

  Orbit index: 0
  Number of atoms in cluster: 1
  Cluster radius: 0.0000
  Distances between atoms in the clusters: 
  Ti      1
  W       1

  Orbit index: 1
  Number of atoms in cluster: 2
  Cluster radius: 1.5000
  Distances between atoms in the clusters: 3.0000
  Ti  Ti      2
  Ti  W       2
  W   W       2

  Orbit index: 2
  Number of atoms in cluster: 2
  Cluster radius: 2.1213
  Distances between atoms in the clusters: 4.2426
  Ti  Ti      2
  Ti  W       8
  W   W       2

The orbits are listed by increasing order (singles, pairs, ...) and
increasing size. For every orbit, the number of clusters containing a certain
set of species is printed. Note that these are the numbers that would be used
to construct a cluster vector for this particular structure.


Source code
-----------

.. container:: toggle

    .. container:: header

       The complete source code is available in
       ``examples/get_cluster_counts.py``

    .. literalinclude:: ../../../examples/advanced_topics/get_cluster_counts.py
