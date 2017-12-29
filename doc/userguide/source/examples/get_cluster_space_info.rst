.. _example_get_cluster_space_info:
.. highlight:: python
.. index::
   single: Tutorial; Cluster space information

Cluster space information
=========================

The purpose of this example is to demonstrate how to use functionalities, which are built into :program:`icet`, to obtain basic information about a cluster
space.

Import modules
--------------

Firstly, one needs to import the class :class:`ClusterSpace <icetdev.ClusterSpace>`, which is used to store information regarding a given cluster space. In addition, one needs two tools, namely :func:`get_singlet_info <icetdev.cluster_space.get_singlet_info>` and :func:`view_singlets <icetdev.cluster_space.view_singlets>`, for extracting specific details regarding the singlet clusters. Moreover, the `ASE <https://wiki.fysik.dtu.dk/ase>`_ function :func:`ase.build.bulk` will be needed to construct the structures.

.. literalinclude:: ../../../../examples/get_cluster_space_info.py
   :start-after: # Import modules
   :end-before: # Create a prototype

Generate prototype Re structure
-------------------------------

The next step is to build a prototype structure, in the form of a rhenium, bulk, unit cell. It is, in addition, decided that the cluster space will be created through substitution of the Re atoms with the following, additional, subelements; titanium, tungsten or molybdenum. Also, the cutoffs for pairs, triplets and quadruplets are set to 10.0 Å , 7.0 Å and 5.0 Å, respectively.

.. literalinclude:: ../../../../examples/get_cluster_space_info.py
   :start-after: # triplets (7.0 A)
   :end-before: # Generate and print

Create the cluster space
------------------------

The cluster space is created by simpling initiating a :class:`ClusterSpace <icetdev.ClusterSpace>` object and providing the prototype structure, cutoffs and subelements, discussed above, as arguments. Next, the :meth:`ClusterSpace.print` method is used to print all relevant information regarding the cluster space in a table format.

.. literalinclude:: ../../../../examples/get_cluster_space_info.py
   :start-after: # Generate and print
   :end-before: # Extract and print

Specifically, the final call should produce the following (partial) output::

  ------------------------- Cluster Space -------------------------
   subelements: Ti Mo W Re
   cutoffs: 10.0 7.0 5.0
   number of orbits: 3768
  -----------------------------------------------------------------
  order |  radius  | multiplicity | index | orbit |    MC vector
  -----------------------------------------------------------------
    1   |   0.0000 |        2     |    0  |    0  |    [0]
    1   |   0.0000 |        2     |    1  |    0  |    [1]
    1   |   0.0000 |        2     |    2  |    0  |    [2]
    2   |   1.3699 |        6     |    3  |    1  |  [0, 0]
  ...
    4   |   2.7094 |       12     | 3767  |  135  | [2, 2, 2, 2]
  -----------------------------------------------------------------

Information regarding singlets
------------------------------

Additonal information regarding the singlets is extracted with help of the :func:`get_singlet_info <icetdev.cluster_space.get_singlet_info>` function. Afterwards, the corresponding clusters are printed by calling :func:`view_singlets <icetdev.cluster_space.view_singlets>`. One should note that both functions take the prototype :class:`ase.Atoms` object created earlier as input argument.

.. literalinclude:: ../../../../examples/get_cluster_space_info.py
   :start-after: # Extract and print

These lines ought to yield the following result::

  Singlets:
   orbit_index            : 0
   sites                  : [[0 : [ 0.  0.  0.]], [1 : [ 0.  0.  0.]]]
   multiplicity           : 2
   representative_site    : [0 : [ 0.  0.  0.]]

Source code
-----------

.. container:: toggle

    .. container:: header

       The complete source code is available in
       ``examples/get_cluster_space_info.py``

    .. literalinclude:: ../../../../examples/get_cluster_space_info.py
