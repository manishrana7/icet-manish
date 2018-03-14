.. _tutorial_construct_cluster_expansion:
.. highlight:: python
.. index::
   single: Tutorial; Constructing a cluster expansion

Construction of a cluster expansion
===================================

This exercise is intended to provide an overview of the basic procedure for
constructing a cluster expansion using structures generated :ref:`previously
<tutorial_prepare_reference_data>`. The main objective here is to map out the
miscibility between silver and gold. This means that the property of interest
is the mixing energy. In the end, predicted and target mixing energies are
plotted as a function of the silver concentration.

General preparations
--------------------

A number of `ASE <https://wiki.fysik.dtu.dk/ase>`_ and :program:`icet`
functions are needed in order to set up and train the cluster expansion.
Specifically, :func:`ase.build.bulk` and :func:`ase.db.connect` are required to
build a primitive structure and import relaxed configurations from the database
that was generated :ref:`previously <tutorial_prepare_reference_data>`. The
:program:`icet` classes
:class:`ClusterSpace <icet.core.cluster_space.ClusterSpace>`,
:class:`StructureContainer <icet.core.structure_container.StructureContainer>`,
:class:`Optimizer <icet.fitting.Optimizer>` and
:class:`ClusterExpansion <icet.core.cluster_expansion.ClusterExpansion>`
are used, in sequence, during preparation, compilation and training of the
cluster expansion followed by the extraction of information in the form of
predicted energies from the latter. In the final step, the function
:func:`enumerate_structures
<icet.tools.structure_enumeration.enumerate_structures>` is employed to
generate a large pool of structures for which the mixing energies can be
calculated with help of the finalized cluster expansion. These data are plotted
using the `matplotlib <https://matplotlib.org>`_ library.

.. literalinclude:: ../../../../tutorial/construct_cluster_expansion.py
   :start-after: # import modules
   :end-before: # step 1

Preparation of cluster space
----------------------------

In order to be able to build a cluster expansion, it is first necessary to
create a :class:`ClusterSpace <icet.core.cluster_space.ClusterSpace>` object
based on a prototype structure, here in the form of a bulk gold unit cell. When
initiating the former, one must also provide cutoffs and a list of elements
that should be considered, in this case gold and silver. Here, the cutoffs are
set to 6, 5, and 4 Å for pairs, triplets and quadruplets.

.. literalinclude:: ../../../../tutorial/construct_cluster_expansion.py
   :start-after: # step 1
   :end-before: # step 2

As with many other :program:`icet` objects, it is possible to print core
information in a tabular format by simply calling the :func:`print` function
with the instance of interest as input argument. For the case at hand, the
output should be the following::

  -------------------------- Cluster Space ---------------------------
   subelements: Ag Au
   cutoffs: 6.0000 5.0000 4.0000
   total number of orbits: 14
   number of orbits by order: 0= 1  1= 1  2= 4  3= 7  4= 1
  --------------------------------------------------------------------
  index | order |   size   | multiplicity | orbit index |  MC vector  
  --------------------------------------------------------------------
     0  |   0   |   0.0000 |        1     |      -1
     1  |   1   |   0.0000 |        1     |       0     |    [0]
     2  |   2   |   1.4425 |        6     |       1     |  [0, 0]
     3  |   2   |   2.0400 |        3     |       2     |  [0, 0]
     4  |   2   |   2.4985 |       12     |       3     |  [0, 0]
     5  |   2   |   2.8850 |        6     |       4     |  [0, 0]
     6  |   3   |   1.6657 |        8     |       5     | [0, 0, 0]
     7  |   3   |   1.8869 |       12     |       6     | [0, 0, 0]
     8  |   3   |   2.0168 |       24     |       7     | [0, 0, 0]
     9  |   3   |   2.3021 |       24     |       8     | [0, 0, 0]
    10  |   3   |   2.4967 |       24     |       9     | [0, 0, 0]
    11  |   3   |   2.7099 |       24     |      10     | [0, 0, 0]
    12  |   3   |   2.8850 |        8     |      11     | [0, 0, 0]
    13  |   4   |   1.7667 |        2     |      12     | [0, 0, 0, 0]
  --------------------------------------------------------------------

Compilation of structure container
----------------------------------

Once a :class:`ClusterSpace <icet.core.cluster_space.ClusterSpace>` has been
prepared, the next step is to compile a :class:`StructureContainer
<icet.core.structure_container.StructureContainer>`. To this end, we first
initialize an empty :class:`StructureContainer
<icet.core.structure_container.StructureContainer>` and then add the
tructures from the database prepared previously including for each structure
the mixing energy in the property dictionary.

.. literalinclude:: ../../../../tutorial/construct_cluster_expansion.py
   :start-after: # step 2
   :end-before: # step 3

By calling the :func:`print` function with the :class:`StructureContainer
<icet.core.structure_container.StructureContainer>` as input argument, one
obtains the following result::

  ----------------------- Structure Container ------------------------
  Total number of structures: 137
  --------------------------------------------------------------------
  index |       user_tag        | natoms | chemical formula |  energy
  --------------------------------------------------------------------
     0  | 0                     |     1  | Ag               |    0.000
     1  | 1                     |     1  | Au               |    0.000
     2  | 2                     |     2  | AgAu             |   -0.010
     3  | 3                     |     2  | AgAu             |   -0.011
     4  | 4                     |     3  | Ag2Au            |   -0.008
     5  | 5                     |     3  | AgAu2            |   -0.008
     6  | 6                     |     3  | Ag2Au            |   -0.009
     7  | 7                     |     3  | AgAu2            |   -0.011
     8  | 8                     |     3  | Ag2Au            |   -0.011
     9  | 9                     |     3  | AgAu2            |   -0.010
   ...
   127  | 127                   |     6  | Ag2Au4           |   -0.010
   128  | 128                   |     6  | AgAu5            |   -0.006
   129  | 129                   |     6  | Ag5Au            |   -0.006
   130  | 130                   |     6  | Ag4Au2           |   -0.009
   131  | 131                   |     6  | Ag4Au2           |   -0.009
   132  | 132                   |     6  | Ag3Au3           |   -0.011
   133  | 133                   |     6  | Ag3Au3           |   -0.012
   134  | 134                   |     6  | Ag2Au4           |   -0.011
   135  | 135                   |     6  | Ag2Au4           |   -0.011
   136  | 136                   |     6  | AgAu5            |   -0.007
  --------------------------------------------------------------------

Training of parameters
----------------------

Since the :class:`StructureContainer
<icet.core.structure_container.StructureContainer>` object created in the
previous section, contains all the information required for constructing a
cluster expansion, the next step is to train the parameters, i.e. to fit the
*effective cluster interactions* (ECIs) using the target data. More precisely,
the goal is to achieve the best possible agreement with set of training
structures, which represent a subset of all the structures in the
:class:`StructureContainer
<icet.core.structure_container.StructureContainer>`. In practice, this is a
two step process that involves the initiation of an :class:`Optimizer
<icet.fitting.Optimizer>` object with the a list of target properties
produced by the :func:`StructureContainer.get_fit_data
<icet.core.structure_container.StructureContainer.get_fit_data>` method as
input argument.

.. literalinclude:: ../../../../tutorial/construct_cluster_expansion.py
   :start-after: # step 3
   :end-before: # step 4

The training process is started by calling the :func:`Optimizer.train
<icet.fitting.Optimizer.train>` method. Once it is finished, the results
can be displayed by providing the :class:`Optimizer
<icet.fitting.Optimizer>` object to the :func:`print` function, which gives
the output shown below::

  ===================== Optimizer ======================
  fit_method                : least-squares
  number_of_target_values   : 137
  number_of_parameters      : 14
  rmse_training             : 0.000177534
  rmse_test                 : 0.000216184
  training_size             : 102
  test_size                 : 35
  ======================================================

Comparison with target data
---------------------------

At this point, the task of constructing the cluster expansion is almost
complete. In fact, the only step that remains is to tie the parameter values
obtained from the optimization to the cluster space. This is achieved through
the initiation of a :class:`ClusterExpansion
<icet.core.cluster_expansion.ClusterExpansion>` object with the previously
created :class:`ClusterSpace <icet.core.cluster_space.ClusterSpace>`
instance together with the list of parameters, available via the
:class:`Optimizer.parameters <icet.fitting.Optimizer.parameters>` attribute,
as input arguments. Using this cluster expansion, it is now possible to predict
the properties of any configuration that is based on a supercell of the
primitive structure. A convenient and illustrative ways of checking the
accuracy of this model is to plot the predicted and target properties for all
available structures. First, however, it is necessary to gather the information
required for constructing such a plot, specifically by, once again, looping
over all rows in the ``structures.db`` database and by compiling the silver
concentration as well as the target and predicted mixing energies into a list.
Note that the latter of the three values is calculated by calling the
:func:`ClusterExpansion.predict
<icet.core.cluster_expansion.ClusterExpansion.predict>` method with the
:class:`ASE Atoms` object that represents the present structure as input
argument. Once the list is complete the predicted and target mixing energies
are plotted as functions of the concentration of silver atoms.

.. literalinclude:: ../../../../tutorial/construct_cluster_expansion.py
   :start-after: # step 3
   :end-before: # step 4

Be aware of the fact that the diagram, which is shown below is saved as an
image file, ``mixing-energy-comparison.png``.

.. figure:: ../../../../tutorial/mixing-energy-comparison.png

  Predicted (crosses) and target (open circles) mixing energies versus silver
  concentration for the structures used in the construction of the cluster
  expansion.

Predictions for enumerated structures
-------------------------------------

Since the construction of a cluster expansion is, generally, only the first
part of a study, equivalent to the development of a model, a natural
continuation is to use the latter to predict the property of interest for other
structures, e.g., generated by systematic enumeration.

In principle, the procedure is the same as in the previous section, with the
exception that the structures are not read from a database but rather generated
with help of the :func:`enumerate_structures
<icet.tools.structure_enumeration.enumerate_structures>` function, which
requires three input arguments. This includes the :class:`ASE Atoms` object
that represents the primitive structure and the list of elements, defined
above, in addition to a list of sizes.

Also, the thermodynamically relevant convex hull for the predicted
structures can be calculated with the help of a :class:`ConvexHull
<icet.tools.convex_hull.ConvexHull>` object.

.. literalinclude:: ../../../../tutorial/construct_cluster_expansion.py
   :start-after: # step 5
   :end-before: # step 6

The result is a file, ``mixing-energy-predicted.png``, that contains an image
of the following diagram:

.. figure:: ../../../../tutorial/mixing-energy-predicted.png

  Predicted mixing energies versus silver concentration for a set of
  systematically enumerated structures.

Extract structures of special interest
--------------------------------------

Now that the energy has been computed for thousands of structures, it is
interesting to identify the ground state structures. Since even small errors in
the cluster expansion can cause the energetic ordering of the structures to
change, it is advisable to pick out structures that are energetically close to
the convex hull. This can be done with
:func:`extract_low_energy_structures
<icet.tools.convex_hull.ConvexHull.extract_low_energy_structures>`.

.. literalinclude:: ../../../../tutorial/construct_cluster_expansion.py
   :start-after: # step 6

All structures that are within 0.5 meV/atom from the convex hull are now
available and can, for example, be fed into another cluster expansion once
their energy has been calculated using the reference method of choice.

Source code
-----------

.. container:: toggle

    .. container:: header

       The complete source code is available in
       ``tutorial/construct_cluster_expansion.py``

    .. literalinclude:: ../../../../tutorial/construct_cluster_expansion.py
