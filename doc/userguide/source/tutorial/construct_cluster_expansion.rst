.. _tutorial_construct_cluster_expansion:
.. highlight:: python
.. index::
   single: Tutorial; Constructing a cluster expansion

Construction of a cluster expansion
===================================

This exercise is intended to provide an overview of the basic procedure for constructing a cluster expansion for a prototype system, using structures found in the databased generated in the preceeding part of the tutorial. In particular, the main objective is to determine the miscibility between silver and gold. This means that the property of interest is the mixing energy, which, therefore, will be used as the cluster property. In the end, the predicted and target mixing energies are plotted as functions of the silver atom concentration.

Importation of modules
----------------------

 A number of `ASE <https://wiki.fysik.dtu.dk/ase>`_ and :program:`icet` functions are needed in order to set up and train the cluster expansion. Specifically, :func:`ase.build.bulk` and :func:`ase.db.connect` are required to build a primitive structure and import relaxed configurations from a database, which was generated in a preceeding part of the tutorial. The :program:`icet` classes :class:`icetdev.ClusterSpace`, :class:`icetdev.StructureContainer`, :class:`icetdev.Optimizer` and :class:`icetdev.ClusterExpansion`, meanwhile, are used, in sequence, during the preparation, compilation and training of the cluster expansion followed by the extraction of information, in the form of predicted energies, from the latter. In the final step, the function :func:`icetdev.enumeration.enumerate_structures` is implemented to generate a large pool of structures for which the mixing energies can be calculated with help of the finalised cluster expansion. These predictions are plotted with defined in the :mod:`matplotlib.pyplot` module, which forms part of the `Matplotlib <https://matplotlib.org>`_ plotting library.

.. literalinclude:: ../../../../tutorial/construct_cluster_expansion.py
   :start-after: # import modules
   :end-before: # step 1

Preparation of cluster space
----------------------------

In order to be able to build a cluster expansion, it is first necessary to create a :class:`icetdev.CluserSpace` object based on a primitive structure, here in the form of a bulk gold unit cell. When initiating the former, one must also provide a :class:`list` of cutoffs and another with the subelements that should be considered, in this case gold and silver. At present, the former indicates that pairs, triplets and quadruplets smaller than 6.0 Å, 5.0 Å and 4.0 Å shall be included.

.. literalinclude:: ../../../../tutorial/construct_cluster_expansion.py
   :start-after: # step 1
   :end-before: # step 2

Like many other :program:`icet` objects, it is possible to print the information contained within, in a tabular format, by simply calling the :func:`print` function with the instance of interest as input argument. For the case at hand, the output should be the following::

  ------------------------- Cluster Space -------------------------
   subelements: Ag Au
   cutoffs: 6.0 5.0 4.0
   number of orbits: 13
  -----------------------------------------------------------------
  order |  radius  | multiplicity | index | orbit |    MC vector
  -----------------------------------------------------------------
    1   |   0.0000 |        1     |    0  |    0  |    [0]
    2   |   1.4425 |        6     |    1  |    1  |  [0, 0]
    2   |   2.0400 |        3     |    2  |    2  |  [0, 0]
    2   |   2.4985 |       12     |    3  |    3  |  [0, 0]
    2   |   2.8850 |        6     |    4  |    4  |  [0, 0]
    3   |   1.6657 |        8     |    5  |    5  | [0, 0, 0]
    3   |   1.8869 |       12     |    6  |    6  | [0, 0, 0]
    3   |   2.0168 |       24     |    7  |    7  | [0, 0, 0]
    3   |   2.3021 |       24     |    8  |    8  | [0, 0, 0]
    3   |   2.4967 |       24     |    9  |    9  | [0, 0, 0]
    3   |   2.7099 |       24     |   10  |   10  | [0, 0, 0]
    3   |   2.8850 |        8     |   11  |   11  | [0, 0, 0]
    4   |   1.7667 |        2     |   12  |   12  | [0, 0, 0, 0]
  -----------------------------------------------------------------

Compililation of structure container
------------------------------------

Once a :class:`icetdev.ClusterSpace` has been prepared, the next step is compile a :class:`icetdev.StructureContainer`. The initation of the latter requires, in addition to the former, a :class:`list` of :class:`ase.Atoms` objects as well as a :class:`list` of target properties, with one item for each structure, as input arguments. Since the property of interest, in the present case, is the mixing energy, it is necessary to start by determining the energy associated with each subelement. This is achieved by extracting the energies, per atom, for each a structure that contains a single atom, with the correct identity. Once these, reference, energies have been obtained, it is possible loop over all entries in the database and calculate the mixing energies before adding these to the :class:`list` of properties and the :class:`ase.Atoms` object to the :class:`list` of structures. Once these lists have been compiled, the :class:`icetdev.StructureContainer` object is initiated.

.. literalinclude:: ../../../../tutorial/construct_cluster_expansion.py
   :start-after: # step 2
   :end-before: # step 3

By calling the :func:`print` function with the :class:`icetdev.StructureContainer` as input argument, one obtains the following result::

  ---------------- Structure Container -----------------
  Total number of structures: 137
  ------------------------------------------------------
  index |           user_tag           | natoms | energy
  ------------------------------------------------------
     0  | None                         |     1  | 0.000
     1  | None                         |     1  | 0.000
     2  | None                         |     2  | -0.010
     3  | None                         |     2  | -0.011
     4  | None                         |     3  | -0.008
     5  | None                         |     3  | -0.008
     6  | None                         |     3  | -0.009
     7  | None                         |     3  | -0.011
     8  | None                         |     3  | -0.011
     9  | None                         |     3  | -0.010
   ...
   127  | None                         |     6  | -0.010
   128  | None                         |     6  | -0.006
   129  | None                         |     6  | -0.006
   130  | None                         |     6  | -0.009
   131  | None                         |     6  | -0.009
   132  | None                         |     6  | -0.011
   133  | None                         |     6  | -0.012
   134  | None                         |     6  | -0.011
   135  | None                         |     6  | -0.011
   136  | None                         |     6  | -0.007
  ------------------------------------------------------

Training of parameters
----------------------

Since the :class:`icetdev.StructureContainer` object created in the previous section, contains all the information required to construct a cluster expansion, the next step is *train the parameters*. In this context the term *training* corresponds to fitting the *effective cluster interactions* (ECI:s) to the values on the target properties. More precisely, the goal is to achieve the best possible agreement with the set of *training* configurations, which is a subset of all the structures in the :class:`icetdev.StructureContainer`. In practice, this is a two step process that involves the initiation of an :class:`icetdev.Optimizer` object, with the a list of target properties, produced by the :func:`icetdev.StructureContainer.get_fit_data` method, as input argument.

.. literalinclude:: ../../../../tutorial/construct_cluster_expansion.py
   :start-after: # step 3
   :end-before: # step 4

When the :func:`icetdev.Optimizer.train` method is called the *training* process starts. Once it is finished, the lines below should have been printed::

  ------------------Training-------------------
  Fit Method least-squares, N_params 14
  Train size 102, Test size 35
  Train RMSE  0.00018
  Test  RMSE  0.00022
  --------------------Done---------------------

Now that the solution to the optimisation problem has been found, the results can displayed by providing the :class:`icetdev.Optimizer` object to the :func:`print` function, which gives the output shown below::

  fit method             : least-squares
  number of target values : 137
  number of parameters   : 14
  training set size      : 102
  testing set size       : 35

Comparison with target data
---------------------------

At this point, the task of constructing the cluster expansion is almost complete. In fact, the only step that remains is to tie the parameter values obtained from the optimisation to the cluster space. This is, more precisely, achieved through the initiation of a :class:`icetdev.ClusterExpansion` object, with the previously created :class:`icetdev.ClusterSpace` instance together with the :class:`list` of parameters, represented by the :class:`icetdev.Optimizer.parameters` attribute, as input arguments. Using this cluster expansion, it is now possible to predict the properties of any configuration that is based on a supercell of the primitve structure. The most convenient and illustrative ways of checking the accuracy of this model, is to plot the predicted and target properties for all available structures. Firstly, however, it is necessary to gather the information required for constructing such a plot, specifically by, once again, looping over all rows in the ``structures.db`` database and collecting the silver atom concentration, as well as the target and predicted mixing energies into a :class:`list`. Note that the latter of the three values is calculated by calling the :func:`icetdev.ClusterExpansion.predict` with the :class:`ase.Atoms` object that represents the present structure as input argument. Once the list is complete the predicted and target mixing energies are plotted as functions of the concentration of silver atoms.

.. literalinclude:: ../../../../tutorial/construct_cluster_expansion.py
   :start-after: # step 3
   :end-before: # step 4

Be aware of the fact that the diagram, which is shown below, will not automatically appear since it is saved as an image file, ``mixing-energy-comparison.pdf``.

.. figure:: ../../../../tutorial/mixing-energy-comparison.pdf

  Plot of the predicted (crosses) and target (open cirecles) mixing energies versus the concentration of silver atoms for the structures used in the construction of the cluster expansion.

Predictions for enumerated structures
-------------------------------------

Since the construction of a cluster expansion is, generally, only the first part of a study, equivalent to the development of a model, a natural continuation is to use it to predict the properties for various randomly generated structures. In principle, the procedure is the same as in the previous section, with the exeception that the structures are not read from a database but rather generated with help of the :func:`icetdev.enumeration.enumerate_structures` function, which requires three input arguments. This includes the :class:`ase.Atoms` object that represents the primitive structure and the :class:`list` of subelements, defined above, in addition to a :class:`list` of sizes, which in the present case ranges from 1 to 16 atoms per cell.

.. literalinclude:: ../../../../tutorial/construct_cluster_expansion.py
   :start-after: # step 5

The result is a file, ``mixing-energy-predicted.pdf``, that contains an image of the following diagram:

.. figure:: ../../../../tutorial/mixing-energy-predicted.pdf

  Plot of the predicted mixing energies versus the concentration of silver atoms for a set of randomly generated, enumerated, structures.

Source code
-----------

.. container:: toggle

    .. container:: header

       The complete source code is available in
       ``tutorial/construct_cluster_expansion.py``

    .. literalinclude:: ../../../../tutorial/construct_cluster_expansion.py
