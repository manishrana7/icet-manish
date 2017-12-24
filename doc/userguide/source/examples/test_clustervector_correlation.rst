.. _example_test_clustervector_correlation:
.. highlight:: python
.. index::
   single: Tutorial; Cluster vector correlations

Cluster vector correlations
===========================

The purpose of this example is to demonstrate how to use functionalities, which are built into :program:`iceT`, to test the correlation between cluster vectors.

Import modules
--------------

In the present case it is necessary to import two :program:`iceT` classes, namely :class:`ClusterSpace <icetdev.cluster_space.ClusterSpace>` and :class:`Structure <icetdev.structure.Structure>`. In particular, these objects are used to store information regarding a specific cluster space and structure, respectively. Additionally, the `ASE <https://wiki.fysik.dtu.dk/ase>`_ function :func:`ase.db.connect` will be needed to extract structures from a previously prepared database.

.. literalinclude:: ../../../../examples/test_clustervector_correlation.py
   :start-after: # Import modules
   :end-before: # Function for generating random

.. _generate-random-structure
Function: Generate random structure
-----------------------------------

Since the intention is to calculate the correlation between the cluster vectors for a number of structures in a database, a number of functions will be defined that will be called in the loop that is described in the :ref:`last section <_test-the-correlation>`. The first function, :func:`generateRandomStructure`, takes a :class:`ase.Atoms` object, ``atoms_prim``, creates a supercell with the dimensions given by the ``repeat`` aregument, which should be a single or a :class:`list` of integers, and populates it randomly with the elements in the ``subelements`` :class:`list`. At the end the corresponding :class:`ase.Atoms` object is returned.

.. literalinclude:: ../../../../examples/test_clustervector_correlation.py
   :start-after: # Function for generating random
   :end-before: # Function for generating cluster

.. _generate-cv-set
Function: Generate CV set
-------------------------

The second function, :func:`generateCVSet`, generates cluster vectors for ``n`` supercells, which have been constructed from the :class:`ase.Atoms` object, ``atoms_prim`` and randomly populated with the elements in the ``subelements`` :class:`list`, based on the :class:`ClusterSpace <icetdev.cluster_space.ClusterSpace>` object ``clusterspace``. Specifically, this is achieved by first calling the `:func:`generateRandomStructure` function, defined in the :ref:`previous section <_generate-random-structure>`, with the :class:`ase.Atoms` object, ``atoms_prim``, the ``subelements`` :class:`list` and size of the supercell ``repeat`` as input arguments. The resulting configuration is then converted to a :class:`Structure <icetdev.structure.Structure>` object, using the :meth:`Structure.from_atoms` method. This structure is, thereafter, provided as an input argument to the :func:`ClusterSpace.get_cluster_vector` method to extract the cluster vectors.

.. literalinclude:: ../../../../examples/test_clustervector_correlation.py
   :start-after: # Function for generating cluster
   :end-before: # Function for calculating

.. _get-column-correlation
Function: Get column correlation
--------------------------------

The third function, :func:`getColumnCorrelation`, calculates the correlation between the columns ``i`` and ``j`` in the matrix ``cv_matrix``, which should be a :class:`numpy.array` object.

.. literalinclude:: ../../../../examples/test_clustervector_correlation.py
   :start-after: # Function for calculating
   :end-before: # Function for asserting

.. _assert-no-correlation
Function: Check correlation
---------------------------

The fourth function, :func:`assertNoCorrelation`, checks that the correlations between the columns in the matrix ``cvs`` are smaller than the tolerance, ``tol``. This is achieved by calculating the correlations between all pairs of columns with help of the :func:`getColumnCorrelation`, which was defined :ref:`earlier <_get-column-correlation>`. If this value is not smaller than the tolerance, the number of the columns and the correlation is printed to the standard output.

.. literalinclude:: ../../../../examples/test_clustervector_correlation.py
   :start-after: # Function for asserting
   :end-before: # Create a list

.. _setup
Setup
-----

Before starting the test a list of the relevant subelements, namely palladium, hydrogen and vanadinum, of which the latter is used to represent vacancies, is defined. Also, the cutoffs for pairs is set to 2.0 Å and it is decided that the cluster vectors for :math:`8\times8\times8` supercells.

.. literalinclude:: ../../../../examples/test_clustervector_correlation.py
   :start-after: # cutoff distance for
   :end-before: # Test the correlation

.. _test-the-correlation
Test the correlation
--------------------

The functions defined above are now implemented to test the cluster vector correlations for the first 10 structures in the database ``PdHVac-fcc.db``, which should already have been created using the script ``enumerate_structures.py``. At each step of the loop, some basic information regarding the structure is first extracted and printed. Thereafter, a :class:`ClusterSpace <icetdev.cluster_space.ClusterSpace>` object, ``clusterspace``, is initiated using the :class:`ase.Atoms` object ``atoms_row`` as well as ``cutoffs`` :class:`list` and ``subelements`` :class:`list`, given in the :ref:`previous section <_setup>`, as input arguments. After that, the :ref:`previously defined function <_generate-cv-set>` is called with ``n`` set to 20 and the ``atoms_row``, ``subelements`` and ``clusterspace`` as, additional input arguments. The :class:`list` of cluster vectors thus obtained is then fed into the :func:`assertNoCorrelation` function, that was :ref:`discussed earlier <_assert-no-correlation>`, to check if correlations between the columns are smaller than the tolerance, 0.99.

.. literalinclude:: ../../../../examples/test_clustervector_correlation.py
   :start-after: # pregenerated database.

The result of running the script should be the following: ::
   Testing structure: PdH(id=1) with cutoffs [2.0]
   size of atoms 1024. len of cv 6
   Testing structure: PdV(id=2) with cutoffs [2.0]
   size of atoms 1024. len of cv 6
   Testing structure: Pd2VH(id=3) with cutoffs [2.0]
   size of atoms 2048. len of cv 6
   Testing structure: Pd2VH(id=4) with cutoffs [2.0]
   size of atoms 2048. len of cv 6
   Testing structure: Pd3VH2(id=5) with cutoffs [2.0]
   size of atoms 3072. len of cv 6
   Testing structure: Pd3V2H(id=6) with cutoffs [2.0]
   size of atoms 3072. len of cv 6
   Testing structure: Pd3VH2(id=7) with cutoffs [2.0]
   size of atoms 3072. len of cv 6
   Testing structure: Pd3V2H(id=8) with cutoffs [2.0]
   size of atoms 3072. len of cv 6
   Testing structure: Pd3VH2(id=9) with cutoffs [2.0]
   size of atoms 3072. len of cv 6
   Testing structure: Pd3V2H(id=10) with cutoffs [2.0]
   size of atoms 3072. len of cv 6

Source code
-----------

.. container:: toggle

    .. container:: header

       The complete source code is available in
       ``examples/test_clustervector_correlation.py``

    .. literalinclude:: ../../../../examples/test_clustervector_correlation.py
