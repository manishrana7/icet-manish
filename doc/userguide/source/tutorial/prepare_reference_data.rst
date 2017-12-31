.. _tutorial_prepare_reference_data:
.. highlight:: python
.. index::
   single: Tutorial; Preparing reference data

Preparation of reference data
=============================

Throughout the tutorial we will be using some reference data for training as
well as comparison. The present section provides a short description of the
code to generate these data.


Importation of modules
----------------------

As will be seen later, the several `ASE <https://wiki.fysik.dtu.dk/ase>`_ functions are required to generate a database with meaningful information. Specifically, :func:`ase.db.connect` and :func:`ase.build.bulk` are needed to be initiate the database and the primitive structure, respectively, while :func:`ase.calculators.emt.EMT` and :func:`ase.optimize.BFGS` will be used relax and optimise the final structures. To construct the latter, however, the :program:`icet` function :func:`icetdev.enumeration.enumerate_structures` will be implemented.

.. literalinclude:: ../../../../tutorial/prepare_reference_data.py
   :start-after: # import modules
   :end-before: # step 1


Preparation
-----------

The first step is to initiate the database, which will be called ``structures.db``, as well as the primitive structure, in the form of an gold bulk unit cell. Additionally, it is decided that the enumerated structures, created in the next step, will be randomly populated with gold and silver atoms.

.. literalinclude:: ../../../../tutorial/prepare_reference_data.py
   :start-after: # step 1
   :end-before: # step 2


Enumeration and relaxation
--------------------------

The second and final step is to generate, relax and then add the final structures to the database. This is achieved by looping over the :class:`ase.Atoms` instances obtained by calling the :func:`icetdev.enumeration.enumerate_structure` function with the primitive structure and the :class:`list` of subelements, defined earlier, as well as the :class:`list` sizes, which corresponds to the number of atoms per cell, as input arguments. Provided that the given object is not already present in the database, the original positions of all atoms are recorded. Next, attaching a :func:`ase.calculators.emt.EMT` calculator is attached and then the structure is relaxed, using the quasi-Newton type minimisation algorithm :class:`ase.optimize.BFGS`, until all forces are smaller than 0.01 eV/atom. Lastly, the relaxed structure is added to the database, together with additional information, regarding the original atomic positions.

.. literalinclude:: ../../../../tutorial/prepare_reference_data.py
   :start-after: # step 2

When running this script, the :func:`ase.optimize.BFGS.run` function prints the time, energy and maximum forces after each step of the structure relaxation. In particular, the output should be a long list of entries, similar to the following::

  ...
        Step     Time          Energy         fmax
  BFGS:    0 17:12:42       -0.024699        0.1726
  BFGS:    1 17:12:42       -0.025870        0.1556
  BFGS:    2 17:12:42       -0.031735        0.0157
  BFGS:    3 17:12:42       -0.031748        0.0107
  BFGS:    4 17:12:42       -0.031755        0.0073
        Step     Time          Energy         fmax
  BFGS:    0 17:12:42       -0.020691        0.0691
  BFGS:    1 17:12:42       -0.020822        0.0636
  BFGS:    2 17:12:42       -0.021549        0.0004
        Step     Time          Energy         fmax
  BFGS:    0 17:12:42       -0.023467        0.0000
  ...


Source code
-----------

.. container:: toggle

    .. container:: header

       The complete source code is available in
       ``tutorial/prepare_reference_data.py``

    .. literalinclude:: ../../../../tutorial/prepare_reference_data.py
