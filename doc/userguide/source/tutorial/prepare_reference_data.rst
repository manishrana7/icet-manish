.. _tutorial_prepare_reference_data:
.. highlight:: python
.. index::
   single: Tutorial; Preparing reference data

Preparing reference data
========================

Throughout the tutorial we will be using some reference data for training as
well as comparison. The present section provides a short description of the
code to generate these data.


Import modules
------------------------

As will be seen later, the several `ASE <https://wiki.fysik.dtu.dk/ase>`_ functions are required to generate a database with meaningful information. Specifically, :func:`ase.db.connect` and :func:`ase.build.bulk` are needed to be initiate the database and the primitive structure, respectively, while :func:`ase.calculators.emt.EMT` and :func:`ase.optimize.BFGS` will be used relax and optimise the final structures. To construct the latter, however, the :program:`iceT` function :func:`icetdev.enumeration.enumerate_structures` will be implemented.

.. literalinclude:: ../../../../tutorial/prepare_reference_data.py
   :start-after: # import modules
   :end-before: # step 1


Preparation
------------------------

The first step is to initiate the database, which will be called `structures.db`, as well as the primitive structure, in the form of an gold bulk unit cell. Additionally, it is decided that the enumerated structures, created in the next step, will be randomly populated with gold and silver atoms.

.. literalinclude:: ../../../../tutorial/prepare_reference_data.py
   :start-after: # step 1
   :end-before: # step 2


Enumeration and relaxation
--------------------------

The second and final step is to generate, relax and then add the final structures to the database. This is achieved by looping over the :class:`ase.Atoms` instances obtained by calling the :func:`icetdev.enumeration.enumerate_structure` function with the primitive structure and the :class:`list` of subelements, defined earlier, as well as the :class:`list` sizes, which corresponds to the number of atoms per cell, as input arguments. Provided that the given object is not already present in the database, the original positions of all atoms are recorded. Next, attaching a :func:`ase.calculators.emt.EMT` calculator is attached and then the structure is relaxed, using the quasi-Newton type minimisation algorithm:func:`ase.optimize.BFGS`, until all forces are smaller than 0.01 eV/atom. Lastly, the relaxed structure is added to the database, together with the, additional, information regarding the original atomic positions.

.. literalinclude:: ../../../../tutorial/prepare_reference_data.py
   :start-after: # step 2


Source code
-----------

.. container:: toggle

    .. container:: header

       The complete source code is available in
       ``tutorial/prepare_reference_data.py``

    .. literalinclude:: ../../../../tutorial/prepare_reference_data.py
