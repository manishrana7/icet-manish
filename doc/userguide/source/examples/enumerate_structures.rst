.. _example_enumerate_structures:
.. highlight:: python
.. index::
   single: Tutorial; Structure enumeration

Structure enumeration
=====================

The purpose of this example is to demonstrate how to use functionalities, which are built into :program:`iceT`, to *enumerate structures*. In the present context, this *structural enumeration* means *the generation of all inequivalent structures derived from a primitive structure up to a certain size*.

Import modules
--------------

The :func:`enumerate_structures <icetdev.enumeration.enumerate_structures>` function, which is used for performing the `structural enumeration` needs to be imported together with some additional functions from `ASE <https://wiki.fysik.dtu.dk/ase>`_.

.. literalinclude:: ../../../../examples/enumerate_structures.py
   :start-after: # Import modules
   :end-before: # Generate all binary

.. _generate-binary-au-pd-structures:
Generate binary Au/Pd structures
--------------------------------

Before being able to perform the *structural enumeration*, it is first necessary to generate a primitive structure. In this case, an Au fcc :class:`ase.Atoms` object is created using the :func:`ase.build.bulk` function. Then a database ``AuPd-fcc.db`` is initiated, in which the enumerated structures will be stored. All possible binary Au/Pd structures, with between 1 and 6 atoms/cell are then generated and stored in this database.

.. literalinclude:: ../../../../examples/enumerate_structures.py
   :start-after: # and save them
   :end-before: # Enumerate all palladium

Generate PdH structures, with vacancies
---------------------------------------

The :ref:`steps above <_generate-binary-au-pd-structures>` are now repeated to enumerate all palladium hydride structures with up to 4 primtive cells, which contain up to 4 Pd atoms and between 0 and 4 H atoms. In addition, vacancies, represented by vanadium, are also included, which results in ternary systems. The structures thus obtained are stored in database named ``PdHVac-fcc.db``.

.. literalinclude:: ../../../../examples/enumerate_structures.py
   :start-after: # either a hydrogen

Source code
-----------

.. container:: toggle

    .. container:: header

       The complete source code is available in
       ``examples/enumerate_structures.py``

    .. literalinclude:: ../../../../examples/enumerate_structures.py
