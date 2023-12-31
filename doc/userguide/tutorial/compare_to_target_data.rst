.. _tutorial_compare_to_target_data:
.. highlight:: python
.. index::
   single: Tutorial; Compare with target data

Comparison with target data
===========================

In this step the performance of the cluster expansion constructed in the
:ref:`previous step <tutorial_construct_cluster_expansion>` will be tested
against the target data. After loading the CE from file, we loop over all
configurations in the database of reference structure and compile the
concentration as well as the target and predicted mixing energies into a
dictionary. Note that the latter of the three values is calculated by calling
the :func:`predict() <icet.ClusterExpansion.predict>` method of the
:class:`ClusterExpansion <icet.ClusterExpansion>` object with the :class:`ASE
Atoms <ase.Atoms>` object that represents the present structure as input
argument.

.. literalinclude:: ../../../examples/tutorial/2_compare_to_target_data.py
   :start-after: # step 1
   :end-before: # step 2

Once this step has completed the predicted and target mixing energies are
plotted as functions of the concentration.

.. literalinclude:: ../../../examples/tutorial/2_compare_to_target_data.py
   :start-after: # step 2

The figure generated by this diagram is shown below.

.. figure:: _static/mixing_energy_comparison.png

  Predicted (orange crosses) and target (blue circles) mixing energies versus
  concentration for the structures used in the construction of the
  cluster expansion.


Source code
-----------

.. container:: toggle

    .. container:: header

       The complete source code is available in
       ``examples/tutorial/2_compare_to_target_data.py``

    .. literalinclude:: ../../../examples/tutorial/2_compare_to_target_data.py
