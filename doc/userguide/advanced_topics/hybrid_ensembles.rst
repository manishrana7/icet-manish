.. _advanced_topics_parallel_monte_carlo_simulations:
.. highlight:: python
.. index::
   single: Advanced topics; Hybrid ensembles

Hybrid ensembles
================

In systems with multiple sublattices, it can sometimes be advantageous to use different thermodynamical ensembles for different sublattices in a single Monte Carlo simulation.
This can for example be relevant in systems in which one sublattice does not exchange atoms with the environment (closed system) while the other does (open system).
Metallic alloys exposed to hydrogen are examples where this is often the case (for an example, see [RahLofFra21]_).
This tutorial demonstrates how sublattice-specific ensembles can be used in :program:`mchammer` using the :class:`HybridEnsemble <mchammer.ensembles.HybridEnsemble>` class.

As in any Monte Carlo simulation with :program:`mchammer`, the first steps are to define a simulation cell and to construct a :class:`ClusterExpansionCalculator <mchammer.calculators.ClusterExpansionCalculator>` object.
To this end, we first construct a toy :class:`ClusterExpansion <icet.ClusterExpansion>` for a system with two sublattices, one occupied by Pd/Au and one occupied by H and vacancies.

.. literalinclude:: ../../../examples/advanced_topics/hybrid_ensembles.py
   :start-after: # step 1
   :end-before: # step 2

We then define the parameters that will enter our Monte Carlo simulation.
Here, we will run a simulation in which the Pd-Au sublattice is sampled in the canonical ensemble and the H-vacancy sublattice in the semi-grand canonical ensemble.
The ensembles are specified via a list of dictionaries, which define the parameters specific to each ensemble.
Moreover, since the concentrations on the Pd-Au sublattice are fixed once a starting configuration is defined, we must also create a supercell with the desired concentration
(note that the concentrations on the H-vacancy sublattice will change during the simulation since concentrations are not conserved in the semi-grand canonical ensemble, hence the choice of starting concentrations on the H-vacancy sublattice is unimportant).

.. literalinclude:: ../../../examples/advanced_topics/hybrid_ensembles.py
   :start-after: # step 2
   :end-before: # step 3

Finally, we define our :class:`HybridEnsemble <mchammer.ensembles.HybridEnsemble>` and commence the simulation.

.. literalinclude:: ../../../examples/advanced_topics/hybrid_ensembles.py
   :start-after: # step 3

The :class:`HybridEnsemble <mchammer.ensembles.HybridEnsemble>` can also be used together with the variance-constrained semi-grand canonical (VCGSC) ensemble.
The following is a valid specification for a simulation with the canonical ensemble on the Pd/Au sublattice and VCSGC on the H/vacancy sublattice:

.. code-block:: python

   ensemble_specs = [{'ensemble': 'canonical', 'sublattice_index': 0},
                     {'ensemble': 'vcsgc', 'sublattice_index': 1,
                      'phis': {'H': phiH}, 'kappa': 200}]


Source code
-----------

.. container:: toggle

    .. container:: header

       The complete source code is available in
       ``examples/advanced_topics/hybrid_ensembles.py``

    .. literalinclude:: ../../../examples/advanced_topics/hybrid_ensembles.py
