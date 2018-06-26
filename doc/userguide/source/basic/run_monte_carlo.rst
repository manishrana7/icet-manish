.. _tutorial_monte_carlo_simulations:
.. highlight:: python
.. index::
   single: Tutorial; Monte Carlo simulations

Monte Carlo simulations
=======================

We are now in a position to carry out a series of Monte Carlo (MC)
simulations to sample the cluster expansion model that was constructed
and validated in the previous steps. To set up the simulation we set
up a supercell and then create a calculator by combining our cluster
expansion model with the supercell.

.. literalinclude:: ../../../../tutorial/basic/6_run_monte_carlo.py
   :end-before: # step 2

In this example the sampling will be carried out in the semi-grand
canonical (SGC) ensemble, which requires a temperature and a set of
chemical potentials as input. Accordingly we set up a loop over
different temperature and chemical potential values. In the body of
the two nested loops we instantiate a SGC ensemble object using the
calculator configured before and then run a MC simulation for a number
of trial steps.

.. literalinclude:: ../../../../tutorial/basic/6_run_monte_carlo.py
   :start-after: # step 2

Here, the results of the simulation are stored in a data container,
which is written to a file named `sgc.dc`.  In the next step these
data will be analyzed to generate e.g., a map of the chemical
potential difference vs composition.


Source code
-----------

.. container:: toggle

    .. container:: header

       The complete source code is available in
       ``tutorial/basic/6_run_monte_carlo.py``

    .. literalinclude:: ../../../../tutorial/basic/6_run_monte_carlo.py
