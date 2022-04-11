.. _constituent_strain_example:
.. highlight:: python
.. index::
   single: Examples; Constituent strain

Constituent strain calculations
===============================

This example demonstrates how strain can be incorporated in cluster
expansions in :program:`icet` using the constituent strain formalism
originally introduced by Laks *et al.* [LakFerFro92]_.

.. warning ::
  
  Since this module has only been tested for FCC crystals, great caution
  should be exercised when using it, in particular for non-cubic systems.

.. warning ::
  
  Due to the requirement to calculate
  :math:`\Delta E_\text{CS}(\hat{\mathbf{k}}, c)` separately
  (see below), this module is only recommended to advanced users.
  We can only provide very limited support.

In a regular cluster expansion, strain is typically not captured properly.
Consider, for example, a phase separated system with
phases 1 and 2. If the interface is coherent (i.e., the atomic planes
line up without defects), the two phases will have to adopt the same
in-plane lattice parameter. This leads to strain, unless the two phases
happen to have the same equilibrium lattice parameter. A cluster expansion,
however, normally has fairly short cutoffs, and will not be able to describe
this situation, since the strain extends into the phases far from the
interface, well beyond the cutoffs of the cluster expansion. Furthermore,
the strain energy is dependent not only on the composition of the phases,
but also the orientation of the interface, since the stiffness tensors are
not isotropic.

To overcome this limitation of cluster expansions, the constituent strain
formalism adds a term to the cluster expansion that is a sum in reciprocal
space,

.. math ::

   \text{strain energy} = \sum_\mathbf{k} \Delta E_\text{CS}
      (\hat{\mathbf{k}}, c) F(\mathbf{k}, \sigma)

Here, the sum runs over all :math:`\mathbf{k}` points in the first Brillouin
zone. The term :math:`\Delta E_\text{CS}(\hat{\mathbf{k}}, c)` expresses
the dependency of strain energy on crystallographic direction and
overall concentration :math:`c`, while :math:`F(\mathbf{k}, \sigma)`
quantifies the variation of concentration with periodicity :math:`\mathbf{k}`
for an atomic structure with occupation vector :math:`\sigma`.

:program:`icet` provides functionality for efficiently evaluating the above
sum when fitting cluster expansions and running Monte Carlo simulations.
It does, however, only provide a framework for evaluating
:math:`\Delta E_\text{CS}(\hat{\mathbf{k}}, c)`, not explicit functionality
for its definition for a specific material. For more information, including
an example, see below.

Here, we will demonstrate the usage of the constituent strain module
using the Ag-Cu alloy as an example. Ag-Cu is immiscible and has a large
size mismatch. This means that it will phase separate, and if it does so
coherently there will be significant strain energy constributions. In
practice, it is likely that dislocations will form that make the interface
bewteen Ag and Cu incoherent. If that happens, the strain energy will be
relaxed, but for the purpose of demonstration we will ignore this aspect
here.

Define :math:`\Delta E_\text{CS} (\hat{\mathbf{k}}, c)`
-------------------------------------------------------

The first step in the constituent strain calculation is to define
:math:`\Delta E_\text{CS}(\hat{\mathbf{k}}, c)`. This module does not
provide explicit functionality to calculate this term, but an example for
how it can be done is provided below.

Ideally, we need to be able to evaluate this term for any :math:`\mathbf{k}`
point (or rather for any direction :math:`\hat{\mathbf{k}}` in reciprocal
space) and for any overall concentration :math:`c`. When we later on create 
instances of the :class:`ConstituentStrain <icet.tools.ConstituentStrain>`
class, we will thus supply it with *functions* (and not floats or arrays)
that evaluate :math:`\Delta E_\text{CS} (\hat{\mathbf{k}}, c)`
for given :math:`\hat{\mathbf{k}}` and :math:`c`. There are at least two
approaches for constructing such functions:

#. Deduce them analytically from the stiffness tensors of the two phases, or

#. Do a long series of DFT calculations with strain in different orientations
   and interpolate the results for remaining :math:`\hat{\mathbf{k}}` and
   :math:`c`.

The first of these approaches is described in [LakFerFro92]_. The second
was outlined in [OzoWolZun98]_, and is demonstrated **in the cubic case**
in :download:`this Jupyter notebook
<../../../examples/advanced_topics/constituent_strain/calculate-constituent-strain-parameters.ipynb>`.
The notebook demonstrates a procedure for obtaining the necessary energies
in the form of expansion coefficients for Redlich-Kister polynomials in six
crystallographic directions. To handle arbitrary directions, these expansion
coefficients are interpolated using the functions in
:download:`this module
<../../../examples/advanced_topics/constituent_strain/custom_functions.py>`.
Functions from this module (which is *specific for the cubic case and the
choice of crystallographic directions in the notebook*) will be supplied
as input arguments to :class:`ConstituentStrain <icet.tools.ConstituentStrain>`. 
For further details, please see [RahLofErh22]_, in particular its Supporting Information.


Fit a cluster expansion with strain
-----------------------------------

Fitting of a cluster expansion with strain is not very different from fitting
a regular cluster expansion. The only thing we need to keep in mind is that
the total energy of a configuration :math:`\sigma` consists of two terms,

.. math ::

    E(\sigma) = E_\text{CE}(\sigma) + \Delta E_\text{CS} (\hat{\mathbf{k}}, c).

Here, the first term is described by a regular cluster expansion and the
second term is added afterwards. When fitting the cluster expansion to
energies calculated with DFT (or any other method), we therefore first need
to subtract the second term from this energy, and fit only the remainder
(:math:`E_\text{CE}(\sigma)`).

The following example illustrates the process. For the sake of simplicity,
we will use the EMT calculator in ASE instead of DFT to calculate energies.
We begin by some imports (including our custom functions from the previous 
section), we define a functon that will allow us to calculate
energies, and we calculate the energies of the pure phases (Ag and Cu)
to be able to define mixing energy properly.

.. literalinclude:: ../../../examples/advanced_topics/constituent_strain/1_fit_cluster_expansion.py
   :end-before: # Define a cluster space

We then define a cluster space just as we normally would. We use short cutoffs
here for the sake of demonstration.

.. literalinclude:: ../../../examples/advanced_topics/constituent_strain/1_fit_cluster_expansion.py
   :start-after: # Define a cluster space
   :end-before: # Fill a structure container with data

We then start to fill a structure container with data using enumeration.
Here, we need to calculate strain by initiating
an instance of :class:`ConstituentStrain
<icet.tools.ConstituentStrain>` for each structure. To this end, we use
the functions defined in ``custom_functions``.

.. literalinclude:: ../../../examples/advanced_topics/constituent_strain/1_fit_cluster_expansion.py
   :start-after: # Fill a structure container with data
   :end-before: # Constrain sensing matrix

Finally, we fit the cluster expansion. Here, we use
:func:`get_mixing_energy_constraints <icet.tools.constraints.get_mixing_energy_constraints>`
to reproduce the pure phases exactly.

.. literalinclude:: ../../../examples/advanced_topics/constituent_strain/1_fit_cluster_expansion.py
   :start-after: # then fit a cluster expansion


Run Monte Carlo with strain cluster expansions
----------------------------------------------

The :program:`mchammer` module includes custom calculators and observers
for running Monte Carlo simulations with cluster expansions with strain.
If these are used, the procedure is not very different from Monte Carlo
simulations with a regular cluster expansion. Instead of creating a
:class:`ClusterExpansionCalculator <mchammer.calculators.ClusterExpansionCalculator>`,
we create a
:class:`ConstituentStrainCalculator <mchammer.calculators.ConstituentStrainCalculator>`
with a :class:`ConstituentStrain <icet.tools.ConstituentStrain>`
instance as input. The latter is constructed in the same way as
demonstrated above. If we want to be able to extract the strain energy
separately, we also need to attach an observer to the ensemble
object.

.. note ::
  
  The
  :class:`ConstituentStrainCalculator <mchammer.calculators.ConstituentStrainCalculator>`
  is currently only compatible with ensembles that flip one site at a time.
  Presently it can therefore not be used together with the canonical ensemble.

.. literalinclude:: ../../../examples/advanced_topics/constituent_strain/2_run_monte_carlo.py

.. note ::
  
  Due to the non-locality of the structure factor, Monte Carlo simulations
  with constituent strain are typically significantly more time-consuming
  than if strain is not included.

Source code
-----------

.. container:: toggle

    .. container:: header

       The complete source code is available in
       ``examples/advanced_topics/constituent_strain/1_fit_cluster_expansion.py``

    .. literalinclude:: ../../../examples/advanced_topics/constituent_strain/1_fit_cluster_expansion.py

.. container:: toggle

    .. container:: header

       The complete source code is available in
       ``examples/advanced_topics/constituent_strain/2_run_monte_carlo.py``

    .. literalinclude:: ../../../examples/advanced_topics/constituent_strain/2_run_monte_carlo.py

.. container:: toggle

    .. container:: header

       The source code for interpolation of
       :math:`\Delta E_\text{CS} (\hat{\mathbf{k}}, c)` in the cubic case
       is available in
       ``examples/advanced_topics/constituent_strain/custom_functions.py``

    .. literalinclude:: ../../../examples/advanced_topics/constituent_strain/custom_functions.py
