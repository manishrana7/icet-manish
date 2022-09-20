.. _training_of_ce:
.. highlight:: python
.. index::
   single: Examples; Training of CE

Training of cluster expansions
==============================
Here, we will explore some more advanced concepts regarding the training of a cluster expansion.
Specifically we will consider

1. Selecting cutoffs
2. Selecting an optimization algorithm and hyper-parameters
3. Ensemble training


Note here that steps 1 and 2 are to some degree dependent on each other, i.e., the optimal cutoffs may depend on what optimization algorithm you use.
But here, for simplicity, we will consider them in succession.

We consider the Ag-Pd system with the same training structures as used in the :ref:`tutorial <tutorial>`.

Note also when doing cross-validation with small datasets one always ends up with some amount of noise, e.g., in the resulting RMSEs.


.. toctree::
   :maxdepth: 2
   :caption: Contents

   training_cutoffs_selection
   training_hyper_parameter_scans
   training_ensemble_of_models
