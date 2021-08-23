.. _advanced_topics_fitting_with_constraints:
.. highlight:: python
.. index::
   single: Advanced topics; Fitting with constraints

Fitting with constraints
========================

To construct cluster expansions that yield physically reasonable results,
it is sometimes advantageous to enforce certain constraints. For example,
if we fit a cluster expansion for the mixing energy of a binary system,
we know that, by definition, the energy should be zero for the pure phases.
Another example is low-symmetry systems that contain a large number of
symmetrically inequivalent orbits, but in which those orbits nevertheless
can be expected to have similar interactions, such as pair clusters in the
sixth and seventh atomic layer of a surface slab. Such insight can be
encoded into the cluster expansion by imposing constraints during fitting.
This can be done in several ways, and this tutorial demonstrates a few
different approaches.




Source code
-----------

.. container:: toggle

    .. container:: header

       The complete source code is available in
       ``examples/sqs_generation.py``

    .. literalinclude:: ../../../../examples/advanced_topics/sqs_generation.py
