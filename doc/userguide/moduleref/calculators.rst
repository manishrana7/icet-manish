.. index::
   single: Function reference; calculators
   single: Class reference; calculators

Calculators
===========

Given a force constant model (FCM) the calculator classes described below
provide functionality for readily evaluating the energy and atomic forces.
Using the :class:`ASECalculator` class one can integrate FCMs with integrators
and optimizers supported by the `atomic simulation environment (ASE)
<https://wiki.fysik.dtu.dk/ase/index.html>`_. In this fashion, it is for
example possible to carry out molecular dynamics (MD) simulations. The
calculators also readily enable integration with other Python libraries and
codes that are interfaced with Python.

ASE calculator
--------------

.. index:: ASE calculator
.. automodule:: ase_calculator
   :members:
   :undoc-members:

General FCM calculator
----------------------

.. index:: FCM calculator
.. automodule:: fcm_calculator
   :members:
   :undoc-members:
