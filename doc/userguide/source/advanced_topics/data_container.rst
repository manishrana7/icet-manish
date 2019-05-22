.. _example_data_container:
.. highlight:: python
.. index::
   single: Examples; Data container

Data container
==============


Getting data
------------
The raw observable data as a function of mctrial steps can be obtained via the ``get_data``
function. This function allows for selecting a start and stop mctrial step. This is particular
useful when e.g. wanting to discard the equilibration part of the simulations. 

.. code-block:: python

    energy = dc.get_data('potential', start=5000)

The ``get_data`` also allows for pairing multiple observerables for the mctrial steps they all exist.

.. code-block:: python

    mctrial, energy, sro = dc.get_data('potential', start=5000)


Trajectory
----------
The trajectory from the MC simulation can be obtained via

.. code-block:: python

    traj = dc.get_trajectory()

The trajectory can also be obtained from the ``get_data`` function, which also allows for pairing the snapshots in the trajectory with observables in the data container.

.. code-block:: python

    E_mix, traj = dc.get_trajectory('potential', 'trajectory')


Updating data container
-----------------------
If you forgot to add an observer when running the MC simulation it can be applied after the fact,
given that the trajectory is stored in the datacontainer. This is done via the ``apply_observer``
function. The datacontainer can then be rewritten to file using the ``write`` function.

.. code-block:: python

    obs = ClusterExpansionObserver(ce, tag='new_obs')
    dc = DataContainer.read('my_dc.dc')
    dc.apply_observer(obs)
    new_obs_data = dc.get_data('')
    dc.write('my_dc.dc')


Data analysis
-------------
The datacontainer also allows for more detailed analysis. The ``analyze_data`` function computes
the average, std, correlation-length and 95% error estimate of the average for a given observable.

.. code-block:: python

    summary = dc.analyze_data('potential')
    print(summary)

The correlation length, :math:`s`, is estimated from the autocorrelation function (ACF).
When the ACF have decayed below :math:`\mathrm{e^{-2}}` the configurations are said to be
uncorrelated and this thus provides an estimate of the correlation length.

.. figure::
    autocorrelation.svg

The error estimate (https://en.wikipedia.org/wiki/Standard_error) of the average can be calculated via

.. math::
    \mathrm{error} = \frac{t \sigma }{\sqrt{Ns}}

where :math:`\sigma` is the standard deviation, :math:`N` the number of samples, :math:`s` the correlation length and :math:`t` is the t-factor (https://en.wikipedia.org/wiki/Student%27s_t-distribution>) which should be adjusted depending on the confidence interval desired. 


Obtaining the autocorrelation function directly or carrying out error estimates can be done via functionality provided in ``mchammer/data_analysis.py``.