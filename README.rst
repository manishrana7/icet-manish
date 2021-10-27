icet
====

**icet** is a tool for the construction and sampling of alloy cluster expansions.
A detailed description of the functionality provided as well as an extensive tutorial can be found in the `user guide <https://icet.materialsmodeling.org/>`_.

**icet** is written in Python, which allows easy integration with countless first-principles codes and analysis tools accessible from Python, and allows for a simple and intuitive user interface.
All computationally demanding parts are, however, written in C++ providing performance while maintaining portability.
The following snippet illustrates how one can train a cluster expansion:

.. code-block:: python

   cs = ClusterSpace(primitive_cell, cutoffs, species)
   sc = StructureContainer(cs)
   for structure in training_structures:
       sc.add_structure(structure)
   opt = Optimizer(sc.get_fit_data())
   opt.train()
   ce = ClusterExpansion(cs, opt.parameters)

Afterwards the cluster expansion can be used, e.g., for finding ground state structures or sampled via Monte Carlo simulations.


Installation
------------

**icet** can be installed using `pip <https://pypi.org/project/icet/>`_::

    pip3 install icet --user

or via `conda <https://anaconda.org/conda-forge/icet>`_::

    conda install -c conda-forge icet

Installation via `pip` requires a C++11 compliant compiler.
Please consult the `installation section of the user guide <https://icet.materialsmodeling.org/installation.html>`_ for details.

**icet** is based on Python3 and invokes functionality from other Python libraries including
`ase <https://wiki.fysik.dtu.dk/ase>`_,
`pandas <https://pandas.pydata.org/>`_,
`numpy <http://www.numpy.org/>`_,
`scipy <https://www.scipy.org/>`_,
`spglib <https://atztogo.github.io/spglib/>`_, and
`trainstation <https://trainstation.materialsmodeling.org/>`_.


Credits
-------

**icet** has been developed at the Department of Physics at Chalmers University of Technology (Gothenburg, Sweden), in collaboration with the Data Analysis group at the `Data Management and Software Center of the European Spallation Source <https://europeanspallationsource.se/data-management-software#data-analysis-modelling>`_ (Copenhagen, Denmark).

When using **icet** in your research please cite

| M. Ångqvist, W. A. Muñoz, J. M. Rahm, E. Fransson, C. Durniak, P. Rozyczko, T. H. Rod, and P. Erhart
| *ICET – A Python Library for Constructing and Sampling Alloy Cluster Expansions*
| Adv. Theory. Sim., 1900015 (2019)
| `doi: 10.1002/adts.201900015 <https://doi.org/10.1002/adts.201900015>`_

Also consult the `credits <https://icet.materialsmodeling.org/credits>`_ page of the documentation for additional references.

For questions and help please use the `icet discussion forum on matsci.org <https://matsci.org/icet>`_.
**icet** and its development are hosted on `gitlab <https://gitlab.com/materials-modeling/icet>`_.
