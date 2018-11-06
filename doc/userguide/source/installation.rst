.. index:: Installation

Installation
************

Requirements
============

:program:`icet` requires a C++11 compiler, eigen3, Boost, Python3 and depends on the `numpy
<http://www.numpy.org/>`_, `spglib <https://atztogo.github.io/spglib/>`_, and
`ase <https://wiki.fysik.dtu.dk/ase>`_ packages. If you want to use more
advanced optimization techniques
(`LASSO <http://scikit-learn.org/stable/modules/linear_model.html#lasso>`_,
`Bayesian ridge regression <http://scikit-learn.org/stable/modules/linear_model.html#bayesian-ridge-regression>`_,
`automatic relevance determination regression  <http://scikit-learn.org/stable/modules/linear_model.html#automatic-relevance-determination-ard>`_ etc.) one
also requires `scitkit-learn <http://scikit-learn.org/>`_.

Setup
=====

Setting up :program:`icet` involves

1. acquiring the sources

2. installing a set of basic libraries,

3. installing a set of Python libraries,

4. compiling the C++ part of :program:`icet`,

5. adding :program:`icet` to the Python path, and

6. running the tests


Acquiring the sources
---------------------

First off, one has to download the :program:`icet` sources.
Most commonly this is achieved by cloning the (public) :program:`icet` repo::

    git clone https://gitlab.com/materials-modeling/icet.git


Basic libraries
---------------

Next one must install a number of basic libraries that are needed for compilation of the C++ components.
This part is system dependent as described in the following sections.

Ubuntu (debian-based)
^^^^^^^^^^^^^^^^^^^^^

To install :program:`icet` on debian-based systems we recommend installing the following packages (tested on bionic)::

    # update repositories
    sudo apt-get update
    # install packages
    sudo apt-get install \
        -qy \
        cmake \
        build-essential \
        libboost-all-dev \
        libgsl-dev \
        libeigen3-dev \
        python3-pip \
        python3-numpy \
        python3-scipy \
        python3-h5py


Python libraries
----------------

:program:`icet` relies on several Python libraries, which we recommend installing via pip.
Here, we use the `--user` option, which implies that the libraries are installed in the home directory of the current user.
If :program:`icet` is to be installed system wide the `--user` option must be omitted and the command be executed as root (super user/admin)::

    python3 -m pip install \
        --user --upgrade \
        pip \
        setuptools \
        ase \
        pandas \
        scikit-learn \
        spglib


Compile :program:`icet`
-----------------------

The compilation of :program:`icet` is configured using cmake.
The following snippet is to be run in the :program:`icet` home directory.
The C++ library will be built in the ``build`` directory::

    mkdir build
    cd build
    cmake ..
    make -j4
    cd ..

Here, the ``-j4`` option instructs ``make`` to use four cores in parallel (if available), which commonly speeds up the build process.


Add :program:`icet` to the path
-------------------------------

Now :program:`icet` must be added to the ``PYTHONPATH`` environment variable.
To this end, when using a bash shell or similar the following command should be added to the ``.bashrc`` file (or equivalent) in the home directory::

    export PYTHONPATH=${PYTHONPATH}:$<icet-directory>i/:$<icet-directory>/build/src/

Here, ``<icet-directory>`` must be replaced with the path to the :program:`icet` root directory.


Testing
-------
Finally, it is strongly recommended to run the test suite in order to ensure that all parts of :program:`icet` function properly.
To this end, the following command should be executed at the command line::

    python3 tests/main.py

Running the test suite will commonly take several minutes on an Intel core-i5 or i7 system.
