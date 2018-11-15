.. index:: Installation

Installation
************

Requirements
============

:program:`icet` requires a C++11 compiler, eigen3, Boost, Python3 and depends on the
`numpy <http://www.numpy.org/>`_,
`scipy <https://www.scipy.org/>`_,
`spglib <https://atztogo.github.io/spglib/>`_, and
`ase <https://wiki.fysik.dtu.dk/ase>`_ packages.
In order to be able to use more advanced optimization techniques
one also requires `scitkit-learn <http://scikit-learn.org/>`_.

Setup
=====

Setting up :program:`icet` involves

  1. Installing basic libraries
  2. Installing Python libraries
  3. Acquiring the sources
  4. Compiling the :program:`icet` C++ module
  5. Adding :program:`icet` to the Python path

These steps are detailed in the sections below.
Afterwards the correctness of the installation :ref:`ought to be tested <testing>`.


Installing basic libraries
--------------------------

Next one must install a number of basic libraries that are needed for compilation of the C++ components.
This part is system dependent as described in the following sections.

Debian, Ubuntu etc.
^^^^^^^^^^^^^^^^^^^

To install :program:`icet` on `Debian <https://en.wikipedia.org/wiki/Debian>`_ and Debian based distros
including e.g., `Ubuntu <https://en.wikipedia.org/wiki/Ubuntu>`_ and `Linux Mint <https://en.wikipedia.org/wiki/Linux_Mint>`_
the basic libraries can be installed using the commands below.
This approach has been tested using the latest `Ubuntu Docker image <https://hub.docker.com/_/ubuntu/>`_::

    # update repositories
    sudo apt update
    # upgrade packages
    sudo apt upgrade -qy
    # install packages
    sudo apt install -qy \
        cmake \
        build-essential \
        libboost-dev \
        libeigen3-dev \
        python3-dev
    # install pip (to be used for installing Python libraries below)
    sudo apt install -qy \
        python3-pip


.. warning::

    The latest versions of Debian (Stretch)
    and `Mint <https://hub.docker.com/r/vcatechnology/linux-mint/>`_
    include only Python 3.5, which is `currently not supported <https://gitlab.com/materials-modeling/icet/issues/269>`_.

Fedora
^^^^^^

To install :program:`icet` on `Fedora <https://getfedora.org/>`_.
the basic libraries can be installed using the commands below.
This approach has been tested using the latest `Fedora Docker image <https://hub.docker.com/_/fedora/>`_::

    # update repositories
    sudo yum update -y
    # upgrade packages
    sudo yum upgrade -y
    # install packages
    sudo yum install -y \
        cmake \
        make \
        gcc-c++ \
        kernel-devel \
        boost-devel \
        eigen3 \
        python3-devel
    # install pip (to be used for installing Python libraries below)
    sudo yum install -y \
        python3-pip

Mac OS
^^^^^^

.. warning::

    :program:`icet` can be compiled on Mac OS but 4 out of 386 tests fail.
    Until this behavior is `resolved <https://gitlab.com/materials-modeling/icet/issues/270>`_
    it is not recommended to run :program:`icet` on Mac OS.

.. comment

    To install :program:`icet` on Mac OS one should employ a package manager such as
    `Homebrew <https://en.wikipedia.org/wiki/Homebrew_(package_management_software)>`_.
    One also requires Apple's compilers, which can be obtained as part of `Xcode <https://en.wikipedia.org/wiki/Xcode>`_.
    Provided the compilers have been installed and using Homebrew one can install the necessary packages as follows::

        brew install \
            cmake \
            make \
            gcc \
            boost \
            eigen


Installing Python libraries
---------------------------

:program:`icet` relies on several Python libraries.
The two most basic ones are `numpy <http://www.numpy.org/>`_ and `scipy <https://www.scipy.org/>`_.
At least the former is often already installed as part of the standard Python environment.
If you need to install any of these packages yourself
we recommend using `pip <https://en.wikipedia.org/wiki/Pip_(package_manager)>`_::

    python3 -m pip install \
        --user --upgrade \
        numpy \
        scipy

Here, we use the ``--user`` option, which implies that the libraries
are installed in the home directory of the current user.
:program:`icet` also invokes several more specialized packages, which can be readily installed as follows::

    python3 -m pip install \
        --user --upgrade \
        ase \
        pandas \
        scikit-learn \
        spglib


Acquiring the sources
---------------------

To begin with, one has to download the :program:`icet` sources.
Most commonly this is achieved by cloning the :program:`icet` repo::

    git clone https://gitlab.com/materials-modeling/icet.git


Compiling the :program:`icet` C++ module
----------------------------------------

The compilation of :program:`icet` is configured using `CMake <https://en.wikipedia.org/wiki/CMake>`_.
If the following snippet is run in the :program:`icet` home directory
the :program:`icet` C++ library will be built in the ``build`` directory::

    mkdir build
    cd build
    cmake ..
    make -j4
    cd ..

Here, the ``-j4`` option instructs ``make`` to use four cores in parallel (if available), which commonly speeds up the build process.


Adding :program:`icet` to the Python path
-----------------------------------------

Now :program:`icet` must be added to the ``PYTHONPATH`` environment variable.
To this end, when using the `Bash shell <https://en.wikipedia.org/wiki/Bash_(Unix_shell)>`_
or similar (bash, ksh) the following command should be added to the ``.bashrc`` file (or equivalent) in the home directory::

    export PYTHONPATH=${PYTHONPATH}:<ICET_PATH>/
    export PYTHONPATH=${PYTHONPATH}:<ICET_PATH>/build/src/

Here, ``ICET_PATH`` must be replaced with the path to the :program:`icet` root directory.
If you are using `C shell <https://en.wikipedia.org/wiki/C_shell>`_ (csh, tcsh) the equivalent line reads::

    setenv PYTHONPATH ${PYTHONPATH}:<ICET_PATH>/
    setenv PYTHONPATH ${PYTHONPATH}:<ICET_PATH>/build/src/

.. _testing:

Testing
=======

Finally, it is strongly recommended to run the test suite in order to ensure that all parts of :program:`icet` function properly.
To this end, the following command should be executed at the command line::

    python3 tests/main.py

Running the test suite will commonly take several minutes.
