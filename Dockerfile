# Base image
FROM ubuntu:17.10

# Install required packages
RUN apt-get update -qy
RUN apt-get upgrade -qy
RUN apt-get install -qy cmake \
                        build-essential \
                        libboost-all-dev \
                        libgsl-dev \
                        libeigen3-dev \
                        python3-dev \
                        python3-pip \
                        python3-numpy \
                        python3-scipy \
                        python3-h5py \
                        python3-sphinx \
                        doxygen \
                        pyflakes \
                        pep8 \
                        graphviz

# Set up some Python3 packages via pip
RUN pip3 install --upgrade pip
RUN pip3 install --upgrade setuptools
RUN pip3 install ase \
                 coverage \
                 flake8 \
                 spglib \
                 scikit-learn \
                 sphinx-rtd-theme \
                 sphinxcontrib-bibtex \
                 sphinx_sitemap
