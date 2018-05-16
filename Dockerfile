# Base image
FROM ubuntu:18.04

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
                        doxygen \
                        pyflakes \
                        pep8 \
                        graphviz

# Set up some Python3 packages via pip
RUN pip3 install --upgrade pip
RUN pip3 install --upgrade setuptools
RUN pip3 install --upgrade \
                 ase \
                 breathe \
                 cloud_sptheme \
                 coverage \
                 flake8 \
                 pandas \
                 scikit-learn \
                 spglib \
		 sphinx \
                 sphinx-rtd-theme \
                 sphinx_autodoc_typehints \
                 sphinx_sitemap \
                 sphinxcontrib-bibtex
