"""
This module provides the Cluster class.
"""
from typing import List

from _icet import Cluster

import numpy as np

__all__ = ['Cluster']


def get_distances(obj) -> List[float]:
    """
    Calculates the distances between all pairs of points in the cluster.
    """
    positions = np.array(obj.positions)
    distances = []
    for i in range(1, positions.shape[0]):
        for j in range(i):
            distances.append(np.linalg.norm(positions[i, :] - positions[j, :]))
    return distances


def _get_string_representation(obj) -> str:
    """
    String representation of the cluster that provides an overview
    (order, radius, positions, distances).
    """
    def repr_lattice_site(lattice_site, position, header=False):
        formats = {'unitcell_index': '{:3}',
                   'unitcell_offset': '{:10.5f} {:10.5f} {:10.5f}',
                   'position': '{:10.5f} {:10.5f} {:10.5f}'}
        data = {'unitcell_index': [lattice_site.index],
                'unitcell_offset': lattice_site.unitcell_offset,
                'position': position}
        s = []
        for name, value in data.items():
            str_repr = formats[name].format(*value)
            n = max(len(name), len(str_repr))
            if header:
                s += ['{s:^{n}}'.format(s=name, n=n)]
            else:
                s += ['{s:^{n}}'.format(s=str_repr, n=n)]
        return ' | '.join(s)

    # basic information
    # (use lattice site to obtain maximum line length)
    prototype_lattice_site = obj.lattice_sites[-1]
    distances_string = ' '.join([f'{d:3.5f}' for d in obj.get_distances()])
    width = max(len(distances_string), len(repr_lattice_site(prototype_lattice_site, [10, 10, 10])))
    s = []  # type: List
    s += ['{s:=^{n}}'.format(s=' Cluster ', n=width)]
    s += [f' {"order":12} : {obj.order}']
    s += [f' {"radius":12} : {obj.radius:.5f}']
    s += [f' {"distances":12} : {distances_string}']

    # table header
    s += [''.center(width, '-')]
    s += [repr_lattice_site(prototype_lattice_site, [10, 10, 10], header=True)]
    s += [''.center(width, '-')]

    # table body
    for lattice_site, position in zip(obj.lattice_sites, obj.positions):
        s += [repr_lattice_site(lattice_site, position)]

    s += [''.center(width, '=')]

    return '\n'.join(s)


# Attach functions to C++ Cluster objects such that they can be used
# even if the Cluster is initialized on the C++ side.
Cluster.get_distances = get_distances
Cluster.__str__ = _get_string_representation
