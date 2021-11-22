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


Cluster.get_distances = get_distances
