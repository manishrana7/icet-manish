"""
This module provides the Orbit class.
"""

from _icet import Orbit

__all__ = ['Orbit']


@property
def distances(self):
    """Returns the distances between all sites in the representative cluster"""
    return self.representative_cluster.distances


@property
def sites(self):
    """Returns the site indices of all sites in the representative cluster"""
    return [site.index for site in self.representative_cluster.lattice_sites]


@property
def site_offsets(self):
    """Returns the unitcell offsets of all sites in the representative cluster"""
    return [site.unitcell_offset.astype(int)
            for site in self.representative_cluster.lattice_sites]


@property
def positions(self):
    """Returns the position of all sites in the representative cluster"""
    return self.representative_cluster.positions


@property
def all_distances(self):
    """Returns the distances between all sites in all clusters"""
    return [cluster.distances for cluster in self.clusters]


@property
def all_sites(self):
    """Returns the site indices of all sites in all clusters"""
    return [[site.index for site in cluster.lattice_sites]
            for cluster in self.clusters]


@property
def all_site_offsets(self):
    """Returns the unitcell offsets of all sites in all clusters"""
    return [[site.unitcell_offset.astype(int) for site in cluster.lattice_sites]
            for cluster in self.clusters]


@property
def all_positions(self):
    """Returns the position of all sites in all clusters"""
    return [cluster.positions for cluster in self.clusters]


Orbit.distances = distances
Orbit.sites = sites
Orbit.site_offsets = site_offsets
Orbit.positions = positions

Orbit.all_distances = all_distances
Orbit.all_sites = all_sites
Orbit.all_site_offsets = all_site_offsets
Orbit.all_positions = all_positions
