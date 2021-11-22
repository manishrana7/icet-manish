import pytest

from icet.core.cluster import Cluster
from icet.core.structure import Structure
from icet.core.lattice_site import LatticeSite

from ase.build import bulk
import numpy as np


@pytest.fixture
def lattice_sites(request):
    """Defines lattice sites used to initialize a cluster."""
    order = request.param
    lattice_sites = []
    if order == 'singlet':
        lattice_sites.append(LatticeSite(0, [1, 0, 0]))
    elif order == 'triplet':
        indices = [0, 1, 2]
        offsets = [[0, 0, 0], [1, 0, 0], [-1, 0, 0]]
        for ind, offset in zip(indices, offsets):
            lattice_site = LatticeSite(ind, offset)
            lattice_sites.append(lattice_site)
    return lattice_sites


@pytest.fixture
def structure():
    """Defines structure used to initialize a cluster."""
    prim = bulk('H', a=4.0, crystalstructure='sc').repeat((3, 1, 1))
    return Structure.from_atoms(prim)


@pytest.mark.parametrize('lattice_sites', [
    'singlet',
    'triplet'],
    indirect=['lattice_sites'])
def test_init(lattice_sites, structure):
    """Test initialization of Cluster object."""
    cluster = Cluster(lattice_sites, structure)
    assert isinstance(cluster, Cluster)
    assert len(cluster) == len(lattice_sites)
    assert cluster.order == len(lattice_sites)
    ret_lattice_sites = cluster.lattice_sites
    assert len(ret_lattice_sites) == len(lattice_sites)
    for ret_site, target_site in zip(ret_lattice_sites, lattice_sites):
        assert ret_site == target_site


@pytest.mark.parametrize('lattice_sites, target_radius', [
    ('singlet', 0),
    ('triplet', 8)],
    indirect=['lattice_sites'])
def test_radius(lattice_sites, target_radius, structure):
    """Test radius property."""
    cluster = Cluster(lattice_sites, structure)
    assert abs(cluster.radius - target_radius) < 1e-6


@pytest.mark.parametrize('lattice_sites, target_positions', [
    ('singlet', [[12, 0, 0]]),
    ('triplet', [[0, 0, 0], [16, 0, 0], [-4, 0, 0]])],
    indirect=['lattice_sites'])
def test_positions(lattice_sites, target_positions, structure):
    """Test positions property."""
    cluster = Cluster(lattice_sites, structure)
    ret_positions = cluster.positions
    assert len(ret_positions) == len(target_positions)
    for ret_pos, target_pos in zip(ret_positions, target_positions):
        assert np.allclose(ret_pos, target_pos)


@pytest.mark.parametrize('lattice_sites, target_distances', [
    ('singlet', []),
    ('triplet', [16, 4, 20])],
    indirect=['lattice_sites'])
def test_get_distances(lattice_sites, target_distances, structure):
    """Test get_distances function."""
    cluster = Cluster(lattice_sites, structure)
    assert np.allclose(cluster.get_distances(), target_distances)
