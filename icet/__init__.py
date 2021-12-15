# -*- coding: utf-8 -*-
"""
Main module of the icet package.
"""

from .core.cluster_space import ClusterSpace
from .core.cluster_expansion import ClusterExpansion
from .core.structure_container import StructureContainer
from .core.variable_basis_cluster_space import VariableBasisClusterSpace
from .core.variable_basis_cluster_expansion import VariableBasisClusterExpansion
from .core.variable_basis_structure_container import VariableBasisStructureContainer

__project__ = 'icet'
__description__ = 'A Pythonic approach to cluster expansions'
__copyright__ = '2021'
__license__ = 'Mozilla Public License 2.0 (MPL 2.0)'
__version__ = '1.4'
__maintainer__ = 'The icet developers team'
__email__ = 'icet@materialsmodeling.org'
__status__ = 'Stable'
__url__ = 'http://icet.materialsmodeling.org/'

__all__ = ['ClusterSpace',
           'ClusterExpansion',
           'StructureContainer',
           'VariableBasisClusterSpace',
           'VariableBasisClusterExpansion',
           'VariableBasisStructureContainer']
