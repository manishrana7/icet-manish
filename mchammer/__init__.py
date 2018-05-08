# -*- coding: utf-8 -*-
from __future__ import print_function, division

from .ensembles.base_ensemble import BaseEnsemble
from .observers.base_observer import BaseObserver
from .data_container import DataContainer

'''
mchammer module of icet.
'''

__project__ = 'icet-mchammer'
__description__ = 'A Pythonic approach to cluster expansions'
__authors__ = ['Mattias Ångqvist',
               'William Armando Muñoz',
               'Thomas Holm Rod',
               'Paul Erhart']
__copyright__ = '2018'
__license__ = 'MIT'
__credits__ = ['Mattias Ångqvist',
               'William Armando Muñoz',
               'Thomas Holm Rod',
               'Paul Erhart']
__version__ = '0.1'
__all__ = ['DataContainer',
           'BaseEnsemble',
           'BaseObserver']
__maintainer__ = 'The icet developers team'
__maintainer_email__ = 'icet@materialsmodeling.org'
__status__ = 'alpha-version'
__url__ = 'http://icet.materialsmodeling.org/'
