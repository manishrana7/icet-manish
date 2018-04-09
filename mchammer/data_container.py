import pandas as pd
from datetime import date as dt
from datetime import datetime
import getpass
import socket
from ase import Atoms
from collections import OrderedDict
from icet.io.logging import logger
from mchammer.observers.base_observer import BaseObserver

logger = logger.getChild('data_container')


class DataContainer:
    """
    Data container class, which serves for storing
    all information concerned with Monte Carlo
    simulations performed with mchammer.
    """

    def __init__(self, atoms, ensemble_name, random_seed):
        """
        Initialize a DataContainer object.

        Paremeter
        ---------
        atoms : ASE Atoms object
            Add a reference atomic structure to the data container.

        name_ensemble : str
            BaseEnsemble name

        random_seed : int
            BaseEnsemble random seed

        Attributes
        ----------
        observables : dict
            dictionary of tag-type pair of added observables.

        parameters : dict
            dictionary of tag-value pair of added parameters.
        
        metadata : dict
            dictionary of tag-value pair of added metadata.

        data : Pandas data frame object
            Runtime data collected during the Monte Carlo simulation. 
        """

        if atoms is not None:
            assert isinstance(atoms, Atoms), \
            'Structure must be provided as ASE Atoms object'
            self.structure = atoms.copy()

        self._observables = OrderedDict()
        self._parameters = OrderedDict()
        self._metadata = OrderedDict()
        self._data = pd.DataFrame()

        self.add_parameter('random-seed', random_seed)

        self._metadata['ensemble-name'] =  ensemble_name
        self._metadata['date-created'] = datetime.now()
        self._metadata['username'] = getpass.getuser()
        self._metadata['hostname'] = socket.gethostname()

    def add_observable(self, tag: str, obs_type):
        """
        Add an observable to the list of observables.

        Parameters
        ----------
        tag : str
            name of observable used as header for the Pandas data frame
        obs_type : type
            type of observable parameter
        """
        assert isinstance(tag, str), \
            'Observable tag has wrong type (str)'
        assert isinstance(obs_type, (type, list)), \
            'Unknown observable type: {}'.format(obs_type)
        self._observables[tag] = obs_type

    def add_parameter(self, tag, value):
        """
        Add parameter required for initial setting
        of the associated ensemble.

        Parameters :
        ----------
        tag : str
            name of the parameter
        value : int, float, list of int or float
            parameter value
        """
        assert isinstance(tag, str), \
            'Parameter tag has the wrong type (str).'
        assert isinstance(value, (int, float, list)), \
            'Unknown parameter type: {}'.format(type(value))

        if isinstance(value, list):
            self._parameters[tag] = value[:]
        else:
            self._parameters[tag] = value


    def append(self, mctrial: int, record: dict):
        """
        Append runtime data to buffer dict.

        Parameters
        ----------
        mctrial : int
            current Monte Carlo step
        record : dict
            dictionary of tag-value pairs representing observations

        Todo
        ----
        Might be a quite expensive way to add data
        to data frame.
        """
        assert isinstance(mctrial, int), \
            'Monte Carlo step has wrong type (int)'
        assert isinstance(record, dict), \
            'Input record has wrong type (dict)'

        row_data = OrderedDict()
        row_data['mctrial'] = mctrial
        row_data.update(record)

        self._data = self._data.append(row_data,
                                       ignore_index=True)

    def get_data(self, tags, fill_missing=False):
        """
        Returns a list of lists with the current
        data stored in the Pandas data frame.

        Parameters
        ----------
        tags : list of str
            list with tags of the required properties.
        """
        import math
        data_list = []
        for row in range(len(self._data)):
            data_row = []
            for tag in tags:
                assert tag in self._data, \
                    'observable is not part of DataContainer: {}'.format(tag)
                data_elem = self._data.loc[row, tag] 
                if math.isnan(data_elem):
                    if fill_missing and row>0:
                        data_elem = self._data.loc[row-1, tag]
                    else:
                        data_elem = None  
                data_row.append(data_elem)
            data_list.append(data_row)
        return data_list

    @property
    def parameters(self):
        """list : simulation parameters"""
        return self._parameters.copy()

    @property
    def observables(self):
        """dict of BaseObservers objects : observables"""
        return self._observables

    @property
    def metadata(self):
        """dict : metadata associated with data container"""
        return self._metadata

    def reset(self):
        """Reset data frame of data container."""
        self._data = pd.DataFrame()

    def __len__(self):
        """number of rows in data frame"""
        return len(self._data)

    def get_average(self, start, stop, tag: str, ):
        """
        Return average of an observable over an interval of trial steps.

        Parameters
        ----------
        start : int
            lower limit of trial step interval
        stop : int
            upper limit of trial step interval
        tag : str
            tag of field over which to average
        """
        pass

    def read(self, infile):
        """
        Read DataContainer object from file.

        Parameters
        ----------
        infile : str or FileObj
            file from which to read
        """
        pass

    def write(self, outfile):
        """
        Write DataContainer object to file.

        Parameters
        ----------
        outfile : str or FileObj
            file to which to write
        """
        pass
