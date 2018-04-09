import pandas as pd
from datetime import date as dt
from datetime import datetime
import getpass
import socket
from ase import Atoms
from collections import OrderedDict
from icet.io.logging import logger

logger = logger.getChild('data_container')



class DataContainer:
    """
    Data container class, which serves for storing
    all information concerned with Monte Carlo
    simulations performed with mchammer.
    """

    def __init__(self, ensemble):
        """
        Initialize a DataContainer object.

        Attributes
        ----------
        ensemble : BaseEnsemble object
        """

        self._observables = OrderedDict()
        self._parameters = OrderedDict()
        self._metadata = OrderedDict()
        self._data = pd.DataFrame()

        # self._add_metadata('ensemble-name', ensemble.name)
        self._add_metadata('date-created', datetime.now())
        self._add_metadata('username', getpass.getuser())
        self._add_metadata('hostname', socket.gethostname())

    def add_structure(self, atoms):
        """
        Add a reference atomic structure to the data container.

        Parameters
        ----------
        atoms : ASE Atoms object
        """
        assert isinstance(atoms, Atoms), \
            'Structure must be provided as ASE Atoms object'
        self.structure = atoms

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
            'Unknown observable type'
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
        self._parameters[tag] = value

    def _add_metadata(self, tag, value):
        """
        Add metadata to the data container.

        Parameters
        ----------
        tag : str
            name of the metadata field
        value : str, date, int
            value of the metadata
        """
        assert isinstance(tag, str), \
            'Metadata tag has wrong type (str)'
        assert isinstance(value, (str, dt, int)), \
            'Unknown metadata type {}'.format(value)
        if tag in self._metadata:
            logger.warning('Cannot overwrite existing metadata ({})'
                           .format(tag))
            return
        self._metadata[tag] = value

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

        self._row_data = OrderedDict()
        self._row_data['mctrial'] = mctrial
        self._row_data.update(record)

        self._data = self._data.append(self._row_data,
                                       ignore_index=True)

    def get_data(self, tags):
        """
        Returns a list of lists with the current
        data stored in the Pandas data frame.

        Parameters
        ----------
        tags : list of str
            list with tags of the required properties.
        """
        data_list = []
        for row in range(len(self._data)):
            data_row = []
            for tag in tags:
                assert tag in self._data, \
                    'observable/parameter is not part of DataContainer'
                data_row.append(self._data.loc[row, tag])
            data_list.append(data_row)
        return data_list

    @property
    def parameters(self):
        """list : ensemble parameters"""
        return self._parameters.copy()

    @property
    def observables(self):
        """list of BaseObserver objects : observables"""
        return self._observables

    @property
    def metadata(self):
        """list : metadata associated with data container"""
        return self._metadata

    def reset(self):
        """Reset data frame of data container."""
        self._data = pd.DataFrame()

    def __len__(self):
        """number of rows in data frame"""
        return len(self._data)

    def get_average(self, tag: str):
        """
        Return average over MonteCarlo steps of parameter
        and/or observable.

        Parameters
        ----------
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
