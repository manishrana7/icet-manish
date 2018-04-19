import tarfile
import tempfile
import getpass
import socket
import pickle
from datetime import datetime
from collections import OrderedDict
import pandas as pd
from ase import Atoms


class DataContainer:
    """
    Data container class, which serves for storing
    all information concerned with Monte Carlo
    simulations performed with mchammer.

    Parameters
    ----------
    atoms : ASE Atoms object
        reference atomic structure associated with the data container

    name_ensemble : str
        name of associated ensemble

    random_seed : int
        random seed used in random number generator

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

    def __init__(self, atoms, ensemble_name: str, random_seed: int):
        """
        Initialize a DataContainer object.
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

        self._metadata['ensemble-name'] = ensemble_name

        self._metadata['date-created'] = datetime.now()
        self._metadata['username'] = getpass.getuser()
        self._metadata['hostname'] = socket.gethostname()

    def add_observable(self, tag: str, obs_type):
        """
        Add an observable to the dict with observables.

        Parameters
        ----------
        tag : str
            name of observable.
        obs_type : type
            type of observable parameter.
        """
        assert isinstance(tag, str), \
            'Observable tag has wrong type (str)'
        assert isinstance(obs_type, (type, list)), \
            'Unknown observable type: {}'.format(obs_type)
        self._observables[tag] = obs_type

    def add_parameter(self, tag: str, value):
        """
        Add parameter of the associated ensemble.

        Parameters
        ----------
        tag : str
            parameter name
        value : int or float or list of int or float
            parameter value
        """
        import copy
        assert isinstance(tag, str), \
            'Parameter tag has the wrong type (str).'
        assert isinstance(value, (int, float, list)), \
            'Unknown parameter type: {}'.format(type(value))
        self._parameters[tag] = copy.deepcopy(value)

    def _update_data(self, data):
        """
        Private method to update data from read function.

        Parameters
        ----------
        data : Pandas DataFrame object
            data of DataContainer object read from file

        Todo
        ----
        Using concatenate here has futher purposes when an implemenation
        of backup will be done.
        """
        self._data = pd.concat([self._data, data])

    def append(self, mctrial: int, record: dict):
        """
        Append data to the data container.

        Parameters
        ----------
        mctrial : int
            current Monte Carlo trial step
        record : dict
            dictionary of tag-value pairs representing observations

        Todo
        ----
        This might be a quite expensive way to add data to the data
        frame. Testing and profiling to be carried out later.
        """
        assert isinstance(mctrial, int), \
            'Monte Carlo trial step has wrong type (int)'
        assert isinstance(record, dict), \
            'Input record has wrong type (dict)'

        row_data = OrderedDict()
        row_data['mctrial'] = mctrial
        row_data.update(record)

        self._data = self._data.append(row_data,
                                       ignore_index=True)

    def get_data(self, tags=None, interval=None, fill_missing=False):
        """
        Returns a list of lists representing the accumulated data for
        the observables specified via tags.

        Parameters
        ----------
        tags : list of str
            tags of the required properties; if None all columns of the data
            frame will be returned

        interval : tuple
            range of trial steps values from which data frame will be returned

        fill_missing : bool
            If True fill missing values forward
        """
        import math

        if tags is None:
            tags = self._data.columns.tolist()
        else:
            for tag in tags:
                assert tag in self._data, \
                    'observable is not part of DataContainer: {}'.format(tag)

        if interval is None:
            data = self._data.loc[:, tags]
        else:
            assert isinstance(interval, tuple), \
                'interval must be a tuple: {}'.format(type(interval))
            assert len(interval) == 2, \
                'interval must contain only a lower and an upper value'
            data = self._data.set_index(self._data.mctrial)
            lower, upper = interval
            data = data.loc[lower:upper, tags]

        if fill_missing:
            data = data.fillna(method='pad')

        data_list = []
        for tag in tags:
            data_column = data.get(tag).tolist()
            data_column = [None if math.isnan(x) else x for x in data_column]
            data_list.append(data_column)

        return data_list

    @property
    def parameters(self):
        """ list : simulation parameters """
        return self._parameters.copy()

    @property
    def observables(self):
        """ dict : observables """
        return self._observables

    @property
    def metadata(self):
        """ dict : metadata associated with data container """
        return self._metadata

    def reset(self):
        """ Reset (clear) data frame of data container """
        self._data = pd.DataFrame()

    def __len__(self):
        """ Number of rows in data frame """
        return len(self._data)

    def get_average(self, start: int, stop: int, tag: str):
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

    @staticmethod
    def restart(filename):
        """
        Restart DataContainer object from file.

        Parameters
        ----------
        filename : str or FileObj
            file from which to read
        """
        pass

    @staticmethod
    def read(filename):
        """
        Read DataContainer object from file.

        Parameters
        ----------
        filename : str or FileObj
            file from which to read
        """
        temp_pkl_file = tempfile.NamedTemporaryFile()
        temp_cvs_file = tempfile.NamedTemporaryFile()

        with tarfile.open(mode='r', name=filename) as tar_file:
            temp_pkl_file.write(tar_file.extractfile('reference-data').read())
            temp_pkl_file.seek(0)
            reference_data = pickle.load(temp_pkl_file)

            dc = DataContainer(reference_data['atoms'],
                               reference_data['metadata']['ensemble-name'],
                               reference_data['parameters']['random-seed'])

            for key in reference_data:
                if key == 'atoms':
                    continue
                if key == 'metadata':
                    for tag, value in reference_data[key].items():
                        dc._metadata[tag] = value
                if key == 'parameters':
                    for tag, value in reference_data[key].items():
                        if tag == 'random-seed':
                            continue
                        dc.add_parameter(tag, value)
                if key == 'observables':
                    for tag, value in reference_data[key].items():
                        dc.add_observable(tag, value)

            temp_cvs_file.write(tar_file.extractfile('runtime-data').read())
            temp_cvs_file.seek(0)
            runtime_data = pd.read_csv(temp_cvs_file)
            dc._update_data(runtime_data)
        return dc

    def write(self, filename):
        """
        Write DataContainer object to file.

        Parameters
        ----------
        filename : str or FileObj
            file to which to write
        """
        self._metadata['date-last-backup'] = datetime.now()
        # Save reference data to a pickle tempfile
        reference_data = {'atoms': self.structure,
                          'observables': self. _observables,
                          'parameters': self._parameters,
                          'metadata': self._metadata}

        temp_pkl_file = tempfile.NamedTemporaryFile()
        with open(temp_pkl_file.name, 'wb') as handle:
            pickle.dump(reference_data, handle,
                        protocol=pickle.HIGHEST_PROTOCOL)

        # Save Pandas data frame as a csv tempfile
        temp_csv_file = tempfile.NamedTemporaryFile()
        self._data.to_csv(temp_csv_file.name)

        with tarfile.open(filename, mode='w') as handle:
            handle.add(temp_csv_file.name, arcname='runtime-data')
            handle.add(temp_pkl_file.name,
                       arcname='reference-data')
