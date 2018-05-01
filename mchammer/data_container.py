import tarfile
import tempfile
import getpass
import socket
import json
from datetime import datetime
from collections import OrderedDict
import pandas as pd
from ase import Atoms
from ase.io import write as ase_write
from ase.io import read as ase_read
from icet import __version__ as icet_version


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

    Todo
    ----
    atoms is always expected to be given among the arguments so stop
    querying whether atoms is None once it has been fixed in BaseEnsemble
    class.
    """

    def __init__(self, atoms, ensemble_name: str, random_seed: int):
        """
        Initialize a DataContainer object.
        """

        if atoms is not None:
            assert isinstance(atoms, Atoms), \
             'Structure must be provided as ASE Atoms object'
            self.structure = atoms.copy()

        self._observables = []
        self._parameters = OrderedDict()
        self._metadata = OrderedDict()
        self._data = pd.DataFrame(columns=['mctrial'])
        self._data = self._data.astype({'mctrial': int})

        self.add_parameter('seed', random_seed)

        self._metadata['ensemble_name'] = ensemble_name

        self._metadata['date_created'] = \
            datetime.now().strftime('%Y-%m-%dT%H:%M:%S')
        self._metadata['username'] = getpass.getuser()
        self._metadata['hostname'] = socket.gethostname()
        self.metadata['icet_version'] = icet_version

    def add_observable(self, tag: str):
        """
        Add an observable to the dict with observables.

        Parameters
        ----------
        tag : str
            name of observable.
        """
        assert isinstance(tag, str), \
            'Observable tag has wrong type (str)'
        self._observables.append(tag)

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
        # dict to DataFrame to avoid dtype conversion when appending data
        temp_data = pd.DataFrame(row_data, index=[0])
        self._data = self._data.append(temp_data,
                                       ignore_index=True)

    def get_data(self, tags=None, interval=None, fill_missing=False):
        """
        Returns a list of lists representing the accumulated data for
        the observables specified via tags.

        Parameters
        ----------
        tags : list of str
            tags of the required properties; if None all columns of the data
            frame will be returned in lexigraphical order

        interval : tuple
            range of trial steps values from which data frame will be returned

        fill_missing : bool
            If True fill missing values backward
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
            data = data.fillna(method='bfill')

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

    def get_number_of_entries(self, tag=None):
        """
        Return the total number of entries in the column labeled with the
        given observable tag.

        Parameters
        ----------
        tag : str
            name of observable. If None the total number of rows in the Pandas
            data frame will be returned.
        """
        if tag is None:
            return len(self._data)
        else:
            assert tag in self._data, \
                'observable is not part of DataContainer: {}'.format(tag)
            return self._data[tag].count()

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
    def read(infile):
        """
        Read DataContainer object from file.

        Parameters
        ----------
        infile : str or FileObj
            file from which to read
        """
        import os

        if not os.path.isfile(infile):
            raise Exception("File cannot be found")

        temp_atoms_file = tempfile.NamedTemporaryFile()
        temp_json_file = tempfile.NamedTemporaryFile()
        temp_cvs_file = tempfile.NamedTemporaryFile()

        with tarfile.open(mode='r', name=infile) as tar_file:
            temp_atoms_file.write(tar_file.extractfile('atoms').read())
            atoms = ase_read(temp_atoms_file.name, format='xyz')

            temp_json_file.write(tar_file.extractfile('reference_data').read())
            temp_json_file.seek(0)
            reference_data = json.load(temp_json_file)

            dc = DataContainer(atoms,
                               reference_data['metadata']['ensemble_name'],
                               reference_data['parameters']['seed'])

            for key in reference_data:
                if key == 'metadata':
                    for tag, value in reference_data[key].items():
                        if tag == 'ensemble_name':
                            continue
                        dc._metadata[tag] = value
                elif key == 'parameters':
                    for tag, value in reference_data[key].items():
                        if tag == 'seed':
                            continue
                        dc.add_parameter(tag, value)
                elif key == 'observables':
                    for value in reference_data[key]:
                        dc.add_observable(value)

            temp_cvs_file.write(tar_file.extractfile('runtime_data').read())
            temp_cvs_file.seek(0)
            runtime_data = pd.read_csv(temp_cvs_file)
            dc._data = runtime_data

        return dc

    def write(self, outfile):
        """
        Write DataContainer object to file.

        Parameters
        ----------
        outfile : str or FileObj
            file to which to write
        """
        self._metadata['date_last_backup'] = \
            datetime.now().strftime('%Y-%m-%dT%H:%M:%S')

        # Save reference atomic structure
        temp_atoms_file = tempfile.NamedTemporaryFile()
        ase_write(temp_atoms_file.name, self.structure, format='xyz')

        # Save reference data to a json tempfile
        reference_data = {'observables': self._observables,
                          'parameters': self._parameters,
                          'metadata': self._metadata}

        temp_json_file = tempfile.NamedTemporaryFile()
        with open(temp_json_file.name, 'w') as handle:
            json.dump(reference_data, handle)

        # Save Pandas data frame as a csv tempfile
        temp_csv_file = tempfile.NamedTemporaryFile()
        self._data.to_csv(temp_csv_file.name)

        with tarfile.open(outfile, mode='w') as handle:
            handle.add(temp_csv_file.name, arcname='runtime_data')
            handle.add(temp_json_file.name, arcname='reference_data')
            handle.add(temp_atoms_file.name, arcname='atoms')
