import pandas as pd
from datetime import date as dt

from ase import Atoms
from collections import OrderedDict


class DataContainer(object):
    """
    Data Container class which serves to store
    all information concerned to Monte Carlo
    calculation performed with mchammer.

    Todo
    ----
    * Improvements required once integrated with
      Base Ensemble class.

    """

    def __init__(self, base_ensemble):
        """
        Initialize a DataContainer object.

        Attributes:
        ----------
        base_ensemble : BaseEnsemble object

        """
        self._base_ensemble = base_ensemble

        self._observables = OrderedDict()
        self._parameters = OrderedDict()
        self._metadata = OrderedDict()

        self._row_data = OrderedDict()
        self._data = pd.DataFrame()

    def add_structure(self, atoms):
        """
        Add a reference atomic structure to
        data container.

        Parameters :
        ----------
        atoms : ASE Atoms object

        """
        atoms_copy = atoms.copy()

        msg = 'Structure has not ASE atoms type'
        assert isinstance(atoms_copy, Atoms), msg

        self.structure = atoms_copy

    def add_observable(self, tag, property_type):
        """
        Add observable to list of observables
        which serves as headers of the data frame
        containing the runtime data.

        Parameters :
        ----------
        tag : str
            tag used to name observable

        property_type : type
            type of observable parameter

        """
        msg = 'Observable tag has wrong type (str)'
        assert isinstance(tag, str), msg

        msg = 'Unknow observable type'
        assert isinstance(property_type, (type, list)), msg

        self._observables[tag] = property_type

    def add_parameter(self, tag, property_value):
        """
        Add parameter required for initial setting
        of the associated ensemble.

        Parameters :
        ----------
        tag : str
            name of the parameter

        property_value : int, float
            value of the parameter

        """
        msg = 'Parameter tag has wrong type (str)'
        assert isinstance(tag, str), msg

        msg = 'Unknow observable type'
        assert isinstance(property_value, (int, float, list)), msg

        self._parameters[tag] = property_value

    def add_metadata(self, tag, input_value):
        """
        Add descriptive metadata (author, hostname, time).

        Parameters :
        ----------
        tag : str
            name of the metadata

        input_value : str, date
            value of the input metadata

        Todo
        ----
        * Include name of ensemble, observer, keywords

        """
        msg = 'Metadata tag has wrong type (str)'
        assert isinstance(tag, str), msg

        msg = 'Unknow metadata type'
        assert isinstance(input_value, (str, dt)), msg

        self._metadata[tag] = input_value

    def append(self, mcstep, tag, input_value):
        """
        Append runtime data to data frame.

        Parameters :
        ----------

        mcstep : int
            current Monte Carlo step

        tag : str
            name of parameter or observable

        input_value: int, float, list
            value of parameter or observable to be appended
            to data container

        """
        msg = 'Monte Carlo step has wrong type (int)'
        assert isinstance(mcstep, int), msg

        msg = 'Input tag has wrong type (str)'
        assert isinstance(tag, str), msg

        msg = 'Input value has unknown type'
        assert isinstance(input_value, (int, float, list)), msg

        if not len(self._row_data) == 0:
            if not mcstep == self._row_data['mcstep']:
                # append row
                self._data = self._data.append(self._row_data,
                                               ignore_index=True)
                # restart row
                self._row_data.clear()

        self._row_data['mcstep'] = mcstep
        self._row_data[tag] = input_value

    def get_data(self, key_tags):
        """
        Returns a list of lists with the current
        data stored in the Data Frame.

        Parameters :
        ----------

        keys : list of str
            list with tags of the required properties.

        """
        msg = 'observable/parameter is not listed in DataContainer'
        data_list = []
        for row in range(len(self._data)):
            data_row = []
            for key in key_tags:
                assert key in self._data, msg
                data_row.append(self._data.loc[row, key])
            data_list.append(data_row)

        return data_list

    @property
    def parameters(self):
        """
        Return added parameters.
        """
        return self._parameters

    @property
    def observables(self):
        """
        Return added observables.
        """
        return self._observables

    @property
    def metadata(self):
        """
        Return added metadata.
        """
        return self._metadata

    def reset(self):
        """
        Restart data frame in data container.
        """
        self._data = pd.DataFrame()

    def __len__(self):
        """
        Return number of existing rows in data frame.
        """
        return len(self._data)

    def get_average(self):
        """
        Return average over MonteCarlo steps of parameter
        and/or observable.
        """
        pass

    def load(self):
        """
        Load data container from file.

        Todo
        ----
        * To be implemented in a different issue.

        """
        pass

    def write(self):
        """
        Write data container to file.

        Todo
        ----
        * To be implemented in a different issue.

        """
        pass
