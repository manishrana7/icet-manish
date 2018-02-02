from functools import total_ordering

import numpy as np


@total_ordering
class LatticeSite(object):
    """
    Representation of a lattice site.
    This is basically a named tuple with
    index and unitcell_offset as attributes but it
    also defines comparing LatticeSites.
    """

    def __init__(self, index, unitcell_offset):
        """
        LatticeSite initializer.

        parameters
        ----------
        index : int
            Refers to an index in a lattice
        unitcell_offset : ndarray with shape (3,)
            Refers to an unitcell_offset of the unit
            cell vectors.

        """
        self._index = index
        self.unitcell_offset = unitcell_offset

    @property
    def index(self):
        """
        The index refers to an index in a lattice.
        """
        return self._index

    @property
    def unitcell_offset(self):
        """
        Returns the unitcell offset with type ndarray with shape (3,).
        """
        return self._unitcell_offset

    @unitcell_offset.setter
    def unitcell_offset(self, offset):
        if isinstance(offset, (list, tuple)):
            self._unitcell_offset = np.array(offset)
        else:
            self._unitcell_offset = offset
        self.as_list = [self.index, self.unitcell_offset[0], self.unitcell_offset[1], self.unitcell_offset[2]]


    def __eq__(self, other):
        """
        Test equality of this and another LatticeSite

        Will first test index, then unitcell_offset element by element.
        It will return True/False at first difference it spots.
        If no difference is found the two objects are equal.
        """
        return self.as_list == other.as_list


    def __lt__(self, other):
        """
        Test if this Lattice Site is less than the
        other LatticeSite

        Will first test index, then unitcell_offset element by element.
        It will return True/False at first difference it spots.
        If no difference is found the two objects are equal
        and False is returned.
        """
        if self.index != other.index:
            return self.index < other.index
        for offset, offset_other in zip(self.unitcell_offset,
                                        other.unitcell_offset):
            if offset != offset_other:
                return offset < offset_other
        return False

    def __hash__(self):
        """
        Generate a hash based on index and unitcell_offset.
        """
        return hash((self.index, self.unitcell_offset[0],
                     self.unitcell_offset[1], self.unitcell_offset[2]))

    def __str__(self):
        """
        Return  index and unitcell_offset in str format
        """
        return '{} : {}'.format(self.index, self.unitcell_offset)

    def __repr__(self):
        """
        Return  index and unitcell_offset in str format
        """

        return self.__str__()


class LatticeSiteLean(list):

    def __init__(self, index, unitcell_offset):
        lst = [index] + list(unitcell_offset)
        super(LatticeSiteLean, self).__init__(lst)

        @property
        def index(self):
            """
            The index refers to an index in a lattice.
            """
            return self[0]

        @property
        def unitcell_offset(self):
            """
            Returns the unitcell offset with type ndarray with shape (3,).
            """
            return np.array(self[1:], 'i')  # integer array

        @unitcell_offset.setter
        def unitcell_offset(self, offset):
            self.__init__(self.index, offset)


    def __str__(self):
        """
        Return  index and unitcell_offset in str format
        """
        return '{} : {}'.format(self.index, self.unitcell_offset)


    def __repr__(self):
        """
        Return  index and unitcell_offset in str format
        """

        return self.__str__()
