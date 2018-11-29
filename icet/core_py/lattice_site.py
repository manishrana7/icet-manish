from functools import total_ordering

import numpy as np
from icet.core.lattice_site import LatticeSite as LatticeSite_cpp

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
        self.as_list = [self.index, self.unitcell_offset[0], self.unitcell_offset[1], self.unitcell_offset[2]]

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
        if isinstance(other, LatticeSite_cpp):
            if self.index != other.index:
                return False
            for offset, offset_other in zip(self.unitcell_offset,
                                            other.unitcell_offset):
                if offset != offset_other:
                    return False
            return True
        elif isinstance(other, LatticeSiteLean):
            return self.as_list == other.lst
        else:
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

@total_ordering
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
        return np.array(self[1:])  # integer array

    @unitcell_offset.setter
    def unitcell_offset(self, offset):
        self.__init__(self.index, offset)


    def __str__(self):
        """
        Return  index and unitcell_offset in str format
        """
        return '{} : {}'.format(self.index, self.unitcell_offset)

    def __eq__(self, other):
        """
        Test equality of this and another LatticeSite

        Will first test index, then unitcell_offset element by element.
        It will return True/False at first difference it spots.
        If no difference is found the two objects are equal.
        """
        if isinstance(other, LatticeSite_cpp) or isinstance(other, LatticeSite_old):
            if self.index != other.index:
                return False
            for offset, offset_other in zip(self.unitcell_offset,
                                            other.unitcell_offset):
                if offset != offset_other:
                    return False
            return True        
        else:
            return self[:] == other[:]


    def __lt__(self, other):
        """
        Test equality of this and another LatticeSite

        Will first test index, then unitcell_offset element by element.
        It will return True/False at first difference it spots.
        If no difference is found the two objects are equal.
        """
        if isinstance(other, LatticeSite_cpp) or isinstance(other, LatticeSite_old):
            if self.index != other.index:
                self.index < other.index
            for offset, offset_other in zip(self.unitcell_offset,
                                            other.unitcell_offset):
                if offset != offset_other:
                    return offset < offset_other
            return False        
        else:
            return self[:] < other[:]



    def __repr__(self):
        """
        Return  index and unitcell_offset in str format
        """

        return self.__str__()

<<<<<<< HEAD
    def __hash__(self):
        """
        Generate a hash based on index and unitcell_offset.
        """
        return hash((self.index, self.unitcell_offset[0],
                     self.unitcell_offset[1], self.unitcell_offset[2]))
=======
    def __add__(self, other):
        """
        Define the add operator.
        Allowed values:
        * type ndarray with shape(3,)
        """
        if not isinstance(other, type(np.array([1, 1, 1]))) or len(other) != 3:
            raise TypeError("Adding orbit with {}".format(other))
        site = LatticeSite(self.index, self.unitcell_offset + other)
        return site

    def __sub__(self, other):
        """
        Define the subtract operator.
        Allowed values:
        * type ndarray with shape(3,)
        """
        if not isinstance(other, type(np.array([1, 1, 1]))) or len(other) != 3:
            raise TypeError("Adding orbit with {}".format(other))
        site = LatticeSite(self.index, self.unitcell_offset - other)
        return site


def cmp_mbnl_lattice_site_list(first, second):
    """
    Comparer for list of lattice sites.
    First compare len of lists then do the normal,
    lexicographical comparing.
    """
    if len(first[0]) != len(second[0]):
        return len(first[0]) < len(second[0])
    else:
        return first < second


def cmp_list_of_lattice_sites(first, second):
    """
    Comparer for list of lattice sites.
    First compare len of lists then do the normal,
    lexicographical comparing.
    """
    if len(first) != len(second):
        return len(first) < len(second)
    else:
        for site1, site2 in zip(first, second):
            if site1 != site2:
                return site1 < site2
    return False


def cmp_to_key(mycmp):
    'Convert a cmp= function into a key= function'
    class K:
        def __init__(self, obj, *args):
            self.obj = obj

        def __lt__(self, other):
            return mycmp(self.obj, other.obj)

        def __gt__(self, other):
            return mycmp(self.obj, other.obj)

        def __eq__(self, other):
            return mycmp(self.obj, other.obj)

        def __le__(self, other):
            return mycmp(self.obj, other.obj)

        def __ge__(self, other):
            return mycmp(self.obj, other.obj)

        def __ne__(self, other):
            return mycmp(self.obj, other.obj)
    return K
>>>>>>> master
