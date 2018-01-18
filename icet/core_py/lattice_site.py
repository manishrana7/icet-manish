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
        unitcell_offset : list of ints with length 3
            Refers to an unitcell_offset of the unit
            cell vectors.

        """
        self._index = index
        if isinstance(unitcell_offset, list):
            self._unitcell_offset = np.array(unitcell_offset)
        else:
            self._unitcell_offset = unitcell_offset

    @property
    def index(self):
        """
        The index refers to an index in a lattice.
        """
        return self._index

    @property
    def unitcell_offset(self):
        """
        The index refers to an index in a lattice.
        """
        return self._unitcell_offset

    def __eq__(self, other):
        """
        Test equality of this and another LatticeSite

        Will first test index, then unitcell_offset element by element.
        It will return True/False at first difference it spots.
        If no difference is found the two objects are equal.
        """
        if self.index != other.index:
            return False
        for i in range(3):
            if self.unitcell_offset[i] != other.unitcell_offset[i]:
                return False
        return True

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
