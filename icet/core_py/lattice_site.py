from functools import total_ordering


@total_ordering
class LatticeSite(object):
    """
    Representation of a lattice site.
    This is basically a named tuple with
    index and offset as attributes but it
    also defines comparing LatticeSites.
    """

    def __init__(self, index, offset):
        """
        LatticeSite initializer.

        parameters
        ----------
        index : int
            Refers to an index in a lattice
        offset : list of ints with length 3
            Refers to an offset of the unit
            cell vectors.

        """
        self._index = index
        self._offset = offset

    @property
    def index(self):
        """
        The index refers to an index in a lattice.
        """
        return self._index

    @property
    def offset(self):
        """
        The index refers to an index in a lattice.
        """
        return self._offset

    def __eq__(self, other):
        """
        Test equality of this and another LatticeSite

        Will first test index, then offset element by element.
        It will return True/False at first difference it spots.
        If no difference is found the two objects are equal.
        """
        if self.index != other.index:
            return False
        for i in range(3):
            if self.offset[i] != other.offset[i]:
                return False
        return True

    def __lt__(self, other):
        """
        Test if this Lattice Site is less than the
        other LatticeSite

        Will first test index, then offset element by element.
        It will return True/False at first difference it spots.
        If no difference is found the two objects are equal
        and False is returned.
        """
        if self.index != other.index:
            return self.index < other.index
        for offset, offset_other in zip(self.offset, other.offset):
            if offset != offset_other:
                return offset < offset_other
        return False

    def __hash__(self):
        """
        Generate a hash based on index and offset.
        """
        return hash((self.index, self.offset[0],
                     self.offset[1], self.offset[2]))

    def __str__(self):
        """
        Return  index and offset in str format
        """
        nice_str = str(self.index) + " "
        nice_str += str(self.offset)
        return nice_str
