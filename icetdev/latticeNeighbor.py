from _icetdev import LatticeSite



def __latNbr_repr(self):
    """
    represent the Lattice Neighbor with index and unitcell offset
    """
    co = self.unitcellOffset
    rep = str(self.index) +" : " + str(co)
    return rep

LatticeSite.__repr__ = __latNbr_repr

