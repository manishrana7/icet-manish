"""
This module provides a Python interface to the ManyBodyNeighborList class.
"""
import numpy as np

from lattice_site import LatticeSite


class ManyBodyNeighborList():
    """
    Design approach:
    input pair neighbors and calculate higher order neighbors
    using set intersection.
    """
    def __init__(self, parent=None):
        self._latticeNeighbors = None

    def build(self, neighborLists, index, saveBothWays):
        """
        @details
            This function uses @a neighbor_lists to construct all possible
            neighbors up to the given order. The output will be:
            @code{.cpp}
            std::vector<std::pair<originalNeighbors, manyNeigbhors>>
            @endcode

            The many body neigbhors can be retrieved by doing:
            @code{.cpp}
            for (const auto nbr : manyBodyNeighborIndices)
            {
                std::vector<std::pair<int,Vector3d>> neighbors = nbr.first;
                // this are the first orignal neighbors
                for(const auto manynbr : nbr.second)
                {
                    many_body_neigbhor = neighbors;
                    many_body_neigbhor.append(manynbr);
                }
            }
            @endcode

            This means that if @a originalNeigbhors.size()==2 then for each lattice site in
                @a manyNeigbhors
            you can combine it with @a originalNeigbhors to get all triplets that have these
                first two original neighbors (lattice indices).

            @param neighborLists list of neighbor lists
            @param index
            @param saveBothWays if true then both @a i,j,k and @a j,i,k etc.. will be saved;
                otherwise only @a i,j,k will be saved if @a i<j<k.
        """
        if not neighborLists:
            raise RuntimeError("Error: neigbhorlist vector is empty in ManyBodyNeighborList::build")
        # if (neighborLists.empty())
        # {
        #     throw std::runtime_error("Error: neigbhorlist vector is empty in ManyBodyNeighborList::build");
        # }
        ### ->manyBodyNeighborIndices = [([], [])]
        # std::vector<std::pair<std::vector<LatticeSite>, std::vector<LatticeSite>>> manyBodyNeighborIndices;
        manyBodyNeighborIndices = []
        self.addSinglet(index, manyBodyNeighborIndices)
        self.addPairs(index, neighborLists[0], manyBodyNeighborIndices, saveBothWays)
        # addSinglet(index, manyBodyNeighborIndices);
        # addPairs(index, neighborLists[0], manyBodyNeighborIndices, saveBothWays);
        currentOriginalNeighbors = []
        for c in range(2, len(neighborLists)+2):
            Ni = neighborLists[c-2][index]
            zeroVector = np.array([0.]*3)
            currentOriginalNeighbors.append(LatticeSite(index, zeroVector))
            self.combineToHigherOrder(neighborLists[c-2], manyBodyNeighborIndices, Ni, currentOriginalNeighbors, saveBothWays, c)

        _latticeNeighbors = manyBodyNeighborIndices
        return manyBodyNeighborIndices

        # for (size_t c = 2; c < neighborLists.size() + 2; c++)
        # {
        #     //auto Ni = neighborLists[c - 2].getNeighbors(index);
        #     auto Ni = neighborLists[c - 2][index];
        #     Vector3d zeroVector = {0.0, 0.0, 0.0};
        #     std::vector<LatticeSite> currentOriginalNeighbors;
        #     currentOriginalNeighbors.push_back(LatticeSite(index, zeroVector));  // index is always first index

        #     combineToHigherOrder(neighborLists[c - 2], manyBodyNeighborIndices, Ni, currentOriginalNeighbors, saveBothWays, c);
        # }
        # _latticeNeighbors = manyBodyNeighborIndices;
        # return manyBodyNeighborIndices;

    def combineToHigherOrder(self, nl, manyBodyNeighborIndices, Ni, currentOriginalNeighbors,
                             saveBothWays, maxOrder):
        """
        This will use the matrix of equivalent sites @a manyBodyNeighborIndices together with
            neighbor list @a nl to construct the distinct and indistinct sets of points.

        The output will be std::vector<std::vector<std::vector<LatticeSite>>>.
        The next outer vector contains the set of indistinct set of lattice neighbors.


        Algorithm
        =========

        The algorithm works by taking the first column of the matrix of equivalent sites
        and it will take Ni (lattice neighbors of i) and find the intersection of Ni and col1,
        intersection(Ni, col1) = Ni_pm
        all j in Ni_c1 are then within the cutoff of site i, then depending on the order all
        the pairs/triplets will be constructed from the lattice neighbors in Ni_pm.

        This will be repeated for all the sites in the neighbor_list. In the end all
        the pair/triplet terms will have been generated from col1.

        Then you will take the first vector<LatticeSites> and find the rows of these LatNbrs
        in col1 of the permutation matrix.
        then you traverse through all the columns in the permutation matrix.
        These vector<latticeNeighbors> are then indistinct from the original.
        Note that here a validity check is made to ensure that atleast one LatticeNeigbhor
        originate in the original lattice (unitcell offset = [0,0,0]) otherwise we are including
        a  "ghost cluster".

        The new vector<latticeNeighbors> found when traversing the columns are likely to have
        been found from the combinations in Ni_pm and these must then
        be removed/overlooked when moving to the next vector<LatticeSite>.

        for each j in Ni construct the intersect of N_j and N_i = N_ij.
        all k in N_ij are then neighbors with i,j
        what is saved is then i,j and N_ij up to the desired order "maxorder"
        """
        # void ManyBodyNeighborList::combineToHigherOrder(const std::vector<std::vector<LatticeSite>> &nl,
        #                                                 std::vector<std::pair<std::vector<LatticeSite>,
        #                                                 std::vector<LatticeSite>>> &manyBodyNeighborIndices,
        #                                                 const std::vector<LatticeSite> &Ni,
        #                                                 std::vector<LatticeSite> &currentOriginalNeighbors,
        #                                                 bool saveBothWays,
        #                                                 const size_t maxOrder)
        # {

        for j in Ni:
            cont = False
            if saveBothWays:
                if len(currentOriginalNeighbors) > 1:
                    if j < currentOriginalNeighbors[-1]:
                        cont = True
            else:
                if j < currentOriginalNeighbors[-1]:
                    cont = True
            if cont:
                continue
            originalNeighborCopy = currentOriginalNeighbors
            originalNeighborCopy.append(j)
            Nj = nl[j.index]
            #translate the neighbors
            self.translateAllNi(Nj, j.unitcell_offset)

            if not saveBothWays:
                Nj = self._getFilteredNj(Nj, j)
            intersection_N_ij = self.getIntersection(Ni, Nj)
            if len(originalNeighborCopy) +1 < maxOrder:
                self.combineToHigherOrder(nl, manyBodyNeighborIndices, intersection_N_ij, originalNeighborCopy, saveBothWays, maxOrder)
            if intersection_N_ij and len(originalNeighborCopy) == (maxOrder-1):
                manyBodyNeighborIndices.append((originalNeighborCopy, intersection_N_ij))


        #     for (const auto &j : Ni)
        #     {
        #         //if j is smaller than last added site then continue
        #         // if bothways = True then don't compare to first
        #         bool cont = false;

        #         if (saveBothWays)
        #         {
        #             if (currentOriginalNeighbors.size() > 1)
        #             {
        #                 if (j < currentOriginalNeighbors.back())
        #                 {
        #                     cont = true;
        #                 }
        #             }
        #         }
        #         else
        #         {
        #             if (j < currentOriginalNeighbors.back())
        #             {
        #                 cont = true;
        #             }
        #         }
        #         if (cont)
        #         {
        #             continue;
        #         }
        #         // if ((!saveBothWays && currentOriginalNeighbors.size() == 1) && j < currentOriginalNeighbors.back())
        #         // {
        #         //     continue;
        #         // }


            #         auto originalNeighborCopy = currentOriginalNeighbors;
            #         originalNeighborCopy.push_back(j); // put j in originalNeigbhors

            #         auto Nj = nl[j.index()];

            #         //translate the neighbors
            #         translateAllNi(Nj, j.unitcellOffset());

            #exclude smaller neighbors

        #         if (!saveBothWays)
        #         {
        #             Nj = getFilteredNj(Nj, j);
        #         }

            # construct the intersection
        #         const auto intersection_N_ij = getIntersection(Ni, Nj);

        #         if (originalNeighborCopy.size() + 1 < maxOrder)
        #         {
        #             combineToHigherOrder(nl, manyBodyNeighborIndices, intersection_N_ij, originalNeighborCopy, saveBothWays, maxOrder);
        #         }

        #         if (intersection_N_ij.size() > 0 && originalNeighborCopy.size() == (maxOrder - 1))
        #         {
        #             manyBodyNeighborIndices.push_back(std::make_pair(originalNeighborCopy, intersection_N_ij));
        #         }
        #     }
        # }

    def getIntersection(self, Ni, Nj):
        """
        @details Return the lattice sites that appear in two list of lattice sites.
        @param Ni list of lattice sites
        @param Nj another list of lattice sites        
        """
        intersection = set(Ni).intersection(set(Nj))

        return list(intersection)

        #    std::vector<LatticeSite> getIntersection(const std::vector<LatticeSite> &Ni, const std::vector<LatticeSite> &Nj)
        #     {
        #         std::vector<LatticeSite> N_intersection;
        #         N_intersection.reserve(Ni.size());
        #         std::set_intersection(Ni.begin(), Ni.end(),
        #                               Nj.begin(), Nj.end(),
        #                               std::back_inserter(N_intersection));
        #         return N_intersection;
        #     }

    def translateAllNi(self, Ni=None, offset=None):
        """
        Offsets all indice, offsets pairs in Ni with the input offset, e.g:
        For all j in Ni:
        offset j.offset with "unitCellOffset"
        """
        # void ManyBodyNeighborList::translateAllNi(std::vector<LatticeSite> &Ni, const Vector3d &offset) const
        # {
        #     for (auto &latticeSite : Ni)
        #     {
        #         latticeSite.addUnitcellOffset(offset);
        #     }
        # }
        for latticeSite in Ni:
            latticeSite.add_unitcell_offset(offset)

    def getNumberOfSites(self, index):
        """
        Returns number of manybodies one can make from _latticeNeighbors[index]
        """
        # size_t ManyBodyNeighborList::getNumberOfSites(const unsigned int index) const
        # {
        #     //std::vector<std::pair<std::vector<LatticeSite>, std::vector<LatticeSite>>> _latticeNeighbors;
        #     return _latticeNeighbors[index].second.size();
        # }        
        return len(self._latticeNeighbors()[index][1])

    def getSites(self, firstIndex=None, secondIndex=None):
        """
        Return the many_body neighbor at 'firstIndex'  and 'secondIndex'
        in _latticeNeighbors and _latticeNeighbors[firstIndex] respectively
        """

            # std::vector<LatticeSite> ManyBodyNeighborList::getSites(const unsigned int &firstIndex,
            #                                                             const unsigned int &secondIndex) const
            # {
            #     std::vector<LatticeSite> sites = _latticeNeighbors[firstIndex].first;
            #     //If zero then this is a singlet
            #     if (getNumberOfSites(firstIndex) > 0)
            #     {
            #         sites.push_back(_latticeNeighbors[firstIndex].second[secondIndex]);
            #     }
            #     return sites;
            # }
        sites = self._latticeNeighbors[firstIndex][0]
        # if zero - this is a single
        if self.getNumberOfSites(firstIndex) > 0:
            sites.append(self._latticeNeighbors[firstIndex][1][secondIndex])

        return sites

    def addSinglet(self, index, manyBodyNeighborIndices):
        """
        Adds singlet from the index to manyBodyNeighborIndices
        """
        # Vector3d zeroVector = {0.0, 0.0, 0.0};
        # LatticeSite latticeNeighborSinglet = LatticeSite(index, zeroVector);
        # std::vector<LatticeSite> singletLatticeSites;
        # singletLatticeSites.push_back(latticeNeighborSinglet);

        # std::vector<LatticeSite> latticeSitesEmpty;
        # manyBodyNeighborIndices.push_back(std::make_pair(singletLatticeSites, latticeSitesEmpty));
        zeroVector = np.array([0.]*3)
        latticeNeighborSinglet = LatticeSite(index, zeroVector)
        singleLatticeSites = [latticeNeighborSinglet]
        manyBodyNeighborIndices.append((singleLatticeSites, []))


    def addPairs(self, lattice_index, neighborList, manyBodyNeighborIndices, saveBothWays):
        """
        Add all pairs originating from index using neighbor_list
        """
        zeroVector = np.array([0.0]*3)
        latticeNeighborIndex = LatticeSite(lattice_index, zeroVector)
        firstSite = [latticeNeighborIndex]
        Ni = neighborList[lattice_index]
        # exclude smaller neighbors
        if not saveBothWays:
            Ni = self._getFilteredNj(Ni, latticeNeighborIndex)
        if not Ni:
            return
        manyBodyNeighborIndices.append((firstSite, Ni))

        # {
        #     Vector3d zeroVector = {0.0, 0.0, 0.0};
        #     LatticeSite latticeNeighborIndex = LatticeSite(index, zeroVector);

        #     std::vector<LatticeSite> firstSite = {latticeNeighborIndex};
        #     std::vector<LatticeSite> Ni = neighborList[index];
        #     //exclude smaller neighbors
        #     if (!saveBothWays)
        #     {
        #         Ni = getFilteredNj(Ni, latticeNeighborIndex);
        #     }

        #     if (Ni.size() == 0)
        #     {
        #         return;
        #     }
        #     manyBodyNeighborIndices.push_back(std::make_pair(firstSite, Ni));
        # } 

    def _getFilteredNj(self, N_j, j):
        """
        Since N_j is always sorted then simply search for first k in N_j that have k>= j
        and then filtered are from indexof(k) to end()
        """
        # std::vector<LatticeSite> ManyBodyNeighborList::getFilteredNj(const std::vector<LatticeSite> &N_j, const LatticeSite &j) const
        # {
        #     auto first = std::upper_bound(N_j.begin(), N_j.end(), j);

        #     std::vector<LatticeSite> ret(first, N_j.end());
        #     return ret;
        # }
        import bisect
        first = bisect.bisect_right(N_j, j)
        ret = N_j[first:]
        return ret


from icet.core.structure import Structure
from icet.core.lattice_site import LatticeSite
from icet.core.neighbor_list import get_neighbor_lists
from ase.build import bulk

structure = bulk('Ni', 'hcp', a=3.0).repeat([2, 2, 1])
cutoffs = [5.0, 5.0]
position_tolerance = 1e-5
mbnl = ManyBodyNeighborList()
structure = Structure.from_atoms(structure)
neighbor_lists = get_neighbor_lists(structure, cutoffs, position_tolerance)
mbnl_size = len(mbnl.build(neighbor_lists, 0, False))

for index in range(len(structure)):
    mbnl.build(neighbor_lists, index, True)
for index in range(len(structure)):
    target = tuple(([LatticeSite(index, [0., 0., 0.])], []))
    singlet = mbnl.build(neighbor_lists, index, False)[0]
    assert (singlet == target)

for index in range(1, len(structure)):
    nl = mbnl.build(neighbor_lists, index, False)
    #assert mbnl_size == len(nl)
    #assert mbnl_size != len(mbnl.build(neighbor_lists, index, False))

for index in range(1, len(structure)):
    nl = mbnl.build(neighbor_lists, index, True)
    assert mbnl_size == len(nl)
    #assert mbnl_size == len(mbnl.build(neighbor_lists, index, True))

index = 0
nl_neighbors = neighbor_lists[0][0]
target = tuple(([LatticeSite(index, [0., 0., 0.])], nl_neighbors))
pairs = mbnl.build(neighbor_lists, index, True)[1]
assert pairs == target

index = 0
high_order_neighbors = mbnl.build(neighbor_lists, index, False)[2]

target = ([LatticeSite(0, [0, 0, 0]), LatticeSite(0, [0, 0, 1])],
            [LatticeSite(1, [0, 0, 0]),
            LatticeSite(3, [0, -1, 0]),
            LatticeSite(5, [-1, -1, 0]),
            LatticeSite(5, [-1, 0, 0]),
            LatticeSite(5, [0, 0, 0]),
            LatticeSite(7, [-1, -1, 0])])

assert target == high_order_neighbors

lattice_sites = []
lattice_sites.append(LatticeSite(0, [0, 0, 0]))
lattice_sites.append(LatticeSite(0, [1, 0, 0]))
lattice_sites.append(LatticeSite(1, [0, 0, 0]))
lattice_sites.append(LatticeSite(3, [0, 0, 0]))

lattice_sites2 = []
lattice_sites2.append(LatticeSite(0, [0, 0, 0]))
lattice_sites2.append(LatticeSite(0, [1, 0, 0]))

intersection = mbnl.calculate_intersection(
    lattice_sites, lattice_sites2)

assert sorted(intersection) == [LatticeSite(0, [0, 0, 0]), LatticeSite(0, [1, 0, 0])]

stru = structure.copy()
stru.set_pbc([False])
neighbor_lists = get_neighbor_lists(
    Structure.from_atoms(stru), cutoffs, position_tolerance)

mbnl = ManyBodyNeighborList()

target = [([LatticeSite(0, [0, 0, 0])], []),
            ([LatticeSite(0, [0, 0, 0])],
            [LatticeSite(1, [0, 0, 0]),
            LatticeSite(2, [0, 0, 0]),
            LatticeSite(4, [0, 0, 0]),
            LatticeSite(5, [0, 0, 0]),
            LatticeSite(6, [0, 0, 0])]),
            ([LatticeSite(0, [0, 0, 0]), LatticeSite(1, [0, 0, 0])],
            [LatticeSite(2, [0, 0, 0]), LatticeSite(4, [0, 0, 0]),
            LatticeSite(5, [0, 0, 0]), LatticeSite(6, [0, 0, 0])]),
            ([LatticeSite(0, [0, 0, 0]), LatticeSite(2, [0, 0, 0])],
            [LatticeSite(6, [0, 0, 0])]),
            ([LatticeSite(0, [0, 0, 0]), LatticeSite(4, [0, 0, 0])],
            [LatticeSite(5, [0, 0, 0]), LatticeSite(6, [0, 0, 0])]),
            ([LatticeSite(0, [0, 0, 0]), LatticeSite(5, [0, 0, 0])],
            [LatticeSite(6, [0, 0, 0])])]

neighbors_non_pbc = mbnl.build(neighbor_lists, 0, False)

for k, latt_neighbors in enumerate(neighbors_non_pbc):
    assert target[k] == latt_neighbors


stru = bulk('Al', 'sc', a=4.0).repeat(4)
stru.set_pbc(False)

neighbor_lists = get_neighbor_lists(
    Structure.from_atoms(stru), cutoffs, position_tolerance)

mbnl = ManyBodyNeighborList()
# atomic indices located at the corner of structure
corner_sites = [0, 3, 12, 15, 48, 51, 60, 63]
for index in corner_sites:
    lattice_neighbor = mbnl.build(neighbor_lists,
                                    index, True)
    # check pairs
    assert len(lattice_neighbor[1][1]) == 3
    # not neighbors besides above pairs
    #with self.assertRaises(IndexError):
    #    lattice_neighbor[2]