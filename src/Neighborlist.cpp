#include "Neighborlist.hpp"
#include "Structure.hpp"
#include "Vector3dCompare.hpp"

Neighborlist::Neighborlist(const double cutoff)
{
    _cutoff = cutoff;
}

void Neighborlist::build(const Structure &conf)
{
    //resize the offsets and indices to number of atoms
    int nbrOfSites = conf.size();
    latticeIndices.resize(nbrOfSites);
    offsets.resize(nbrOfSites);
    _neighbors.resize(nbrOfSites);
    Matrix3d cellInverse = conf.getCell().inverse();
    std::vector<int> unitCellExpanse(3);
    for (int i = 0; i < 3; i++)
    {
        if (conf.hasPBC(i))
        {
            auto v = cellInverse.col(i);
            double dotProduct = v.dot(v);
            double h = 1.0 / sqrt(dotProduct);
            int n = (int)(1.0 * _cutoff / h) + 1;
            unitCellExpanse[i] = n;
        }
        else
        {
            unitCellExpanse[i] = 0;
        }
    }

    for (int n1 = 0; n1 < unitCellExpanse[0] + 1; n1++)
    {

        for (int n2 = -unitCellExpanse[1]; n2 < unitCellExpanse[1] + 1; n2++)
        {
            for (int n3 = -unitCellExpanse[2]; n3 < unitCellExpanse[2] + 1;
                 n3++)
            {

                if (n1 == 0 and (n2 < 0 or (n2 == 0 and n3 < 0)))
                {
                    //continue;
                }

                for (int m = -1; m < 2; m += 2)
                {
                    Vector3d extVector(n1 * m, n2 * m, n3 * m);

                    for (int i = 0; i < nbrOfSites; i++)
                    {

                        for (int j = 0; j < nbrOfSites; j++)
                        {
                            Vector3d noOffset(0, 0, 0);

                            double distance_ij = conf.getDistance2(i, noOffset, j, extVector);

                            if (distance_ij <= _cutoff && distance_ij > 2 * DISTTOL)
                            {
                                LatticeSite neighbor = LatticeSite(j, extVector);
                                auto find_neighbor = std::find(_neighbors[i].begin(),_neighbors[i].end(), neighbor);
                                if(find_neighbor == _neighbors[i].end())
                                {
                                    _neighbors[i].push_back(neighbor);
                                }
                            }
                        }
                    }
                } // end m loop
            }
        }
    } // end n loop
    for(auto &nbr : _neighbors )
    {
        std::sort(nbr.begin(), nbr.end());
    }

}
