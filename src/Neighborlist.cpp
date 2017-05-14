#include "Neighborlist.hpp"
#include "Structure.hpp"

Neighborlist::Neighborlist(const double cutoff)
{
    _cutoff = cutoff;
}

void Neighborlist::build(const Structure &conf)
{

    //resize the offsets and indices to number of atoms
    int nbrOfSites = conf.getElements().size();
    latticeIndices.resize(nbrOfSites);
    offsets.resize(nbrOfSites);

    Matrix3d cellInverse = conf.get_cell().inverse();
    std::vector<int> unitCellExpanse(3);
    for (int i = 0; i < 3; i++)
    {
        if (conf.has_pbc(i))
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
                    continue;
                }

                for (int m = -1; m < 2; m += 2)
                {
                    Vector3d extVector(n1 * m, n2 * m, n3 * m);
                    Vector3d displacement = extVector.transpose() * conf.get_cell();

                    for (int i = 0; i < nbrOfSites; i++)
                    {

                        for (int j = i; j < nbrOfSites; j++)
                        {
                            Vector3d noOffset(0, 0, 0);

                            double distance_ij = conf.getDistance2(i, noOffset, j,extVector);

                            if (distance_ij <= _cutoff && distance_ij > 2 * DISTTOL)
                            {

                                //i have a neighbour j, where atom j sits in cell (n1,n2,3)
                                //add j and offset to i'th neighborlist
                                auto find_indice = std::find(latticeIndices[i].begin(),latticeIndices[i].end(),j);
                                auto find_offsets = std::find(offsets[i].begin(), offsets[i].end(), extVector);
                                if(find_indice == latticeIndices[i].end() || find_offsets == offsets[i].end())
                                {
                                    latticeIndices[i].push_back(j);
                                    offsets[i].push_back(extVector);
                                }
                            }
                        }
                    }
                } // end m loop
            }
        }
    } // end n loop
}
