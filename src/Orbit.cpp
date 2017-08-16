#include "Orbit.hpp"

Orbit::Orbit(const Cluster &cluster)
{
    _representativeCluster = cluster;
}

/** 
Returns the number of exactly equal sites in equivalent sites vector

This is used among other things to debug orbits when duplicates is not expected
*/
int Orbit::getNumberOfDuplicates(int verbosity) const
{
    int numberOfEquals = 0;
    for (size_t i = 0; i < _equivalentSites.size(); i++)
    {
        for (size_t j = i + 1; j < _equivalentSites.size(); j++)
        {

            auto i_sites = _equivalentSites[i];
            auto j_sites = _equivalentSites[j];
            //compare the sorted sites
            std::sort(i_sites.begin(), i_sites.end());
            std::sort(j_sites.begin(), j_sites.end());
            if (i_sites == j_sites)
            {
                if (verbosity > 1)
                {
                    std::cout << "Duplicate in orbit: " << i << " " << j << std::endl;
                    if (verbosity > 2)
                    {
                        std::cout << "sites on " << i << std::endl;
                        for (auto i_latnbr : i_sites)
                        {
                            i_latnbr.print();
                        }
                        std::cout << "sites on " << j << std::endl;
                        for (auto j_latnbr : j_sites)
                        {
                            j_latnbr.print();
                        }
                    }
                }
                numberOfEquals++;
            }
        }
    }
    return numberOfEquals;
}