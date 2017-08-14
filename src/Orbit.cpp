    #include "Orbit.hpp"
    
    ///Compare operator for automatic sorting in containers
    friend bool operator<(const Orbit &orbit1, const Orbit &orbit2)
    {
        if (orbit1.getRepresentativeCluster() < orbit2.getRepresentativeCluster())
        {
            return true;
        }
        if (orbit1.getRepresentativeCluster() > orbit2.getRepresentativeCluster())
        {
            return false;
        }
        //representative cluster is equal
        //Try comparing length of equivalent sites
        if (orbit1.size() < orbit2.size())
        {
            return true;
        }
        if (orbit1.size() > orbit2.size())
        {
            return false;
        }
        //Both representative cluster and size of equivalent sites are equal.
        //throw error to see if this ever happens

        throw std::runtime_error("Both representative cluster and size of equivalent sites are equal in orbit < comparison");
    }