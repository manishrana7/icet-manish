#pragma once

#include "LatticeSite.hpp"
#include "Structure.hpp"

namespace icet {

    /// Returns the geometrical radius of a set of lattice sites.
    double getClusterRadius(const std::vector<LatticeSite> &, const Structure &);

}
