#pragma once

#include <vector>
#include "LatticeSite.hpp"
#include "Structure.hpp"

namespace icet {

    /// Returns the geometrical radius of a cluster.
    double getGeometricalRadius(const std::vector<LatticeSite> &, const Structure &);

}
