#ifndef LB_WALBERLA_INTERFACE_HPP
#define LB_WALBERLA_INTERFACE_HPP

#include "config.hpp"

#ifdef LB_WALBERLA

#include "boost/optional.hpp"

namespace Walberla {

boost::optional<Utils::Vector3d> get_node_velocity(Utils::Vector3i ind);
boost::optional<double> get_node_density(Utils::Vector3i ind);

void set_node_velocity(Utils::Vector3i ind, Utils::Vector3d u);
void set_node_density(Utils::Vector3i ind, double density);

} // namespace Walberla
#endif
#endif
