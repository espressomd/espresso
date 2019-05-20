#ifndef LB_WALBERLA_INTERFACE_HPP
#define LB_WALBERLA_INTERFACE_HPP

#include "config.hpp"

#ifdef LB_WALBERLA

#include "boost/optional.hpp"

namespace Walberla {

boost::optional<Utils::Vector3d> get_node_velocity(Utils::Vector3i ind);
boost::optional<double> get_node_density(Utils::Vector3i ind);
boost::optional<double> get_node_is_boundary(Utils::Vector3i ind);
boost::optional<Utils::Vector19d> get_node_pop(Utils::Vector3i ind);

void set_node_velocity(Utils::Vector3i ind, Utils::Vector3d u);
void set_node_density(Utils::Vector3i ind, double density);
void set_node_pop(Utils::Vector3i ind, Utils::Vector19d pop);

} // namespace Walberla
#endif
#endif
