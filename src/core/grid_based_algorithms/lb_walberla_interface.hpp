#ifndef LB_WALBERLA_INTERFACE_HPP
#define LB_WALBERLA_INTERFACE_HPP

#include "config.hpp"

#ifdef LB_WALBERLA

#include "boost/optional.hpp"

namespace Walberla {

boost::optional<Utils::Vector3d> get_node_velocity(Utils::Vector3i ind);

void set_node_velocity(Utils::Vector3i ind, Utils::Vector3d u);

} // namespace Walberla
#endif
#endif
