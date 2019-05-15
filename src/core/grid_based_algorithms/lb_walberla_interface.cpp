#include "config.hpp"

#ifdef LB_WALBERLA

#include "MpiCallbacks.hpp"
#include "lb_walberla_instance.hpp"
#include "utils/Vector.hpp"

#include "boost/optional.hpp"

namespace Walberla {

boost::optional<Utils::Vector3d> get_node_velocity(Utils::Vector3i ind) {
  auto res = lb_walberla()->get_node_velocity(ind);
  return res;
}

REGISTER_CALLBACK_ONE_RANK(get_node_velocity);


boost::optional<double> get_node_density(Utils::Vector3i ind) {
  auto res = lb_walberla()->get_node_density(ind);
  return res;
}

REGISTER_CALLBACK_ONE_RANK(get_node_density);

boost::optional<bool> get_node_is_boundary(Utils::Vector3i ind) {
  auto res = lb_walberla()->get_node_is_boundary(ind);
  return res;
}

REGISTER_CALLBACK_ONE_RANK(get_node_is_boundary);


void set_node_velocity(Utils::Vector3i ind, Utils::Vector3d u) {
  lb_walberla()->set_node_velocity(ind, u);
}

REGISTER_CALLBACK(set_node_velocity);


void set_node_density(Utils::Vector3i ind, double density) {
  lb_walberla()->set_node_density(ind, density);
}

REGISTER_CALLBACK(set_node_density);

} // namespace Walberla
#endif
