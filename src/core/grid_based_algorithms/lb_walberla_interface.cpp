#include "config.hpp"

#ifdef LB_WALBERLA

#include "MpiCallbacks.hpp"
#include "lb_interface.hpp"
#include "lb_walberla_instance.hpp"
#include "utils/Vector.hpp"

#include "boost/optional.hpp"

namespace Walberla {

boost::optional<Utils::Vector3d> get_node_velocity(Utils::Vector3i ind) {
  auto res = lb_walberla()->get_node_velocity(ind);
  return res;
}

REGISTER_CALLBACK_ONE_RANK(get_node_velocity)

boost::optional<double> get_node_density(Utils::Vector3i ind) {
  auto res = lb_walberla()->get_node_density(ind);
  return res;
}

REGISTER_CALLBACK_ONE_RANK(get_node_density)

boost::optional<bool> get_node_is_boundary(Utils::Vector3i ind) {
  auto res = lb_walberla()->get_node_is_boundary(ind);
  return res;
}

REGISTER_CALLBACK_ONE_RANK(get_node_is_boundary)

boost::optional<Utils::Vector19d> get_node_pop(Utils::Vector3i ind) {
  auto res = lb_walberla()->get_node_pop(ind);
  return res;
}

REGISTER_CALLBACK_ONE_RANK(get_node_pop)

boost::optional<Utils::Vector6d> get_node_pressure_tensor(Utils::Vector3i ind) {
  auto res = lb_walberla()->get_node_pressure_tensor(ind);
  return res;
}

REGISTER_CALLBACK_ONE_RANK(get_node_pressure_tensor)

void set_node_velocity(Utils::Vector3i ind, Utils::Vector3d u) {
  lb_walberla()->set_node_velocity(ind, u);
}

REGISTER_CALLBACK(set_node_velocity)

void set_ext_force_density(Utils::Vector3d f) {
  lb_walberla()->set_external_force(f);
}
REGISTER_CALLBACK(set_ext_force_density)

void set_node_density(Utils::Vector3i ind, double density) {
  lb_walberla()->set_node_density(ind, density);
}

REGISTER_CALLBACK(set_node_density)

void set_node_pop(Utils::Vector3i ind, Utils::Vector19d pop) {
  lb_walberla()->set_node_pop(ind, pop);
}

REGISTER_CALLBACK(set_node_pop)

Utils::Vector3d get_momentum() { return lb_walberla()->get_momentum(); }

REGISTER_CALLBACK_REDUCTION(get_momentum, std::plus<>())

boost::optional<Utils::Vector3d> get_velocity_at_pos(Utils::Vector3d pos) {
  auto res = lb_walberla()->get_velocity_at_pos(pos);
  return res;
}

REGISTER_CALLBACK_ONE_RANK(get_velocity_at_pos)

} // namespace Walberla
#endif
