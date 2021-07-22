#include "config.hpp"

#ifdef EK_WALBERLA
#include "ekin_walberla_interface.hpp"

#include "boost/optional.hpp"
#include "ekin_walberla_instance.hpp"
#include "utils/Vector.hpp"

#include "MpiCallbacks.hpp"

namespace walberla {

double ek_get_diffusion() { return ekin_walberla()->get_diffusion(); }
void ek_set_diffusion(double diffusion) {
  ekin_walberla()->set_diffusion(diffusion);
}

REGISTER_CALLBACK(ek_set_diffusion)

double ek_get_kT() { return ekin_walberla()->get_kT(); }
void ek_set_kT(double kT) { ekin_walberla()->set_kT(kT); }

REGISTER_CALLBACK(ek_set_kT)

double ek_get_tau() { return ek_walberla_params()->get_tau(); }

boost::optional<double> ek_get_node_density(const Utils::Vector3i &ind) {
  auto res = ekin_walberla()->get_node_density(ind);
  return res;
}

REGISTER_CALLBACK_ONE_RANK(ek_get_node_density)

void ek_set_node_density(const Utils::Vector3i &ind, double density) {
  ekin_walberla()->set_node_density(ind, density);
  ekin_walberla()->ghost_communication();
}

REGISTER_CALLBACK(ek_set_node_density)

boost::optional<bool> ek_get_node_is_boundary(const Utils::Vector3i &ind) {
  return ekin_walberla()->get_node_is_boundary(ind);
}

REGISTER_CALLBACK_ONE_RANK(ek_get_node_is_boundary)

Utils::Vector3i ek_get_shape() {
  return ekin_walberla()->get_blockforest()->get_grid_dimensions();
}

void ek_create_vtk(unsigned delta_N, unsigned initial_count,
                   unsigned flag_observables, std::string const &identifier,
                   std::string const &base_folder, std::string const &prefix) {
  ekin_walberla()->create_vtk(delta_N, initial_count, flag_observables,
                              identifier, base_folder, prefix);
}

REGISTER_CALLBACK(ek_create_vtk)

void ek_write_vtk(std::string const &vtk_uid) {
  ekin_walberla()->write_vtk(vtk_uid);
}

REGISTER_CALLBACK(ek_write_vtk)

void ek_switch_vtk(std::string const &vtk_uid, int status) {
  ekin_walberla()->switch_vtk(vtk_uid, status);
}

REGISTER_CALLBACK(ek_switch_vtk)

void ek_propagate() { ekin_walberla()->integrate(); }

} // namespace walberla

#endif // EK_WALBERLA
