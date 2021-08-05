#include "config.hpp"

#ifdef EK_WALBERLA
#include "ekin_walberla_interface.hpp"

#include "boost/optional.hpp"
#include "ekin_walberla_instance.hpp"
#include "utils/Vector.hpp"

#include "MpiCallbacks.hpp"

namespace walberla {

double ek_get_diffusion(uint id) {
  return get_ek_walberla(id)->get_diffusion();
}
void ek_set_diffusion(uint id, double diffusion) {
  get_ek_walberla(id)->set_diffusion(diffusion);
}

REGISTER_CALLBACK(ek_set_diffusion)

double ek_get_kT(uint id) { return get_ek_walberla(id)->get_kT(); }
void ek_set_kT(uint id, double kT) { get_ek_walberla(id)->set_kT(kT); }

REGISTER_CALLBACK(ek_set_kT)

double ek_get_tau(uint id) { return get_ek_instance_walberla(id).get_tau(); }

boost::optional<double> ek_get_node_density(uint id,
                                            const Utils::Vector3i &ind) {
  auto res = get_ek_walberla(id)->get_node_density(ind);
  return res;
}

REGISTER_CALLBACK_ONE_RANK(ek_get_node_density)

void ek_set_node_density(uint id, const Utils::Vector3i &ind, double density) {
  get_ek_walberla(id)->set_node_density(ind, density);
  get_ek_walberla(id)->ghost_communication();
}

REGISTER_CALLBACK(ek_set_node_density)

boost::optional<bool> ek_get_node_is_boundary(uint id,
                                              const Utils::Vector3i &ind) {
  return get_ek_walberla(id)->get_node_is_boundary(ind);
}

REGISTER_CALLBACK_ONE_RANK(ek_get_node_is_boundary)

Utils::Vector3i ek_get_shape(uint id) {
  return get_ek_walberla(id)->get_blockforest()->get_grid_dimensions();
}

void ek_create_vtk(uint id, unsigned delta_N, unsigned initial_count,
                   unsigned flag_observables, std::string const &identifier,
                   std::string const &base_folder, std::string const &prefix) {
  get_ek_walberla(id)->create_vtk(delta_N, initial_count, flag_observables,
                                  identifier, base_folder, prefix);
}

REGISTER_CALLBACK(ek_create_vtk)

void ek_write_vtk(uint id, std::string const &vtk_uid) {
  get_ek_walberla(id)->write_vtk(vtk_uid);
}

REGISTER_CALLBACK(ek_write_vtk)

void ek_switch_vtk(uint id, std::string const &vtk_uid, int status) {
  get_ek_walberla(id)->switch_vtk(vtk_uid, status);
}

REGISTER_CALLBACK(ek_switch_vtk)

void ek_propagate() {
  std::for_each(get_eks_walberla().begin(), get_eks_walberla().end(),
                [](const EKWalberlaInstance &ek_instance) {
                  ek_instance.get_ek()->integrate();
                });
}
} // namespace walberla

#endif // EK_WALBERLA
