#include "ekin_walberla_interface.hpp"

#include "boost/optional.hpp"
#include "ekin_walberla_instance.hpp"
#include "utils/Vector.hpp"

#include "MpiCallbacks.hpp"

namespace EK {
ActiveEK lattice_switch = ActiveEK::NONE;
}

void ek_propagate() {
  if (EK::lattice_switch == EK::ActiveEK::WALBERLA) {
    ekin_walberla()->integrate();
  }
}

namespace walberla {
namespace detail {
bool node_is_index_valid(const Utils::Vector3i &ind,
                         const Utils::Vector3i &shape) {
  return ind < shape && ind >= Utils::Vector3i::broadcast(0);
}
} // namespace detail

double ek_get_diffusion() { return ekin_walberla()->get_diffusion(); }
void ek_set_diffusion(double diffusion) {
  ekin_walberla()->set_diffusion(diffusion);
}

double ek_get_kT() { return ekin_walberla()->get_kT(); }
Utils::Vector3i ek_get_shape() {
  return ekin_walberla()->get_blockforest()->get_grid_dimensions();
}
void ek_set_kT(double kT) { ekin_walberla()->set_kT(kT); }

double ek_get_tau() { return ek_walberla_params()->get_tau(); }

bool ek_node_is_index_valid(const Utils::Vector3i &ind) {
  return detail::node_is_index_valid(ind, ek_get_shape());
}

boost::optional<double> ek_get_node_density(const Utils::Vector3i &ind) {
  auto res = ekin_walberla()->get_node_density(ind);
  return res;
}

REGISTER_CALLBACK_ONE_RANK(ek_get_node_density)

double ek_get_density(const Utils::Vector3i &ind) {
  return ::Communication::mpiCallbacks().call(::Communication::Result::one_rank,
                                              ek_get_node_density, ind);
}

boost::optional<bool> get_node_is_boundary(Utils::Vector3i ind) {
  return ekin_walberla()->get_node_is_boundary(ind);
}

REGISTER_CALLBACK_ONE_RANK(get_node_is_boundary)

void ek_set_node_density(const Utils::Vector3i &ind, double density) {
  ekin_walberla()->set_node_density(ind, density);
  ekin_walberla()->ghost_communication();
}

REGISTER_CALLBACK(ek_set_node_density)

bool ek_get_node_is_boundary(const Utils::Vector3i &ind) {
  return ::Communication::mpiCallbacks().call(::Communication::Result::one_rank,
                                              get_node_is_boundary, ind);
}
} // namespace walberla

void mpi_set_ek_lattice_switch_local(EK::ActiveEK lattice_switch) {
  ::EK::lattice_switch = lattice_switch;
}

REGISTER_CALLBACK(mpi_set_ek_lattice_switch_local)

void mpi_set_ek_lattice_switch(EK::ActiveEK lattice_switch) {
  mpi_call_all(mpi_set_ek_lattice_switch_local, lattice_switch);
}

EK::ActiveEK ek_get_lattice_switch() { return EK::lattice_switch; }