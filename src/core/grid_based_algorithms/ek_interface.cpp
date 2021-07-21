#include "config.hpp"

#include "ek_interface.hpp"

#include "communication.hpp"
#include "ekin_walberla_interface.hpp"
#include "utils/Vector.hpp"

namespace EK {

ActiveEK lattice_switch = ActiveEK::NONE;

struct NoEKActive : public std::exception {
  const char *what() const noexcept override { return "EK not activated"; }
};

namespace detail {
bool node_is_index_valid(const Utils::Vector3i &ind,
                         const Utils::Vector3i &shape) {
  return ind < shape && ind >= Utils::Vector3i::broadcast(0);
}
} // namespace detail

double get_density(const Utils::Vector3i &ind) {
#ifdef EK_WALBERLA
  if (EK::get_lattice_switch() == EK::ActiveEK::WALBERLA) {
    return ::Communication::mpiCallbacks().call(
        ::Communication::Result::one_rank, walberla::ek_get_node_density, ind);
  }
#endif // EK_WALBERLA
  throw NoEKActive();
}

void set_density(const Utils::Vector3i &ind, double density) {
#ifdef EK_WALBERLA
  if (EK::get_lattice_switch() == EK::ActiveEK::WALBERLA) {
    mpi_call_all(walberla::ek_set_node_density, ind, density);
  } else
#endif // EK_WALBERLA
  {
    throw NoEKActive();
  }
}

bool get_is_boundary(const Utils::Vector3i &ind) {
#ifdef EK_WALBERLA
  if (EK::get_lattice_switch() == EK::ActiveEK::WALBERLA) {
    return ::Communication::mpiCallbacks().call(
        ::Communication::Result::one_rank, walberla::ek_get_node_is_boundary,
        ind);
  }
#endif // EK_WALBERLA
  throw NoEKActive();
}

double get_diffusion() {
#ifdef EK_WALBERLA
  if (EK::get_lattice_switch() == EK::ActiveEK::WALBERLA) {
    return walberla::ek_get_diffusion();
  }
#endif // EK_WALBERLA
  throw NoEKActive();
}

void set_diffusion(double diffusion) {
#ifdef EK_WALBERLA
  if (EK::get_lattice_switch() == EK::ActiveEK::WALBERLA) {
    mpi_call_all(walberla::ek_set_diffusion, diffusion);
  } else
#endif // EK_WALBERLA
  {
    throw NoEKActive();
  }
}

double get_kT() {
#ifdef EK_WALBERLA
  if (EK::get_lattice_switch() == EK::ActiveEK::WALBERLA) {
    return walberla::ek_get_kT();
  }
#endif // EK_WALBERLA
  throw NoEKActive();
}

void set_kT(double kT) {
#ifdef EK_WALBERLA
  if (EK::get_lattice_switch() == EK::ActiveEK::WALBERLA) {
    mpi_call_all(walberla::ek_set_kT, kT);
  } else
#endif // EK_WALBERLA
  {
    throw NoEKActive();
  }
}

double get_tau() {
#ifdef EK_WALBERLA
  if (EK::get_lattice_switch() == EK::ActiveEK::WALBERLA) {
    return walberla::ek_get_tau();
  }
#endif // EK_WALBERLA
  throw NoEKActive();
}

Utils::Vector3i get_shape() {
#ifdef EK_WALBERLA
  if (EK::get_lattice_switch() == EK::ActiveEK::WALBERLA) {
    return walberla::ek_get_shape();
  }
#endif // EK_WALBERLA
  throw NoEKActive();
}

bool node_is_index_valid(const Utils::Vector3i &ind) {
  return detail::node_is_index_valid(ind, get_shape());
}

void create_vtk(unsigned delta_N, unsigned initial_count,
                unsigned flag_observables, std::string const &identifier,
                std::string const &base_folder, std::string const &prefix) {
#ifdef LB_WALBERLA
  if (EK::get_lattice_switch() == EK::ActiveEK::WALBERLA) {
    mpi_call_all(walberla::ek_create_vtk, delta_N, initial_count,
                 flag_observables, identifier, base_folder, prefix);
    return;
  }
#endif
  throw NoEKActive();
}

EK::ActiveEK get_lattice_switch() { return EK::lattice_switch; }

void propagate() {
#ifdef EK_WALBERLA
  if (EK::get_lattice_switch() == EK::ActiveEK::WALBERLA) {
    walberla::ek_propagate();
  }
#endif
}

} // namespace EK

void mpi_set_ek_lattice_switch_local(EK::ActiveEK lattice_switch) {
  ::EK::lattice_switch = lattice_switch;
}

REGISTER_CALLBACK(mpi_set_ek_lattice_switch_local)

void mpi_set_ek_lattice_switch(EK::ActiveEK lattice_switch) {
  mpi_call_all(mpi_set_ek_lattice_switch_local, lattice_switch);
}
