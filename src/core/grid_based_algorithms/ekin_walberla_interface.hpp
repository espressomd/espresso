#ifndef ESPRESSO_EKIN_WALBERLA_INTERFACE_HPP
#define ESPRESSO_EKIN_WALBERLA_INTERFACE_HPP
#include "config.hpp"

#ifdef EK_WALBERLA
#include "boost/optional.hpp"
#include "utils/Vector.hpp"

namespace EK {
/** @brief EK implementation currently active. */
enum class ActiveEK : int { NONE, WALBERLA };
} // namespace EK

void ek_propagate();

namespace walberla {

double ek_get_density(const Utils::Vector3i &ind);
void ek_set_node_density(const Utils::Vector3i &ind, double density);

bool ek_get_node_is_boundary(const Utils::Vector3i &ind);

double ek_get_diffusion();
void ek_set_diffusion(double diffusion);

double ek_get_kT();
void ek_set_kT(double kT);

double ek_get_tau();

Utils::Vector3i ek_get_shape();

bool ek_node_is_index_valid(const Utils::Vector3i &ind);

} // namespace walberla

void mpi_set_ek_lattice_switch(EK::ActiveEK lattice_switch);

EK::ActiveEK ek_get_lattice_switch();

#endif // EK_WALBERLA
#endif // ESPRESSO_EKIN_WALBERLA_INTERFACE_HPP
