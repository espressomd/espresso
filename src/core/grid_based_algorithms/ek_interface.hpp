#ifndef ESPRESSO_EK_INTERFACE_HPP
#define ESPRESSO_EK_INTERFACE_HPP

#include "utils/Vector.hpp"

namespace EK {

/** @brief EK implementation currently active. */
enum class ActiveEK : int { NONE, WALBERLA };

double get_density(const Utils::Vector3i &ind);
void set_density(const Utils::Vector3i &ind, double density);

bool get_is_boundary(const Utils::Vector3i &ind);

double get_diffusion();
void set_diffusion(double diffusion);

double get_kT();
void set_kT(double kT);

double get_tau();

Utils::Vector3i get_shape();

bool node_is_index_valid(const Utils::Vector3i &ind);

/**
 * @brief Create a VTK observable.
 */
void create_vtk(unsigned delta_N, unsigned initial_count,
                unsigned flag_observables, std::string const &identifier,
                std::string const &base_folder, std::string const &prefix);

EK::ActiveEK get_lattice_switch();

void propagate();
} // namespace EK

void mpi_set_ek_lattice_switch(EK::ActiveEK lattice_switch);

#endif // ESPRESSO_EK_INTERFACE_HPP
