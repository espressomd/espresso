#ifndef ESPRESSO_EK_INTERFACE_HPP
#define ESPRESSO_EK_INTERFACE_HPP

#include "utils/Vector.hpp"

namespace EK {

/** @brief EK implementation currently active. */
enum class ActiveEK : int { NONE, WALBERLA };

double get_density(uint id, const Utils::Vector3i &ind);
void set_density(uint id, const Utils::Vector3i &ind, double density);

bool get_is_boundary(uint id, const Utils::Vector3i &ind);

double get_diffusion(uint id);
void set_diffusion(uint id, double diffusion);

double get_kT(uint id);
void set_kT(uint id, double kT);

double get_tau(uint id);

Utils::Vector3i get_shape(uint id);

bool node_is_index_valid(uint id, const Utils::Vector3i &ind);

/**
 * @brief Create a VTK observable.
 */
void create_vtk(uint id, unsigned delta_N, unsigned initial_count,
                unsigned flag_observables, std::string const &identifier,
                std::string const &base_folder, std::string const &prefix);

/**
 * @brief Write a VTK observable to disk.
 */
void write_vtk(uint id, std::string const &vtk_uid);

/**
 * @brief Toggle a VTK observable on/off.
 */
void switch_vtk(uint id, std::string const &vtk_uid, int status);

EK::ActiveEK get_lattice_switch();

void propagate();
} // namespace EK

void mpi_set_ek_lattice_switch(uint id, EK::ActiveEK lattice_switch);

#endif // ESPRESSO_EK_INTERFACE_HPP
