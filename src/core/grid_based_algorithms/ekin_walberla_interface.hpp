#ifndef ESPRESSO_EKIN_WALBERLA_INTERFACE_HPP
#define ESPRESSO_EKIN_WALBERLA_INTERFACE_HPP
#include "config.hpp"

#ifdef EK_WALBERLA
#include "boost/optional.hpp"
#include "utils/Vector.hpp"

namespace walberla {

boost::optional<double> ek_get_node_density(const Utils::Vector3i &ind);
void ek_set_node_density(const Utils::Vector3i &ind, double density);

boost::optional<bool> ek_get_node_is_boundary(const Utils::Vector3i &ind);

double ek_get_diffusion();
void ek_set_diffusion(double diffusion);

double ek_get_kT();
void ek_set_kT(double kT);

double ek_get_tau();

Utils::Vector3i ek_get_shape();

void ek_create_vtk(unsigned delta_N, unsigned initial_count,
                   unsigned flag_observables, std::string const &identifier,
                   std::string const &base_folder, std::string const &prefix);
void ek_write_vtk(std::string const &vtk_uid);
void ek_switch_vtk(std::string const &vtk_uid, int status);

void ek_propagate();

} // namespace walberla

#endif // EK_WALBERLA
#endif // ESPRESSO_EKIN_WALBERLA_INTERFACE_HPP
