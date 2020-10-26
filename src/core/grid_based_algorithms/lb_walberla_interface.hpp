#ifndef LB_WALBERLA_INTERFACE_HPP
#define LB_WALBERLA_INTERFACE_HPP

#include "config.hpp"

#ifdef LB_WALBERLA

#include "boost/optional.hpp"

#include <vector>

namespace Walberla {

boost::optional<Utils::Vector3d> get_node_velocity(Utils::Vector3i ind);
boost::optional<Utils::Vector3d>
get_node_last_applied_force(Utils::Vector3i ind);
boost::optional<double> get_node_density(Utils::Vector3i ind);
boost::optional<bool> get_node_is_boundary(Utils::Vector3i ind);
void create_vtk(unsigned delta_N, unsigned initial_count,
                unsigned flag_observables, std::string const &identifier,
                std::string const &base_folder, std::string const &prefix);
void write_vtk(std::string const &vtk_uid);
void switch_vtk(std::string const &vtk_uid, int status);
boost::optional<std::vector<double>> get_node_pop(Utils::Vector3i ind);
boost::optional<Utils::Vector6d> get_node_pressure_tensor(Utils::Vector3i ind);

void set_node_last_applied_force(Utils::Vector3i ind, Utils::Vector3d f);
void set_node_velocity(Utils::Vector3i ind, Utils::Vector3d u);
void set_ext_force_density(Utils::Vector3d f);
void set_node_density(Utils::Vector3i ind, double density);
void set_node_pop(Utils::Vector3i ind, std::vector<double> pop);

Utils::Vector3d get_momentum();

boost::optional<Utils::Vector3d> get_velocity_at_pos(Utils::Vector3d pos);

void add_force_at_pos(Utils::Vector3d pos, Utils::Vector3d f);

} // namespace Walberla
#endif
#endif
