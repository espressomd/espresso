/*
 * Copyright (C) 2014-2022 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "config/config.hpp"

#include "GpuParticleData.hpp"
#include "ResourceCleanup.hpp"

#include "electrostatics/solver.hpp"
#include "magnetostatics/solver.hpp"

#include "ek/Solver.hpp"
#include "lb/Solver.hpp"

#include "cell_system/CellStructureType.hpp"

#include <utils/Vector.hpp>

#include <memory>

class BoxGeometry;
class LocalBox;
struct CellStructure;
class InteractionsNonBonded;
class ComFixed;
class Galilei;
class Observable_stat;

namespace System {

class System {
public:
  System();
#ifdef CUDA
  GpuParticleData gpu;
#endif
  ResourceCleanup cleanup_queue;

  Utils::Vector3d box() const;
  void init();

  /** @brief Get @ref time_step. */
  auto get_time_step() { return time_step; }

  /** @brief Set @ref time_step. */
  void set_time_step(double value);

  /** @brief Get @ref force_cap. */
  auto get_force_cap() { return force_cap; }

  /** @brief Set @ref force_cap. */
  void set_force_cap(double value);

  /** @brief Get @ref min_global_cut. */
  auto get_min_global_cut() const { return min_global_cut; }

  /** @brief Set @ref min_global_cut. */
  void set_min_global_cut(double value);

  /** @brief Get @ref verlet_skin. */
  auto get_verlet_skin() const { return verlet_skin; }

  /** @brief Get @ref verlet_skin. */
  auto is_verlet_skin_set() const { return verlet_skin_set; }

  /** @brief Set @ref verlet_skin. */
  void set_verlet_skin(double value);

  /** @brief Get @ref verlet_reuse. */
  auto get_verlet_reuse() const { return verlet_reuse; }

  /** @brief Update @ref verlet_reuse. */
  void update_verlet_stats(int n_steps, int n_verlet_updates);

  /**
   * @brief Tune @ref verlet_skin.
   * Choose the optimal Verlet list skin between @p min_skin and @p max_skin
   * by bisection to tolerance @p tol.
   */
  void tune_verlet_skin(double min_skin, double max_skin, double tol,
                        int int_steps, bool adjust_max_skin);

  /** @brief Change cell structure topology. */
  void set_cell_structure_topology(CellStructureType topology);

  /** @brief Rebuild cell lists. Use e.g. after a skin change. */
  void rebuild_cell_structure();

  /** @brief Calculate the maximal cutoff of all interactions. */
  double maximal_cutoff() const;

  double get_interaction_range() const;

  /** Check electrostatic and magnetostatic methods are properly initialized.
   *  @return true if sanity checks failed.
   */
  bool long_range_interactions_sanity_checks() const;

  /** @brief Calculate the total energy. */
  std::shared_ptr<Observable_stat> calculate_energy();

  /** @brief Calculate the pressure from a virial expansion. */
  std::shared_ptr<Observable_stat> calculate_pressure();

#ifdef DIPOLE_FIELD_TRACKING
  /** @brief Calculate dipole fields. */
  void calculate_long_range_fields();
#endif

  /**
   * @brief Compute the short-range energy of a particle.
   *
   * Iterate through particles inside cell and neighboring cells and compute
   * energy contribution for a specific particle.
   *
   * @param pid    Particle id
   * @return Non-bonded energy of the particle.
   */
  double particle_short_range_energy_contribution(int pid);

  /** \name Hook procedures
   *  These procedures are called if several significant changes to
   *  the system happen which may make a reinitialization of subsystems
   *  necessary.
   */
  /**@{*/
  /**
   * @brief Called when the box length has changed. This routine is relatively
   * fast, and changing the box length every time step as for example necessary
   * for NpT is more or less ok.
   *
   * @param skip_method_adaption skip the long-range methods adaptions
   */
  void on_boxl_change(bool skip_method_adaption = false);
  void on_node_grid_change();
  void on_periodicity_change();
  void on_cell_structure_change();
  void on_thermostat_param_change();
  void on_temperature_change();
  void on_verlet_skin_change();
  void on_timestep_change();
  void on_integration_start();
  void on_short_range_ia_change();
  void on_non_bonded_ia_change();
  void on_coulomb_change();
  void on_dipoles_change();
  /** @brief Called every time a constraint is changed. */
  void on_constraint_change();
  /** @brief Called when the LB boundary conditions change
   *  (geometry, slip velocity, or both).
   */
  void on_lb_boundary_conditions_change();
  /** @brief Called every time a particle local property changes. */
  void on_particle_local_change();
  /** @brief Called every time a particle property changes. */
  void on_particle_change();
  /** @brief Called every time a particle charge changes. */
  void on_particle_charge_change();
  /** @brief Update particles with properties depending on other particles,
   *  namely virtual sites and ICC charges.
   */
  void update_dependent_particles();
  /** called before calculating observables, i.e. energy, pressure or
   *  the integrator (forces). Initialize any methods here which are not
   *  initialized immediately (P3M etc.).
   */
  void on_observable_calc();
  /**@}*/

  Coulomb::Solver coulomb;
  Dipoles::Solver dipoles;
  LB::Solver lb;
  EK::Solver ek;
  std::shared_ptr<BoxGeometry> box_geo;
  std::shared_ptr<LocalBox> local_geo;
  std::shared_ptr<CellStructure> cell_structure;
  std::shared_ptr<InteractionsNonBonded> nonbonded_ias;
  std::shared_ptr<ComFixed> comfixed;
  std::shared_ptr<Galilei> galilei;

protected:
  /** @brief Whether the thermostat has to be reinitialized before integration.
   */
  bool reinit_thermo;
  /** @brief Molecular dynamics integrator time step. */
  double time_step;
  /** @brief Molecular dynamics integrator force capping. */
  double force_cap;
  /**
   * @brief Minimal global interaction cutoff.
   * Particles with a distance smaller than this are guaranteed
   * to be available on the same node (through ghosts).
   */
  double min_global_cut;
  /** @brief Verlet list skin. */
  double verlet_skin;
  /** @brief Average number of integration steps the Verlet list was re-used. */
  double verlet_reuse;
  /** @brief Whether the Verlet list skin was set. */
  bool verlet_skin_set;

  void update_local_geo();
};

System &get_system();
void set_system(std::shared_ptr<System> new_instance);
void reset_system();
bool is_system_set();

} // namespace System
