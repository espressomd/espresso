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
class Propagation;
class InteractionsNonBonded;
class BondedInteractionsMap;
namespace Thermostat {
class Thermostat;
}
class ComFixed;
class Galilei;
class Observable_stat;
class OifGlobal;
class ImmersedBoundaries;
namespace CollisionDetection {
class CollisionDetection;
}
namespace BondBreakage {
class BondBreakage;
}
namespace LeesEdwards {
class LeesEdwards;
}
namespace Accumulators {
class AutoUpdateAccumulators;
}
namespace Constraints {
class Constraints;
}

namespace System {

/**
 * @brief Main system class.
 *
 * Most components follow the composite pattern and the opaque pointer pattern.
 * See @ref SystemClassDesign for more details.
 */
class System : public std::enable_shared_from_this<System> {
private:
  struct Private {};
  void initialize();

public:
  System(Private);

  static std::shared_ptr<System> create();

#ifdef CUDA
  GpuParticleData gpu;
#endif
  ResourceCleanup cleanup_queue;

  /** @brief Get @ref time_step. */
  auto get_time_step() const { return time_step; }

  /** @brief Set @ref time_step. */
  void set_time_step(double value);

  /** @brief Get @ref sim_time. */
  auto get_sim_time() const { return sim_time; }

  /** @brief Set @ref sim_time. */
  void set_sim_time(double value);

  /** @brief Get @ref force_cap. */
  auto get_force_cap() const { return force_cap; }

  /** @brief Set @ref force_cap. */
  void set_force_cap(double value);

  /** @brief Get @ref min_global_cut. */
  auto get_min_global_cut() const { return min_global_cut; }

  /** @brief Set @ref min_global_cut. */
  void set_min_global_cut(double value);

  /** @brief Change the box dimensions. */
  void set_box_l(Utils::Vector3d const &box_l);

  /**
   * @brief Tune the Verlet skin.
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

  /** @brief Get the interaction range. */
  double get_interaction_range() const;

  unsigned get_global_ghost_flags() const;

  /** Check electrostatic and magnetostatic methods are properly initialized.
   *  @return true if sanity checks failed.
   */
  bool long_range_interactions_sanity_checks() const;

  /** @brief Calculate the total energy. */
  std::shared_ptr<Observable_stat> calculate_energy();

  /** @brief Calculate the pressure from a virial expansion. */
  std::shared_ptr<Observable_stat> calculate_pressure();

  /** @brief Calculate all forces. */
  void calculate_forces();

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

  /** Integrate equations of motion
   *  @param n_steps       Number of integration steps, can be zero
   *  @param reuse_forces  Decide when to re-calculate forces
   *
   *  @details This function calls two hooks for propagation kernels such as
   *  velocity verlet, velocity verlet + npt box changes, and steepest_descent.
   *  One hook is called before and one after the force calculation.
   *  It is up to the propagation kernels to increment the simulation time.
   *
   *  This function propagates the system according to the choice of integrator
   *  stored in @ref Propagation::integ_switch. The general structure is:
   *  - if reuse_forces is zero, recalculate the forces based on the current
   *    state of the system
   *  - Loop over the number of simulation steps:
   *    -# initialization (e.g., RATTLE)
   *    -# First hook for propagation kernels
   *    -# Update dependent particles and properties (RATTLE, virtual sites)
   *    -# Calculate forces for the current state of the system. This includes
   *       forces added by the Langevin thermostat and the
   *       Lattice-Boltzmann-particle coupling
   *    -# Second hook for propagation kernels
   *    -# Update dependent properties (Virtual sites, RATTLE)
   *    -# Run single step algorithms (Lattice-Boltzmann propagation, collision
   *       detection, NpT update)
   *  - Final update of dependent properties and statistics/counters
   *
   *  High-level documentation of the integration and thermostatting schemes
   *  can be found in doc/sphinx/system_setup.rst and /doc/sphinx/running.rst
   *
   *  @return number of steps that have been integrated, or a negative error
   * code
   */
  int integrate(int n_steps, int reuse_forces);

  int integrate_with_signal_handler(int n_steps, int reuse_forces,
                                    bool update_accumulators);

  /** @brief Calculate initial particle forces from active thermostats. */
  void thermostat_force_init();
  /** @brief Calculate particle-lattice interactions. */
  void lb_couple_particles();

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
  /** called before calculating observables, i.e. energy, pressure or
   *  the integrator (forces). Initialize any methods here which are not
   *  initialized immediately (P3M etc.).
   */
  void on_observable_calc();
  void on_lees_edwards_change();
  void veto_boxl_change(bool skip_particle_checks = false) const;
  /**@}*/

  /**
   * @brief Update particles with properties depending on other particles,
   * namely virtual sites and ICC charges.
   */
  void update_dependent_particles();
  /**
   * @brief Update the global propagation bitmask.
   */
  void update_used_propagations();
  /**
   * @brief Veto temperature change.
   */
  void check_kT(double value) const;

  Coulomb::Solver coulomb;
  Dipoles::Solver dipoles;
  LB::Solver lb;
  EK::Solver ek;
  std::shared_ptr<BoxGeometry> box_geo;
  std::shared_ptr<LocalBox> local_geo;
  std::shared_ptr<CellStructure> cell_structure;
  std::shared_ptr<Propagation> propagation;
  std::shared_ptr<BondedInteractionsMap> bonded_ias;
  std::shared_ptr<InteractionsNonBonded> nonbonded_ias;
  std::shared_ptr<Thermostat::Thermostat> thermostat;
  std::shared_ptr<ComFixed> comfixed;
  std::shared_ptr<Galilei> galilei;
  std::shared_ptr<OifGlobal> oif_global;
  std::shared_ptr<ImmersedBoundaries> immersed_boundaries;
#ifdef COLLISION_DETECTION
  std::shared_ptr<CollisionDetection::CollisionDetection> collision_detection;
#endif
  std::shared_ptr<BondBreakage::BondBreakage> bond_breakage;
  std::shared_ptr<LeesEdwards::LeesEdwards> lees_edwards;
  std::shared_ptr<Accumulators::AutoUpdateAccumulators>
      auto_update_accumulators;
  std::shared_ptr<Constraints::Constraints> constraints;

protected:
  /** @brief Whether the thermostat has to be reinitialized before integration.
   */
  bool reinit_thermo;
  /** @brief Molecular dynamics integrator time step. */
  double time_step;
  /** @brief Molecular dynamics integrator simulation time. */
  double sim_time;
  /** @brief Molecular dynamics integrator force capping. */
  double force_cap;
  /**
   * @brief Minimal global interaction cutoff.
   * Particles with a distance smaller than this are guaranteed
   * to be available on the same node (through ghosts).
   */
  double min_global_cut;

  void update_local_geo();
#ifdef ELECTROSTATICS
  void update_icc_particles();
#endif // ELECTROSTATICS

private:
  /**
   * @brief Check integrator parameters and incompatibilities between
   * the integrator and the currently active thermostat(s).
   */
  void integrator_sanity_checks() const;
};

System &get_system();
void set_system(std::shared_ptr<System> new_instance);
void reset_system();

} // namespace System
