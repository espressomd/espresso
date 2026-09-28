/*
 * Copyright (C) 2010-2026 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

// Force initialization and Langevin thermostat live in their own translation
// unit. The [[gnu::flatten]] init tree is expensive for gcc to inline; keeping
// it out of forces.cpp leaves the whole-TU inline-growth budget of forces.cpp
// for the hot pair kernel (forces_cabana.hpp) and the Verlet-list build.

#include "forces_init.hpp"

#include <config/config.hpp>

#include "Particle.hpp"
#include "PropagationMode.hpp"
#include "cell_system/CellStructure.hpp"
#include "cell_system/for_each_particle.hpp"
#include "integrators/Propagation.hpp"
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
#include "magnetostatics/dipoles.hpp"
#endif
#include "rotation.hpp"
#include "system/System.hpp"
#include "thermostat.hpp"
#include "thermostats/langevin_inline.hpp"

#ifdef ESPRESSO_CALIPER
#include "caliper_utils.hpp"

#include <caliper/cali.h>
#endif

/** External particle forces */
static ParticleForce external_force(Particle const &p) {
  ParticleForce f = {};

#ifdef ESPRESSO_EXTERNAL_FORCES
  f.f += p.ext_force();
#ifdef ESPRESSO_ROTATION
  f.torque += p.ext_torque();
#endif
#endif

#ifdef ESPRESSO_ENGINE
  // apply a swimming force in the direction of
  // the particle's orientation axis
  if (p.swimming().swimming and !p.swimming().is_engine_force_on_fluid) {
    f.f += p.swimming().f_swim * p.calc_director();
  }
#endif

  return f;
}

/** Combined force initialization and Langevin noise application */
// [[gnu::flatten]]: the per-particle lambda below exceeds gcc's bottom-up
// inlining budget, which turns trivial helpers (Utils::hadamard_product,
// Utils::Vector constructors) into per-particle PLT calls inside the Langevin
// friction path. Flattening forces the whole call tree inline; same
// arithmetic, bitwise-identical trajectories.
[[gnu::flatten]] void init_forces_and_thermostat(System::System const &system) {
#ifdef ESPRESSO_CALIPER
  ESPRESSO_CALI_MARK_FUNCTION;
#endif

  auto &cell_structure = *system.cell_structure;
  auto const &propagation = *system.propagation;
  auto const &thermostat = *system.thermostat;
  auto const kT = thermostat.kT;
  auto const time_step = system.get_time_step();

  // Check if Langevin thermostat is active
  bool const langevin_active =
      thermostat.langevin &&
      (propagation.used_propagations &
       (PropagationMode::TRANS_LANGEVIN | PropagationMode::ROT_LANGEVIN));

  // Hoist the Langevin column-view handles ONCE outside the parallel_for. The
  // Langevin friction is a column kernel (velocity / omega / id (+ gamma,
  // quaternion) read by row); external-force init and the body->space torque
  // rotation stay on the view path (external_force reads the engine sidecar;
  // the rotation needs the quaternion).
  auto &store = cell_structure.particle_store();
  auto vel_view = store.velocity_view();
  auto id_view = store.id_view();
#ifdef ESPRESSO_ROTATION
  auto omega_view = store.angular_velocity_view();
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
  auto gamma_view = store.gamma_view();
#ifdef ESPRESSO_ROTATION
  auto gamma_rot_view = store.gamma_rot_view();
#endif
#endif
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
  auto quat_view = store.quaternion_view();
#endif

  // Single pass over all local particles
  cell_structure.for_each_local_particle_row([&](int const row) {
    Particle p;
    p.attach_to_store(store, row);
    // Initialize force with external forces
    auto const external = external_force(p);
    auto force = p.force();
    force = external.f;
#ifdef ESPRESSO_ROTATION
    auto torque = p.torque();
    torque = external.torque;
#endif

    // Apply Langevin noise if thermostat is active
    if (langevin_active) {
      auto const &langevin = *thermostat.langevin;
      // Read the propagation bitfield once per particle (each
      // should_propagate_with call would otherwise re-read the store column).
      int const prop = p.propagation();
      if (propagation.should_propagate_with(prop,
                                            PropagationMode::TRANS_LANGEVIN))
        force += friction_thermo_langevin(langevin, vel_view, id_view,
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
                                          gamma_view,
#endif
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
                                          quat_view,
#endif
                                          row, time_step, kT);
#ifdef ESPRESSO_ROTATION
      if (propagation.should_propagate_with(prop,
                                            PropagationMode::ROT_LANGEVIN))
        torque += convert_vector_body_to_space(
            p, friction_thermo_langevin_rotation(langevin, omega_view, id_view,
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
                                                 gamma_rot_view,
#endif
                                                 row, time_step, kT));
#endif
    }
  });
  // The local force/torque scatter buffers were already zeroed this force
  // call by prepare_verlet_list_cabana() (via update_verlet_state).

  // Initialize ghost forces
  cell_structure.ghosts_reset_forces();
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  if (system.dipoles.impl->solver.has_value()) {
    cell_structure.ghosts_reset_dipole_fields();
  }
#endif
}
