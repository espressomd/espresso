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
// unit. Keeping the init tree out of forces.cpp leaves the whole-TU
// inline-growth budget of forces.cpp for the hot pair kernel
// (forces_cabana.hpp) and the Verlet-list build.

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
void init_forces_and_thermostat(System::System const &system) {
#ifdef ESPRESSO_CALIPER
  CALI_CXX_MARK_FUNCTION;
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

  // Single pass over all local particles
  cell_structure.for_each_local_particle([&](Particle &p) {
    // Initialize force with external forces
    p.force_and_torque() = external_force(p);

    // Apply Langevin noise if thermostat is active
    if (langevin_active) {
      auto const &langevin = *thermostat.langevin;
      if (propagation.should_propagate_with(p, PropagationMode::TRANS_LANGEVIN))
        p.force() += friction_thermo_langevin(langevin, p, time_step, kT);
#ifdef ESPRESSO_ROTATION
      if (propagation.should_propagate_with(p, PropagationMode::ROT_LANGEVIN))
        p.torque() += convert_vector_body_to_space(
            p, friction_thermo_langevin_rotation(langevin, p, time_step, kT));
#endif
    }
  });
  cell_structure.reset_local_force_and_torque();

  // Initialize ghost forces
  cell_structure.ghosts_reset_forces();
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  if (system.dipoles.impl->solver.has_value()) {
    cell_structure.ghosts_reset_dipole_field();
  }
#endif
}
