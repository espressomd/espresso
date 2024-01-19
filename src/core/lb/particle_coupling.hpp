/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "BoxGeometry.hpp"
#include "LocalBox.hpp"
#include "Particle.hpp"
#include "ParticleRange.hpp"
#include "PropagationMode.hpp"
#include "cell_system/CellStructure.hpp"
#include "lb/Solver.hpp"
#include "system/System.hpp"
#include "thermostat.hpp"

#include <utils/Vector.hpp>
#include <utils/math/sqr.hpp>

#include <cmath>
#include <unordered_set>
#include <vector>

/**
 * @brief Check if a position is within the local LB domain plus halo.
 *
 * @return True iff the point is inside of the domain.
 */
bool in_local_halo(LocalBox const &local_box, Utils::Vector3d const &pos,
                   double agrid);

/**
 * @brief Add a force to the lattice force density.
 * @param lb Hydrodynamics solver
 * @param pos Position of the force in LB units.
 * @param force Force in MD units.
 * @param time_step MD time step.
 */
void add_md_force(LB::Solver &lb, Utils::Vector3d const &pos,
                  Utils::Vector3d const &force, double time_step);

// internal function exposed for unit testing
std::vector<Utils::Vector3d> positions_in_halo(Utils::Vector3d const &pos,
                                               BoxGeometry const &box,
                                               LocalBox const &local_geo,
                                               double agrid);

/** @brief Calculate drag force on a single particle.
 *
 *  See section II.C. @cite ahlrichs99a
 *
 *  @param[in] lb          The coupled fluid
 *  @param[in] lb_gamma    The friction coefficient
 *  @param[in] p           The coupled particle
 *  @param[in] shifted_pos The particle position in LB units with optional shift
 *  @param[in] vel_offset  Velocity offset in MD units to be added to
 *                         interpolated LB velocity before calculating the force
 *
 *  @return The viscous coupling force
 */
Utils::Vector3d lb_drag_force(LB::Solver const &lb, double lb_gamma,
                              Particle const &p,
                              Utils::Vector3d const &shifted_pos,
                              Utils::Vector3d const &vel_offset);

namespace LB {

class ParticleCoupling {
  LBThermostat const &m_thermostat;
  LB::Solver &m_lb;
  BoxGeometry const &m_box_geo;
  LocalBox const &m_local_box;
  double m_time_step;
  double m_noise_pref_wo_gamma;
  bool m_thermalized;

public:
  ParticleCoupling(LBThermostat const &thermostat, LB::Solver &lb,
                   BoxGeometry const &box_geo, LocalBox const &local_box,
                   double time_step)
      : m_thermostat{thermostat}, m_lb{lb}, m_box_geo{box_geo},
        m_local_box{local_box}, m_time_step{time_step} {
    /* Eq. (16) @cite ahlrichs99a, without the gamma term.
     * The factor 12 comes from the fact that we use random numbers
     * from -0.5 to 0.5 (equally distributed) which have variance 1/12.
     * The time step comes from the discretization.
     */
    auto constexpr variance_inv = 12.;
    auto const kT = lb.get_kT() * Utils::sqr(lb.get_lattice_speed());
    m_thermalized = (kT != 0.);
    m_noise_pref_wo_gamma = std::sqrt(variance_inv * 2. * kT / time_step);
  }

  Utils::Vector3d get_noise_term(Particle const &p) const;
  void kernel(Particle &p);

  /**
   * @brief Calculate particle drift velocity offset due to ENGINE and
   * ELECTROHYDRODYNAMICS.
   */
  auto lb_drift_velocity_offset(Particle const &p) const {
    Utils::Vector3d vel_offset{};
#ifdef LB_ELECTROHYDRODYNAMICS
    vel_offset += p.mu_E();
#endif
    return vel_offset;
  }
};

/**
 * @brief Keep track of ghost particles that have already been coupled once.
 * In certain cases, there may be more than one ghost for the same particle.
 * To make sure these are only coupled once, ghosts' ids are recorded.
 */
class CouplingBookkeeping {
  std::unordered_set<int> m_coupled_ghosts;
  CellStructure const &m_cell_structure;

  /** @brief Check if there is locally a real particle for the given ghost. */
  bool is_ghost_for_local_particle(Particle const &p) const {
    return not m_cell_structure.get_local_particle(p.id())->is_ghost();
  }

public:
  explicit CouplingBookkeeping(CellStructure const &cell_structure)
      : m_cell_structure{cell_structure} {}

  /** @brief Determine if a given particle should be coupled. */
  bool should_be_coupled(Particle const &p) {
    auto const propagation = p.propagation();
    if ((propagation & PropagationMode::TRANS_LB_MOMENTUM_EXCHANGE) == 0 and
        (propagation & PropagationMode::SYSTEM_DEFAULT) == 0 and
        (propagation & PropagationMode::TRANS_LB_TRACER) == 0) {
      return false;
    }
#ifdef VIRTUAL_SITES_RELATIVE
    if ((propagation & PropagationMode::TRANS_LB_MOMENTUM_EXCHANGE) == 0 and
        propagation & (PropagationMode::TRANS_VS_RELATIVE |
                       PropagationMode::ROT_VS_RELATIVE)) {
      return false;
    }
#endif
    // real particles: always couple
    if (not p.is_ghost()) {
      return true;
    }
    // ghosts: check we don't have the corresponding real particle on the same
    // node, and that a ghost for the same particle hasn't been coupled already
    if (m_coupled_ghosts.count(p.id()) != 0 or is_ghost_for_local_particle(p)) {
      return false;
    }
    m_coupled_ghosts.insert(p.id());
    return true;
  }
};

inline bool is_tracer(Particle const &p) {
  return (p.propagation() & PropagationMode::TRANS_LB_TRACER) != 0;
}

} // namespace LB
