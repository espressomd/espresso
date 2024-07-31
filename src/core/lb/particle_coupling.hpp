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

// internal function exposed for unit testing
std::vector<Utils::Vector3d> positions_in_halo(Utils::Vector3d const &pos,
                                               BoxGeometry const &box_geo,
                                               LocalBox const &local_box,
                                               double agrid);

/** @brief Calculate drag force on a single particle.
 *
 *  See section II.C. and eq. 9 in @cite ahlrichs99a.
 *
 *  @param[in] lb          The coupled fluid
 *  @param[in] lb_gamma    The friction coefficient
 *  @param[in] p           The coupled particle
 *  @param[in] shifted_pos The particle position in MD units with optional shift
 *
 *  @return The viscous coupling force
 */
Utils::Vector3d lb_drag_force(LB::Solver const &lb, double lb_gamma,
                              Particle const &p,
                              Utils::Vector3d const &shifted_pos);

namespace LB {

class ParticleCoupling {
  LBThermostat const &m_thermostat;
  LB::Solver &m_lb;
  BoxGeometry const &m_box_geo;
  LocalBox const &m_local_box;
  double m_noise_pref_wo_gamma;
  bool m_thermalized;

public:
  ParticleCoupling(LBThermostat const &thermostat, LB::Solver &lb,
                   BoxGeometry const &box_geo, LocalBox const &local_box)
      : m_thermostat{thermostat}, m_lb{lb}, m_box_geo{box_geo},
        m_local_box{local_box} {
    /* Eq. (16) @cite ahlrichs99a, without the gamma term.
     * The factor 12 comes from the fact that we use random numbers
     * from -0.5 to 0.5 (equally distributed) which have variance 1/12.
     * The time step comes from the discretization.
     */
    auto constexpr variance_inv = 12.;
    auto const kT = lb.get_kT() * Utils::sqr(lb.get_lattice_speed());
    m_thermalized = (kT != 0.);
    m_noise_pref_wo_gamma = std::sqrt(variance_inv * 2. * kT / lb.get_tau());
  }

  Utils::Vector3d get_noise_term(Particle const &p) const;
  void kernel(std::vector<Particle *> const &particles);
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
