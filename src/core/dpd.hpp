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

#pragma once

/** \file
 *  Routines to use DPD as thermostat or pair force @cite soddemann03a
 *
 *  Implementation in @ref dpd.cpp.
 */

#include <config/config.hpp>

#ifdef ESPRESSO_DPD

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "system/System.hpp"
#include "thermostat.hpp"

#include <utils/Vector.hpp>

#include <cmath>

// Forward declaration
namespace boost::mpi {
class communicator;
}

inline Utils::Vector3d dpd_pair_force(DPDParameters const &params,
                                      Utils::Vector3d const &v, double dist,
                                      Utils::Vector3d const &noise) {
  if (dist < params.cutoff) {
    auto const r_cut = params.cutoff;
    auto const omega = params.wf ? 1. - std::pow((dist / r_cut), params.k) : 1.;
    auto const f_d = params.gamma * (omega * omega) * v;
    auto const f_r = params.pref * omega * noise;
    return f_r - f_d;
  }
  return {};
}

Utils::Vector3d
dpd_pair_force(Utils::Vector3d const &p1_position,
               Utils::Vector3d const &p1_velocity, int const &p1_id,
               Utils::Vector3d const &p2_position,
               Utils::Vector3d const &p2_velocity, int const &p2_id,
               DPDThermostat const &dpd, BoxGeometry const &box_geo,
               IA_parameters const &ia_params, Utils::Vector3d const &d,
               double dist, double dist2);

Utils::Vector9d dpd_stress(System::System &system,
                           boost::mpi::communicator const &comm);

/**
 * @brief Local contribution to the pressure tensor.
 * Needs to be rescaled by the box volume.
 */
Utils::Vector9d dpd_pressure_local(System::System &system);

#endif // ESPRESSO_DPD
