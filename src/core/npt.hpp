/*
 * Copyright (C) 2010-2022 The ESPResSo project
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
/** \file
 *  Exports for the NpT code.
 */

#pragma once

#include "config/config.hpp"

#ifdef ESPRESSO_NPT

#include <utils/Vector.hpp>

#include <array>
#include <cstddef>
#include <vector>

namespace System {
class System;
} // namespace System

/** Parameters of the isotropic NpT-integration scheme. */
struct NptIsoParameters {
  NptIsoParameters() = default;
  NptIsoParameters(double ext_pressure, double piston,
                   Utils::Vector<bool, 3> const &rescale, bool cubic_box);
  /** mass of a virtual piston representing the shaken box */
  double piston = 0.;
  /** inverse of \ref piston */
  double inv_piston = 0.;
  /** isotropic volume. Note that we use the term volume throughout,
   *  although for a 2d or 1d system we mean Area and Length respectively
   */
  double volume = 0.;
  /** list of particle mass */
  std::vector<double> mass_list;
  /** number of particles */
  std::size_t particle_number = 0ul;
  /** desired pressure to which the algorithm strives to */
  double p_ext = 0.;
  /** conjugate momentum of volume */
  double p_epsilon = 0.;
  /** the coefficient of propagation of position for MTK equation */
  double half_dt_inv_piston = 0.;
  /** the coefficient of propagation of momentum for MTK equation */
  double half_dt_inv_piston_and_Nf = 0.;
  /** geometry information for the NpT integrator. Holds the vector
   *  \< dir, dir, dir \> where a positive value for dir indicates that
   *  box movement is allowed in that direction. To check whether a
   *  given direction is turned on, use bitwise comparison with \ref
   *  nptgeom_dir
   */
  int geometry = 0;
  /** The number of dimensions in which NpT boxlength motion is coupled to
   *  particles */
  int dimension = 0;
  /** Set this flag if you want all box dimensions to be identical. Needed for
   *  electrostatics and magnetostatics. If the value of \ref dimension is
   *  less than 3, then box length motion in one or more directions will
   *  be decoupled from the particle motion
   */
  bool cubic_box = false;
  /** An index to one of the non-constant dimensions. Handy if you just want
   *  the variable box_l
   */
  int non_const_dim = -1;
  void coulomb_dipole_sanity_checks(System::System const &system) const;
  Utils::Vector<bool, 3> get_direction() const;
  static constexpr std::array<int, 3> nptgeom_dir = {1, 2, 4};
};

/** Instantaneous pressure during force calculation for NPT integration*/
struct InstantaneousPressure {
  /** instantaneous pressure for p_inst[0] and virial pressure for p_inst[1]
   * the system currently has */
  Utils::Vector2d p_inst = {0., 0.};
  /** virial (short-range) components of \ref p_inst */
  Utils::Vector3d p_vir = {0., 0., 0.};
  /** ideal gas components of \ref p_inst, derived from the velocities */
  Utils::Vector3d p_vel = {0., 0., 0.};
};

#endif // ESPRESSO_NPT
