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

#ifndef NPT_H
#define NPT_H

#include "config.hpp"

#ifdef NPT
#include "BoxGeometry.hpp"

#include <utils/Vector.hpp>

/** Parameters of the isotropic NpT-integration scheme. */
struct NptIsoParameters {
  /** mass of a virtual piston representing the shaken box */
  double piston;
  /** inverse of \ref piston */
  double inv_piston;
  /** isotropic volume. Note that we use the term volume throughout,
   *  although for a 2d or 1d system we mean Area and Length respectively
   */
  double volume;
  /** desired pressure to which the algorithm strives to */
  double p_ext;
  /** instantaneous pressure the system currently has */
  double p_inst;
  /** difference between \ref p_ext and \ref p_inst */
  double p_diff;
  /** virial (short-range) components of \ref p_inst */
  Utils::Vector3d p_vir;
  /** ideal gas components of \ref p_inst, derived from the velocities */
  Utils::Vector3d p_vel;
  /** geometry information for the NpT integrator. Holds the vector
   *  \< dir, dir, dir \> where a positive value for dir indicates that
   *  box movement is allowed in that direction. To check whether a
   *  given direction is turned on, use bitwise comparison with \ref
   *  nptgeom_dir
   */
  int geometry;
  /** bitwise comparison values corresponding to different directions */
  Utils::Vector3i nptgeom_dir;
  /** The number of dimensions in which NpT boxlength motion is coupled to
   *  particles */
  int dimension;
  /** Set this flag if you want all box dimensions to be identical. Needed for
   *  electrostatics and magnetostatics. If the value of \ref dimension is
   *  less than 3, then box length motion in one or more directions will
   *  be decoupled from the particle motion
   */
  bool cubic_box;
  /** An index to one of the non-constant dimensions. Handy if you just want
   *  the variable box_l
   */
  int non_const_dim;
};

extern NptIsoParameters nptiso;

/** @brief NpT initializer.
 *
 *  @param ext_pressure  Reference pressure
 *  @param piston        Piston mass
 *  @param xdir_rescale  Enable box rescaling in the *x*-direction
 *  @param ydir_rescale  Enable box rescaling in the *y*-direction
 *  @param zdir_rescale  Enable box rescaling in the *z*-direction
 *  @param cubic_box     Determines if the volume fluctuations should also
 *                       apply to dimensions which are switched off by the
 *                       above flags and which do not contribute to the
 *                       pressure (3D) or tension (2D, 1D)
 */
void nptiso_init(double ext_pressure, double piston, bool xdir_rescale,
                 bool ydir_rescale, bool zdir_rescale, bool cubic_box);

/** @name NpT geometry bitmasks.
 *  Allowed values for @ref NptIsoParameters::geometry.
 */
/**@{*/
#define NPTGEOM_XDIR 1
#define NPTGEOM_YDIR 2
#define NPTGEOM_ZDIR 4
/**@}*/

/** @brief Synchronizes NpT state such as instantaneous and average pressure
 */
void synchronize_npt_state();
void npt_ensemble_init(const BoxGeometry &box);
void integrator_npt_sanity_checks();
void npt_reset_instantaneous_virials();
void npt_add_virial_contribution(double energy);
void npt_add_virial_contribution(const Utils::Vector3d &force,
                                 const Utils::Vector3d &d);

#endif // NPT
#endif
