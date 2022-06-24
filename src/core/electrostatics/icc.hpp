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

#ifndef ESPRESSO_SRC_CORE_ELECTROSTATICS_ICC_HPP
#define ESPRESSO_SRC_CORE_ELECTROSTATICS_ICC_HPP

/**
 * @file
 *
 * ICC is a method that allows to take into account the influence
 * of arbitrarily shaped dielectric interfaces. The dielectric
 * properties of a dielectric medium in the bulk of the simulation
 * box are taken into account by reproducing the jump in the electric
 * field at the interface with charge surface segments. The charge
 * density of the surface segments have to be determined
 * self-consistently using an iterative scheme. It can at present
 * be used with P3M, ELCP3M and MMM1D. For details see: @cite tyagi10a
 *
 * To set up ICC, first the dielectric boundary has to be modeled
 * by ESPResSo particles n_0...n_0+n where n_0 and n have to be passed
 * as a parameter to ICC.
 *
 * For the determination of the induced charges, only the forces acting
 * on the induced charges have to be determined. As P3M and the other
 * Coulomb solvers calculate all mutual forces, the force calculation
 * was modified to avoid the calculation of the short-range part
 * of the source-source force calculation. For different particle
 * data organisation schemes, this is performed differently.
 */

#include "config.hpp"

#ifdef ELECTROSTATICS

#include "ParticleRange.hpp"
#include "cell_system/CellStructure.hpp"

#include <utils/Vector.hpp>

#include <vector>

/** ICC data structure */
struct icc_data {
  /** First id of ICC particle */
  int n_icc;
  /** maximum number of iterations */
  int max_iterations;
  /** bulk dielectric constant */
  double eps_out;
  /** areas of the particles */
  std::vector<double> areas;
  /** dielectric constants of the particles */
  std::vector<double> epsilons;
  /** surface charge density of the particles */
  std::vector<double> sigmas;
  /** convergence criteria */
  double convergence;
  /** surface normal vectors */
  std::vector<Utils::Vector3d> normals;
  /** external electric field */
  Utils::Vector3d ext_field;
  /** relaxation parameter */
  double relaxation;
  /** last number of iterations */
  int citeration;
  /** first ICC particle id */
  int first_id;

  void sanity_checks() const;
};

struct ICCStar {
  /** ICC parameters */
  icc_data icc_cfg;

  ICCStar(icc_data data);

  /**
   * The main iterative scheme, where the surface element charges are calculated
   * self-consistently.
   */
  void iteration(CellStructure &cell_structure, ParticleRange const &particles,
                 ParticleRange const &ghost_particles);

  void on_activation() const;
  void sanity_checks_active_solver() const;
  void sanity_check() const;
};

void update_icc_particles();

#endif // ELECTROSTATICS
#endif
