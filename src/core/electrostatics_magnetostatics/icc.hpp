/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
 *
 *  ICC is a method that allows to take into account the influence
 *  of arbitrarily shaped dielectric interfaces. The dielectric
 *  properties of a dielectric medium in the bulk of the simulation
 *  box are taken into account by reproducing the jump in the electric
 *  field at the interface with charge surface segments. The charge
 *  density of the surface segments have to be determined
 *  self-consistently using an iterative scheme. It can at present
 *  be used with P3M, ELCP3M and MMM1D. For details see: @cite tyagi10a
 *
 *  To set up ICC, first the dielectric boundary has to be modeled
 *  by ESPResSo particles n_0...n_0+n where n_0 and n have to be passed
 *  as a parameter to ICC.
 *
 *  For the determination of the induced charges only the forces
 *  acting on the induced charges has to be determined. As P3M and the
 *  other Coulomb solvers calculate all mutual forces, the force
 *  calculation was modified to avoid the calculation of the short
 *  range part of the source-source force calculation. For different
 *  particle data organisation schemes this is performed differently.
 */

#ifndef CORE_ICC_HPP
#define CORE_ICC_HPP

#include "config.hpp"

#if defined(ELECTROSTATICS)

#include "CellStructure.hpp"
#include "ParticleRange.hpp"

#include <utils/Vector.hpp>

#include <vector>

/** ICC data structure */
struct icc_struct {
  /** First id of ICC particle */
  int n_icc;
  /** maximum number of iterations */
  int num_iteration = 30;
  /** bulk dielectric constant */
  double eout;
  /** areas of the particles */
  std::vector<double> areas;
  /** dielectric constants of the particles */
  std::vector<double> ein;
  /** surface charge density of the particles */
  std::vector<double> sigma;
  /** convergence criteria */
  double convergence = 1e-2;
  /** surface normal vectors */
  std::vector<Utils::Vector3d> normals;
  /** external electric field */
  Utils::Vector3d ext_field = {0, 0, 0};
  /** relaxation parameter */
  double relax;
  /** last number of iterations */
  int citeration = 0;
  /** first ICC particle id */
  int first_id = 0;

  template <typename Archive>
  void serialize(Archive &ar, long int /* version */) {
    ar &n_icc;
    ar &num_iteration;
    ar &first_id;
    ar &convergence;
    ar &eout;
    ar &relax;
    ar &areas;
    ar &ein;
    ar &normals;
    ar &sigma;
    ar &ext_field;
    ar &citeration;
  }
};

/** ICC parameters */
extern icc_struct icc_cfg;

/** The main iterative scheme, where the surface element charges are calculated
 *  self-consistently.
 */
void icc_iteration(CellStructure &cell_structure,
                   const ParticleRange &particles,
                   const ParticleRange &ghost_particles);

/** Perform ICC initialization.
 *  @return non-zero value on error
 */
int mpi_icc_init();

/** Set ICC parameters
 */
void icc_set_params(int n_ic, double convergence, double relaxation,
                    Utils::Vector3d const &ext_field, int max_iterations,
                    int first_id, double eps_out, std::vector<double> &areas,
                    std::vector<double> &e_in, std::vector<double> &sigma,
                    std::vector<Utils::Vector3d> &normals);

/** clear ICC vector allocations
 */
void icc_deactivate();

#endif /* ELECTROSTATICS */
#endif /* CORE_ICC_HPP */
