/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
//

/** \file

    ICCP3M is a method that allows to take into account the influence
    of arbitrarily shaped dielectric interfaces.  The dielectric
    properties of a dielectric medium in the bulk of the simulation
    box are taken into account by reproducing the jump in the electric
    field at the inface with charge surface segments. The charge
    density of the surface segments have to be determined
    self-consistently using an iterative scheme.  It can at presently
    - despite its name - be used with P3M, ELCP3M, MMM2D and MMM1D.
    For details see:<br> S. Tyagi, M. Suzen, M. Sega, C. Holm,
    M. Barbosa: A linear-scaling method for computing induced charges
    on arbitrary dielectric boundaries in large system simulations
    (Preprint)

    To set up ICCP3M first the dielectric boundary has to be modelled
    by espresso particles 0..n where n has to be passed as a parameter
    to ICCP3M. This is still a bit inconvenient, as it forces the user
    to reserve the first n particle ids to wall charges, but as the
    other parts of espresso do not suffer from a limitation like this,
    it can be tolerated.

    For the determination of the induced charges only the forces
    acting on the induced charges has to be determined. As P3M and the
    other Coulomb solvers calculate all mutual forces, the force
    calculation was modified to avoid the calculation of the short
    range part of the source-source force calculation.  For different
    particle data organisation schemes this is performed differently.
    */

#ifndef CORE_ICCP3M_HPP
#define CORE_ICCP3M_HPP

#include "config.hpp"

#if defined(ELECTROSTATICS)

#include <utils/Vector.hpp>

/* iccp3m data structures*/
struct iccp3m_struct {
  int n_ic;                  /* Last induced id (can not be smaller then 2) */
  int num_iteration = 30;    /* Number of max iterations                    */
  double eout = 1;           /* Dielectric constant of the bulk             */
  std::vector<double> areas; /* Array of area of the grid elements          */
  std::vector<double>
      ein; /* Array of dielectric constants at each surface element */
  std::vector<double> sigma; /* Surface Charge density */
  double convergence = 1e-2; /* Convergence criterion                       */
  std::vector<Utils::Vector3d> normals;  /* Surface normal vectors */
  Utils::Vector3d ext_field = {0, 0, 0}; /* External field */
  double relax = 0.7; /* relaxation parameter for iterative */
  int citeration = 0; /* current number of iterations*/
  int first_id = 0;

  template <typename Archive>
  void serialize(Archive &ar, long int /* version */) {
    ar &n_ic;
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
extern iccp3m_struct iccp3m_cfg; /* global variable with ICCP3M configuration */

/** The main iterative scheme, where the surface element charges are calculated
 * self-consistently.
 */
int iccp3m_iteration();

/** The allocation of ICCP3M lists for python interface
 */
void iccp3m_alloc_lists();

/** check sanity of parameters for use with ICCP3M
 */
int iccp3m_sanity_check();

#endif /* ELECTROSTATICS */

#endif /* ICCP3M_H */
