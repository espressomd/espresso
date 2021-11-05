/*
 * Copyright (C) 2010-2021 The ESPResSo project
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
#ifndef RIGID_BOND_HPP
#define RIGID_BOND_HPP
/** \file
 *  Definition of the rigid bond data type. It is utilized by the
 *  Rattle algorithm.
 *
 *  Implementation in \ref rigid_bond.cpp.
 */

#include <utils/Vector.hpp>

#include <cmath>

/** Number of rigid bonds. */
extern int n_rigidbonds;

/** Parameters for the rigid_bond/SHAKE/RATTLE ALGORITHM */
struct RigidBond {
  /** Square of the length of Constrained Bond */
  double d2;
  /** Positional Tolerance/Accuracy value for termination of RATTLE/SHAKE
   *  iterations during position corrections
   */
  double p_tol;
  /** Velocity Tolerance/Accuracy for termination of RATTLE/SHAKE iterations
   *  during velocity corrections
   */
  double v_tol;

  double cutoff() const { return std::sqrt(d2); }

  static constexpr int num = 1;

  RigidBond(double d, double p_tol, double v_tol);

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, long int /* version */) {
    ar &d2;
    ar &p_tol;
    ar &v_tol;
  }
};

#endif
