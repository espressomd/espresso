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
 *
 *  Implementation of \ref gay_berne.hpp
 */
#include "gay_berne.hpp"

#ifdef GAY_BERNE
#include "nonbonded_interaction_data.hpp"

#include <cmath>

GayBerne_Parameters::GayBerne_Parameters(double eps, double sig, double cut,
                                         double k1, double k2, double mu,
                                         double nu)
    : eps{eps}, sig{sig}, cut{cut}, k1{k1}, k2{k2}, mu{mu}, nu{nu},
      chi1{((k1 * k1) - 1.) / ((k1 * k1) + 1.)},
      chi2{(std::pow(k2, 1. / mu) - 1.) / (std::pow(k2, 1. / mu) + 1.)} {}

#endif // GAY_BERNE
