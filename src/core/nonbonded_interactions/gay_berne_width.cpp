/*
 * Copyright (C) 2010-2026 The ESPResSo project
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/** \file
 *
 * Implementation of \ref gay_berne_width.hpp
 */

#include "gay_berne_width.hpp"

#ifdef ESPRESSO_GAY_BERNE_WIDTH

#include "nonbonded_interaction_data.hpp"

#include <cmath>

GayBerneWidth_Parameters::GayBerneWidth_Parameters(double eps, double sig,
                                                   double wid, double cut,
                                                   double k1, double k2,
                                                   double mu, double nu)
    : eps{eps}, sig{sig}, wid{wid}, cut{cut}, k1{k1}, k2{k2}, mu{mu}, nu{nu},
      chi1{((k1 * k1) - 1.) / ((k1 * k1) + 1.)},
      chi2{(std::pow(k2, 1. / mu) - 1.) / (std::pow(k2, 1. / mu) + 1.)} {}

#endif // ESPRESSO_GAY_BERNE_WIDTH
