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
 * This file contains the c-type wrapper interface to the (oop-) scafacos
 * interface. */

#ifndef ES_CORE_ELECTROSTATICS_MAGNETOSTATICS_SCAFACOS_HPP
#define ES_CORE_ELECTROSTATICS_MAGNETOSTATICS_SCAFACOS_HPP

#include "config.hpp"

#include <utils/Vector.hpp>

#include <list>
#include <string>

namespace Scafacos {
#if defined(SCAFACOS)

namespace Dipoles {
/** Near-field pair force */
void add_pair_force(double q1q2, Utils::Vector3d const &d, double dist,
                    Utils::Vector3d &force);
/** Near-field pair energy */
double pair_energy(double q1q2, double dist);
/** Long range part */
void add_long_range_force();
/** Calculate long range energy contribution */
double long_range_energy();
/** Get parameters */
std::string get_method_and_parameters();
/** Set parameters */
void set_parameters(const std::string &method, const std::string &params);
double get_r_cut();

void update_system_params();
} // namespace Dipoles

namespace Coulomb {
/** Near-field pair force */
void add_pair_force(double q1q2, Utils::Vector3d const &d, double dist,
                    Utils::Vector3d &force);
/** Near-field pair energy */
double pair_energy(double q1q2, double dist);
/** Long range part */
void add_long_range_force();
/** Calculate long range energy contribution */
double long_range_energy();
double get_r_cut();

void update_system_params();

void set_r_cut_and_tune_local(double r_cut);
} // namespace Coulomb

#endif /* SCAFACOS */

std::list<std::string> available_methods();

std::string get_method_and_parameters(bool dipolar);
void set_parameters(const std::string &method, const std::string &params,
                    bool dipolar);
void free_handle(bool dipolar);
} // namespace Scafacos
#endif
