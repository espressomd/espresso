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
 * interface.
 */

#ifndef ES_CORE_ELECTROSTATICS_MAGNETOSTATICS_SCAFACOS_HPP
#define ES_CORE_ELECTROSTATICS_MAGNETOSTATICS_SCAFACOS_HPP

#include "config.hpp"

#if defined(SCAFACOS)

#include "electrostatics_magnetostatics/ScafacosContextBase.hpp"

#include <utils/Vector.hpp>

#include <list>
#include <string>

namespace Scafacos {

/** @brief Access the per-MPI-node ScaFaCoS Coulomb instance */
ScafacosContextBase *fcs_coulomb();
#ifdef SCAFACOS_DIPOLES
/** @brief Access the per-MPI-node ScaFaCoS dipoles instance */
ScafacosContextBase *fcs_dipoles();
#endif

std::string get_method_and_parameters(bool dipolar);
void set_parameters(const std::string &method, const std::string &params,
                    bool dipolar);
void free_handle(bool dipolar);

void set_r_cut_and_tune(double r_cut);

std::list<std::string> available_methods();

} // namespace Scafacos
#endif // SCAFACOS
#endif
