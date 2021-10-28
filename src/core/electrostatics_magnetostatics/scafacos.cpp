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

#include "config.hpp"

#if defined(SCAFACOS)

#include "electrostatics_magnetostatics/scafacos.hpp"

#include "Scafacos.hpp"
#include "communication.hpp"
#include "electrostatics_magnetostatics/ScafacosContext.hpp"
#include "electrostatics_magnetostatics/coulomb.hpp"
#include "electrostatics_magnetostatics/dipole.hpp"
#include "errorhandling.hpp"
#include "event.hpp"
#include "grid.hpp"

#include <utils/Vector.hpp>
#include <utils/matrix.hpp>

#include <list>
#include <stdexcept>
#include <string>

namespace Scafacos {

/** Get available ScaFaCoS methods */
std::list<std::string> available_methods() {
  return Scafacos::available_methods();
}

namespace {
#ifdef SCAFACOS_DIPOLES
ScafacosContextDipoles *dipoles_instance = nullptr;
#endif
ScafacosContextCoulomb *coulomb_instance = nullptr;
} // namespace

#ifdef SCAFACOS_DIPOLES
ScafacosContextBase *fcs_dipoles() {
  if (!dipoles_instance) {
    throw std::runtime_error(
        "Attempted access to uninitialized Scafacos Dipoles instance.");
  }
  return dipoles_instance;
}
#endif

ScafacosContextBase *fcs_coulomb() {
  if (!coulomb_instance) {
    throw std::runtime_error(
        "Attempted access to uninitialized Scafacos Coulomb instance.");
  }
  return coulomb_instance;
}

#ifdef SCAFACOS_DIPOLES
static void set_parameters_dipoles_local(const std::string &method,
                                         const std::string &params) {
  delete dipoles_instance;
  dipoles_instance = nullptr;

  auto *instance = new ScafacosContextDipoles(method, comm_cart, params);
  if (!instance) {
    runtimeErrorMsg() << "Scafacos Dipoles failed to initialize";
    return;
  }
  dipoles_instance = instance;

  instance->set_dipolar(true);
  instance->update_system_params();

  dipole.method = DIPOLAR_SCAFACOS;
  on_coulomb_change();
}

REGISTER_CALLBACK(set_parameters_dipoles_local)
#endif

static void set_parameters_coulomb_local(const std::string &method,
                                         const std::string &params) {
  delete coulomb_instance;
  coulomb_instance = nullptr;

  auto *instance = new ScafacosContextCoulomb(method, comm_cart, params);
  if (!instance) {
    runtimeErrorMsg() << "Scafacos Coulomb failed to initialize";
    return;
  }
  coulomb_instance = instance;

  instance->set_dipolar(false);
  instance->update_system_params();

  coulomb.method = COULOMB_SCAFACOS;
  on_coulomb_change();
  instance->tune();
}

REGISTER_CALLBACK(set_parameters_coulomb_local)

void set_r_cut_and_tune(double r_cut) {
  coulomb_instance->set_r_cut_and_tune(r_cut);
}

void free_handle(bool dipolar) {
  if (this_node == 0)
    mpi_call(free_handle, dipolar);
  if (dipolar) {
#ifdef SCAFACOS_DIPOLES
    delete dipoles_instance;
    dipoles_instance = nullptr;
#endif
  } else {
    delete coulomb_instance;
    coulomb_instance = nullptr;
  }
}

REGISTER_CALLBACK(free_handle)

void set_parameters(const std::string &method, const std::string &params,
                    bool dipolar) {
  if (dipolar) {
#ifdef SCAFACOS_DIPOLES
    mpi_call_all(set_parameters_dipoles_local, method, params);
#endif
  } else {
    mpi_call_all(set_parameters_coulomb_local, method, params);
  }
}

std::string get_method_and_parameters(bool dipolar) {
  if (dipolar) {
#ifdef SCAFACOS_DIPOLES
    return fcs_dipoles()->get_method_and_parameters();
#else
    return std::string();
#endif
  }
  return fcs_coulomb()->get_method_and_parameters();
}

} // namespace Scafacos
#endif /* SCAFACOS */
