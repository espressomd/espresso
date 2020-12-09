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
/** @file
 *  Provide a C-like interface for ScaFaCoS.
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

#include <list>
#include <stdexcept>
#include <string>

namespace Scafacos {

/** Get available ScaFaCoS methods */
std::list<std::string> available_methods() {
  return Scafacos::available_methods();
}

namespace Dipoles {

static ScafacosContextDipoles *scafacos = nullptr;

void add_pair_force(double q1q2, Utils::Vector3d const &d, double dist,
                    Utils::Vector3d &force) {
  if (!scafacos)
    throw std::runtime_error("Scafacos Dipoles not initialized");
  return scafacos->add_pair_force(q1q2, d, dist, force);
}

double pair_energy(double q1q2, double dist) {
  if (!scafacos)
    throw std::runtime_error("Scafacos Dipoles not initialized");
  return scafacos->pair_energy(q1q2, dist);
}

void add_long_range_force() {
  if (!scafacos)
    throw std::runtime_error("Scafacos Dipoles not initialized");
  scafacos->add_long_range_force();
}

double long_range_energy() {
  if (!scafacos)
    throw std::runtime_error("Scafacos Dipoles not initialized");
  return scafacos->long_range_energy();
}

double get_r_cut() {
  if (!scafacos)
    throw std::runtime_error("Scafacos Dipoles not initialized");
  return scafacos->r_cut();
}

void update_system_params() {
  if (!scafacos)
    throw std::runtime_error("Scafacos Dipoles not initialized");
  scafacos->update_system_params();
}

static void set_parameters_worker(const std::string &method,
                                  const std::string &params) {
  delete scafacos;
  scafacos = nullptr;

  scafacos = new ScafacosContextDipoles(method, comm_cart, params);
  if (!scafacos) {
    runtimeErrorMsg() << "Scafacos Dipoles failed to initialize";
    return;
  }

  scafacos->set_dipolar(true);
  scafacos->update_system_params();

  dipole.method = DIPOLAR_SCAFACOS;
  on_coulomb_change();
}

REGISTER_CALLBACK(set_parameters_worker)

} // namespace Dipoles

namespace Coulomb {

static ScafacosContextCoulomb *scafacos = nullptr;

void add_pair_force(double q1q2, Utils::Vector3d const &d, double dist,
                    Utils::Vector3d &force) {
  if (!scafacos)
    throw std::runtime_error("Scafacos Coulomb not initialized");
  return scafacos->add_pair_force(q1q2, d, dist, force);
}

double pair_energy(double q1q2, double dist) {
  if (!scafacos)
    throw std::runtime_error("Scafacos Coulomb not initialized");
  return scafacos->pair_energy(q1q2, dist);
}

void add_long_range_force() {
  if (!scafacos)
    throw std::runtime_error("Scafacos Coulomb not initialized");
  scafacos->add_long_range_force();
}

double long_range_energy() {
  if (!scafacos)
    throw std::runtime_error("Scafacos Coulomb not initialized");
  return scafacos->long_range_energy();
}

double get_r_cut() {
  if (!scafacos)
    throw std::runtime_error("Scafacos Coulomb not initialized");
  return scafacos->r_cut();
}

void update_system_params() {
  if (!scafacos)
    throw std::runtime_error("Scafacos Coulomb not initialized");
  scafacos->update_system_params();
}

static void set_parameters_worker(const std::string &method,
                                  const std::string &params) {
  delete scafacos;
  scafacos = nullptr;

  scafacos = new ScafacosContextCoulomb(method, comm_cart, params);
  if (!scafacos) {
    runtimeErrorMsg() << "Scafacos Coulomb failed to initialize";
    return;
  }

  scafacos->set_dipolar(false);
  scafacos->update_system_params();

  coulomb.method = COULOMB_SCAFACOS;
  on_coulomb_change();
  scafacos->tune();
}

REGISTER_CALLBACK(set_parameters_worker)

void set_r_cut_and_tune_local(double r_cut) {
  scafacos->set_r_cut_and_tune_local(r_cut);
}

} // namespace Coulomb

void free_handle(bool dipolar) {
  if (this_node == 0)
    mpi_call(free_handle, dipolar);
  if (dipolar) {
    delete Dipoles::scafacos;
    Dipoles::scafacos = nullptr;
  } else {
    delete Coulomb::scafacos;
    Coulomb::scafacos = nullptr;
  }
}

REGISTER_CALLBACK(free_handle)

void set_parameters(const std::string &method, const std::string &params,
                    bool dipolar) {
  if (dipolar) {
    mpi_call_all(Dipoles::set_parameters_worker, method, params);
  } else {
    mpi_call_all(Coulomb::set_parameters_worker, method, params);
  }
}

std::string get_method_and_parameters(bool dipolar) {
  if (dipolar) {
    if (!Dipoles::scafacos)
      throw std::runtime_error("Scafacos Dipoles not initialized");
    return Dipoles::scafacos->get_method_and_parameters();
  } else {
    if (!Coulomb::scafacos)
      throw std::runtime_error("Scafacos Coulomb not initialized");
    return Coulomb::scafacos->get_method_and_parameters();
  }
}

} // namespace Scafacos
#endif /* SCAFACOS */
