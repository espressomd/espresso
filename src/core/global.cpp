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
 *  Implementation of \ref global.hpp "global.hpp".
 */
#include "global.hpp"

#include "bonded_interactions/rigid_bond.hpp"
#include "bonded_interactions/thermalized_bond.hpp"
#include "communication.hpp"
#include "config.hpp"
#include "errorhandling.hpp"
#include "event.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "integrate.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "object-in-fluid/oif_global_forces.hpp"
#include "rattle.hpp"
#include "thermostat.hpp"

#include <utils/mpi/all_compare.hpp>

#include <boost/functional/hash.hpp>
#include <boost/mpi/collectives.hpp>

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <utility>

extern double force_cap;

namespace {

/** Type describing global variables. These are accessible from the
 *  front end, and are distributed to all compute nodes.
 */
struct Datafield {
  enum class Type { INT = 0, DOUBLE = 1, BOOL = 2, UNSIGNED_LONG = 3 };

  Datafield(int *data, size_t dimension, const char *name)
      : data(data), type(Type::INT), dimension(dimension), name(name) {}

  Datafield(double *data, size_t dimension, const char *name)
      : data(data), type(Type::DOUBLE), dimension(dimension), name(name) {}

  Datafield(bool *data, size_t dimension, const char *name)
      : data(data), type(Type::BOOL), dimension(dimension), name(name) {}

  Datafield(unsigned long *data, size_t dimension, const char *name)
      : data(data), type(Type::UNSIGNED_LONG), dimension(dimension),
        name(name) {}

  /** Physical address of the variable. */
  void *data;
  /** Type of the variable. */
  Type type;
  /** Dimension of the variable. Typically in the range 1-3. */
  size_t dimension;
  /** Name of the variable, mainly used for the front end and debugging */
  const char *name;
};

/** This array contains the description of all global variables.
 *
 *  Please declare where the variables come from.
 */
const std::unordered_map<int, Datafield> fields{
#ifndef PARTICLE_ANISOTROPY
    {FIELD_LANGEVIN_GAMMA,
     {&langevin.gamma, 1, "langevin.gamma"}}, /* 5  from thermostat.cpp */
#else
    {FIELD_LANGEVIN_GAMMA,
     {langevin.gamma.data(), 3, "langevin.gamma"}}, /* 5  from thermostat.cpp */
#endif // PARTICLE_ANISOTROPY
    {FIELD_INTEG_SWITCH,
     {&integ_switch, 1, "integ_switch"}}, /* 7  from integrate.cpp */
    {FIELD_RIGIDBONDS,
     {&n_rigidbonds, 1, "n_rigidbonds"}}, /* 19 from rigid_bond.cpp */
    {FIELD_NODEGRID, {node_grid.data(), 3, "node_grid"}}, /* 20 from grid.cpp */
#ifdef NPT
    {FIELD_NPTISO_G0,
     {&npt_iso.gamma0, 1, "npt_iso.gamma0"}}, /* 21 from thermostat.cpp */
    {FIELD_NPTISO_GV,
     {&npt_iso.gammav, 1, "npt_iso.gammav"}}, /* 22 from thermostat.cpp */
#endif
    {FIELD_SKIN, {&skin, 1, "skin"}}, /* 29 from integrate.cpp */
    {FIELD_TEMPERATURE,
     {&temperature, 1, "temperature"}}, /* 30 from thermostat.cpp */
    {FIELD_THERMO_SWITCH,
     {&thermo_switch, 1, "thermo_switch"}},         /* 31 from thermostat.cpp */
    {FIELD_SIMTIME, {&sim_time, 1, "time"}},        /* 32 from integrate.cpp */
    {FIELD_TIMESTEP, {&time_step, 1, "time_step"}}, /* 33 from integrate.cpp */
    {FIELD_LATTICE_SWITCH,
     {reinterpret_cast<std::underlying_type_t<ActiveLB> *>(&lattice_switch), 1,
      "lattice_switch"}}, /* 37 from lattice.cpp */
    {FIELD_MIN_GLOBAL_CUT,
     {&min_global_cut, 1, "min_global_cut"}}, /* 43 from interaction_data.cpp */
#ifndef PARTICLE_ANISOTROPY
    {FIELD_LANGEVIN_GAMMA_ROTATION,
     {&langevin.gamma_rotation, 1,
      "langevin.gamma_rotation"}}, /* 55 from thermostat.cpp */
#else
    {FIELD_LANGEVIN_GAMMA_ROTATION,
     {langevin.gamma_rotation.data(), 3,
      "langevin.gamma_rotation"}}, /* 55 from thermostat.cpp */
#endif
    {FIELD_MAX_OIF_OBJECTS, {&max_oif_objects, 1, "max_oif_objects"}},
    {FIELD_THERMALIZEDBONDS,
     {&n_thermalized_bonds, 1,
      "n_thermalized_bonds"}}, /* 56 from thermalized_bond.cpp */
    {FIELD_FORCE_CAP, {&force_cap, 1, "force_cap"}},
    {FIELD_THERMO_VIRTUAL, {&thermo_virtual, 1, "thermo_virtual"}},
#ifndef PARTICLE_ANISOTROPY
    {FIELD_BROWNIAN_GAMMA_ROTATION,
     {&brownian.gamma_rotation, 1,
      "brownian.gamma_rotation"}}, /* 57 from thermostat.cpp */
#else
    {FIELD_BROWNIAN_GAMMA_ROTATION,
     {brownian.gamma_rotation.data(), 3,
      "brownian.gamma_rotation"}}, /* 57 from thermostat.cpp */
#endif
#ifndef PARTICLE_ANISOTROPY
    {FIELD_BROWNIAN_GAMMA,
     {&brownian.gamma, 1, "brownian.gamma"}}, /* 58  from thermostat.cpp */
#else
    {FIELD_BROWNIAN_GAMMA,
     {brownian.gamma.data(), 3,
      "brownian.gamma"}}, /* 58  from thermostat.cpp */
#endif // PARTICLE_ANISOTROPY
};

std::size_t hash_value(Datafield const &field) {
  using boost::hash_range;

  switch (field.type) {
  case Datafield::Type::INT: {
    auto ptr = reinterpret_cast<int *>(field.data);
    return hash_range(ptr, ptr + field.dimension);
  }
  case Datafield::Type::UNSIGNED_LONG: {
    auto ptr = reinterpret_cast<unsigned long *>(field.data);
    return hash_range(ptr, ptr + field.dimension);
  }
  case Datafield::Type::BOOL: {
    auto ptr = reinterpret_cast<char *>(field.data);
    return hash_range(ptr, ptr + 1);
  }
  case Datafield::Type::DOUBLE: {
    auto ptr = reinterpret_cast<double *>(field.data);
    return hash_range(ptr, ptr + field.dimension);
  }
  default:
    throw std::runtime_error("Unknown type.");
  }
}

void common_bcast_parameter(int i) {
  switch (fields.at(i).type) {
  case Datafield::Type::INT:
    boost::mpi::broadcast(comm_cart, reinterpret_cast<int *>(fields.at(i).data),
                          fields.at(i).dimension, 0);
    break;
  case Datafield::Type::UNSIGNED_LONG:
    boost::mpi::broadcast(comm_cart,
                          reinterpret_cast<unsigned long *>(fields.at(i).data),
                          fields.at(i).dimension, 0);
    break;
  case Datafield::Type::BOOL:
    static_assert(sizeof(bool) == sizeof(char),
                  "bool datatype does not have the expected size");
    boost::mpi::broadcast(comm_cart,
                          reinterpret_cast<char *>(fields.at(i).data), 1, 0);
    break;
  case Datafield::Type::DOUBLE:
    boost::mpi::broadcast(comm_cart,
                          reinterpret_cast<double *>(fields.at(i).data),
                          fields.at(i).dimension, 0);
    break;
  default:
    throw std::runtime_error("Unknown type.");
  }

  on_parameter_change(i);
}
} // namespace

void check_global_consistency() {
  using Utils::Mpi::all_compare;

  /* Local hash */
  auto const hash = boost::hash_range(fields.begin(), fields.end());

  bool is_same = all_compare(comm_cart, hash);

  /* If the hash does not check out, test the individual fields to find out
     which ones are async. */
  if (not is_same) {
    for (auto const &field : fields) {
      if (not all_compare(comm_cart, hash_value(field.second))) {
        runtimeErrorMsg() << "Global field '" << field.second.name << "' ("
                          << field.first
                          << ") is not synchronized between all nodes.";
      }
    }
  }
}

/*************** BCAST PARAMETER ************/

REGISTER_CALLBACK(common_bcast_parameter)

void mpi_bcast_parameter(int i) { mpi_call_all(common_bcast_parameter, i); }
