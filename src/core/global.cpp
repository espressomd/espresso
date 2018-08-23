/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
/** \file global.cpp
    Implementation of \ref global.hpp "global.h".
*/
#include "global.hpp"

#include "communication.hpp"
#include "domain_decomposition.hpp"
#include "errorhandling.hpp"
#include "ghmc.hpp"
#include "grid.hpp"
#include "initialize.hpp"
#include "interaction_data.hpp"
#include "lattice.hpp"
#include "layered.hpp"
#include "npt.hpp"
#include "object-in-fluid/oif_global_forces.hpp"
#include "rattle.hpp"
#include "thermalized_bond.hpp"
#include "tuning.hpp"
#include "utils/mpi/all_compare.hpp"

#include <boost/functional/hash.hpp>

#include <functional>
#include <unordered_map>

extern double force_cap;

namespace {

/** Type describing global variables. These are accessible from the
    front end, and are distributed to all compute nodes. */
typedef struct {
  enum class Type { INT = 0, DOUBLE = 1, BOOL = 2 };
  /** Physical address of the variable. */
  void *data;
  /** Type of the variable, either \ref TYPE_INT or \ref TYPE_DOUBLE.*/
  Type type;
  /** Dimension of the variable. Limited to \ref MAX_DIMENSION */
  int dimension;
  /** Name of the variable, mainly used for the front end and debugging */
  const char *name;
} Datafield;

/** This array contains the description of all global variables.

    Please declare where the variables come from.
*/

const std::unordered_map<int, Datafield> fields{
    {FIELD_BOXL,
     {box_l, Datafield::Type::DOUBLE, 3, "box_l"}}, /* 0  from grid.cpp */
    {FIELD_CELLGRID,
     {dd.cell_grid, Datafield::Type::INT, 3,
      "cell_grid"}}, /* 1  from cells.cpp */
    {FIELD_CELLSIZE,
     {dd.cell_size, Datafield::Type::DOUBLE, 3,
      "cell_size"}}, /* 2  from cells.cpp */
#ifndef PARTICLE_ANISOTROPY
    {FIELD_LANGEVIN_GAMMA,
     {&langevin_gamma, Datafield::Type::DOUBLE, 1,
      "gamma"}}, /* 5  from thermostat.cpp */
#else
    {FIELD_LANGEVIN_GAMMA,
     {langevin_gamma.data(), Datafield::Type::DOUBLE, 3,
      "gamma"}}, /* 5  from thermostat.cpp */
#endif // PARTICLE_ANISOTROPY
    {FIELD_INTEG_SWITCH,
     {&integ_switch, Datafield::Type::INT, 1,
      "integ_switch"}}, /* 7  from integrate.cpp */
    {FIELD_LBOXL,
     {local_box_l, Datafield::Type::DOUBLE, 3,
      "local_box_l"}}, /* 8  from global.cpp */
    {FIELD_MCUT,
     {&max_cut, Datafield::Type::DOUBLE, 1,
      "max_cut"}}, /* 9  from interaction_data.cpp */
    {FIELD_MAXNUMCELLS,
     {&max_num_cells, Datafield::Type::INT, 1,
      "max_num_cells"}}, /* 10 from cells.cpp */
    {FIELD_MAXPART,
     {&max_seen_particle, Datafield::Type::INT, 1,
      "max_part"}}, /* 11 from particle_data.cpp */
    {FIELD_MAXRANGE,
     {&max_range, Datafield::Type::DOUBLE, 1,
      "max_range"}}, /* 12 from integrate.cpp */
    {FIELD_MAXSKIN,
     {&max_skin, Datafield::Type::DOUBLE, 1,
      "max_skin"}}, /* 13 from integrate.cpp */
    {FIELD_MINNUMCELLS,
     {&min_num_cells, Datafield::Type::INT, 1,
      "min_num_cells"}}, /* 14  from cells.cpp */
    {FIELD_NLAYERS,
     {&n_layers, Datafield::Type::INT, 1,
      "n_layers"}}, /* 15 from layered.cpp */
    {FIELD_NNODES,
     {&n_nodes, Datafield::Type::INT, 1,
      "n_nodes"}}, /* 16 from communication.cpp */
    {FIELD_NPART,
     {&n_part, Datafield::Type::INT, 1, "n_part"}}, /* 17 from particle.cpp */
    {FIELD_NPARTTYPE,
     {&max_seen_particle_type, Datafield::Type::INT, 1,
      "max_seen_particle_type"}}, /* 18 from interaction_data.cpp */
    {FIELD_RIGIDBONDS,
     {&n_rigidbonds, Datafield::Type::INT, 1,
      "n_rigidbonds"}}, /* 19 from rattle.cpp */
    {FIELD_NODEGRID,
     {node_grid, Datafield::Type::INT, 3, "node_grid"}}, /* 20 from grid.cpp */
    {FIELD_NPTISO_G0,
     {&nptiso_gamma0, Datafield::Type::DOUBLE, 1,
      "nptiso_gamma0"}}, /* 21 from thermostat.cpp */
    {FIELD_NPTISO_GV,
     {&nptiso_gammav, Datafield::Type::DOUBLE, 1,
      "nptiso_gammav"}}, /* 22 from thermostat.cpp */
    {FIELD_NPTISO_PEXT,
     {&nptiso.p_ext, Datafield::Type::DOUBLE, 1,
      "npt_p_ext"}}, /* 23 from pressure.cpp */
    {FIELD_NPTISO_PINST,
     {&nptiso.p_inst, Datafield::Type::DOUBLE, 1,
      "npt_p_inst"}}, /* 24 from pressure.cpp */
    {FIELD_NPTISO_PINSTAV,
     {&nptiso.p_inst_av, Datafield::Type::DOUBLE, 1,
      "npt_p_inst_av"}}, /* 25 from pressure.cpp */
    {FIELD_NPTISO_PDIFF,
     {&nptiso.p_diff, Datafield::Type::DOUBLE, 1,
      "npt_p_diff"}}, /* 26 from pressure.cpp */
    {FIELD_NPTISO_PISTON,
     {&nptiso.piston, Datafield::Type::DOUBLE, 1,
      "npt_piston"}}, /* 27 from pressure.cpp */
    {FIELD_PERIODIC,
     {&periodic, Datafield::Type::BOOL, 3,
      "periodicity"}}, /* 28 from grid.cpp */
    {FIELD_SKIN,
     {&skin, Datafield::Type::DOUBLE, 1, "skin"}}, /* 29 from integrate.cpp */
    {FIELD_TEMPERATURE,
     {&temperature, Datafield::Type::DOUBLE, 1,
      "temperature"}}, /* 30 from thermostat.cpp */
    {FIELD_THERMO_SWITCH,
     {&thermo_switch, Datafield::Type::INT, 1,
      "thermo_switch"}}, /* 31 from thermostat.cpp */
    {FIELD_SIMTIME,
     {&sim_time, Datafield::Type::DOUBLE, 1,
      "time"}}, /* 32 from integrate.cpp */
    {FIELD_TIMESTEP,
     {&time_step, Datafield::Type::DOUBLE, 1,
      "time_step"}}, /* 33 from integrate.cpp */
    {FIELD_TIMINGSAMP,
     {&timing_samples, Datafield::Type::INT, 1,
      "timings"}}, /* 34 from tuning.cpp */
    {FIELD_MCUT_NONBONDED,
     {&max_cut_nonbonded, Datafield::Type::DOUBLE, 1,
      "max_cut_nonbonded"}}, /* 35 from interaction_data.cpp */
    {FIELD_VERLETREUSE,
     {&verlet_reuse, Datafield::Type::DOUBLE, 1,
      "verlet_reuse"}}, /* 36 from integrate.cpp */
    {FIELD_LATTICE_SWITCH,
     {&lattice_switch, Datafield::Type::INT, 1,
      "lattice_switch"}}, /* 37 from lattice.cpp */
    {FIELD_MCUT_BONDED,
     {&max_cut_bonded, Datafield::Type::DOUBLE, 1,
      "max_cut_bonded"}}, /* 42 from interaction_data.cpp */
    {FIELD_MIN_GLOBAL_CUT,
     {&min_global_cut, Datafield::Type::DOUBLE, 1,
      "min_global_cut"}}, /* 43 from interaction_data.cpp */
    {FIELD_GHMC_NMD,
     {&ghmc_nmd, Datafield::Type::INT, 1,
      "ghmc_nmd"}}, /* 44 from thermostat.cpp */
    {FIELD_GHMC_PHI,
     {&ghmc_phi, Datafield::Type::DOUBLE, 1,
      "ghmc_phi"}}, /* 45 from thermostat.cpp */
    {FIELD_GHMC_RES,
     {&ghmc_mc_res, Datafield::Type::INT, 1,
      "ghmc_mc_res"}}, /* 46 from ghmc.cpp */
    {FIELD_GHMC_FLIP,
     {&ghmc_mflip, Datafield::Type::INT, 1,
      "ghmc_mflip"}}, /* 47 from ghmc.cpp */
    {FIELD_GHMC_SCALE,
     {&ghmc_tscale, Datafield::Type::INT, 1,
      "ghmc_tscale"}}, /* 48 from ghmc.cpp */
    {FIELD_WARNINGS,
     {&warnings, Datafield::Type::INT, 1, "warnings"}}, /* 50 from global.cpp */
#ifndef PARTICLE_ANISOTROPY
    {FIELD_LANGEVIN_GAMMA_ROTATION,
     {&langevin_gamma_rotation, Datafield::Type::DOUBLE, 1,
      "gamma_rot"}}, /* 55 from thermostat.cpp */
#else
    {FIELD_LANGEVIN_GAMMA_ROTATION,
     {langevin_gamma_rotation.data(), Datafield::Type::DOUBLE, 3,
      "gamma_rot"}}, /* 55 from thermostat.cpp */
#endif
#ifdef OIF_GLOBAL_FORCES
    {FIELD_MAX_OIF_OBJECTS,
     {&max_oif_objects, Datafield::Type::INT, 1, "max_oif_objects"}},
#endif
    {FIELD_THERMALIZEDBONDS,
     {&n_thermalized_bonds, Datafield::Type::INT, 1,
      "n_thermalized_bonds"}}, /* 56 from thermalized_bond.cpp */
    {FIELD_FORCE_CAP, {&force_cap, Datafield::Type::DOUBLE, 1, "force_cap"}},
    {FIELD_THERMO_VIRTUAL,
     {&thermo_virtual, Datafield::Type::BOOL, 1, "thermo_virtual"}}};

std::size_t hash_value(Datafield const &field) {
  using boost::hash_range;

  switch (field.type) {
  case Datafield::Type::INT: {
    auto ptr = reinterpret_cast<int *>(field.data);
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
  };
}

void common_bcast_parameter(int i) {
  switch (fields.at(i).type) {
  case Datafield::Type::INT:
    MPI_Bcast((int *)fields.at(i).data, fields.at(i).dimension, MPI_INT, 0,
              comm_cart);
    break;
  case Datafield::Type::BOOL:
    static_assert(sizeof(bool) == sizeof(char),
                  "bool datatype does not have the expected size");
    MPI_Bcast((char *)fields.at(i).data, 1, MPI_CHAR, 0, comm_cart);
    break;
  case Datafield::Type::DOUBLE:
    MPI_Bcast((double *)fields.at(i).data, fields.at(i).dimension, MPI_DOUBLE,
              0, comm_cart);
    break;
  default:
    throw std::runtime_error("Unknown type.");
    break;
  }

  on_parameter_change(i);
}
}

int warnings = 1;

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

/*************** REQ_BCAST_PAR ************/

void mpi_bcast_parameter_slave(int node, int i) {
  common_bcast_parameter(i);
  check_runtime_errors();
}

int mpi_bcast_parameter(int i) {
  mpi_call(mpi_bcast_parameter_slave, -1, i);

  common_bcast_parameter(i);

  return check_runtime_errors();
}
