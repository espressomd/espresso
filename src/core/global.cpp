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
#include "utils.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "domain_decomposition.hpp"
#include "dpd.hpp"
#include "layered.hpp"
#include "lees_edwards.hpp"
#include "npt.hpp"
#include "tuning.hpp"
#include "rattle.hpp"
#include "imd.hpp"
#include "ghmc.hpp"
#include "lb.hpp"
#include "integrate_sd.hpp"

/** This array contains the description of all global variables.

    Please declare where the variables come from.
*/
const Datafield fields[] = {
  {box_l,            TYPE_DOUBLE, 3, "box_l",             1 },         /* 0  from grid.cpp */
  {dd.cell_grid,        TYPE_INT, 3, "cell_grid",         6 },         /* 1  from cells.cpp */
  {dd.cell_size,     TYPE_DOUBLE, 3, "cell_size",         6 },         /* 2  from cells.cpp */
  {&dpd_gamma,       TYPE_DOUBLE, 1, "dpd_gamma",         5 },         /* 3  from thermostat.cpp */
  {&dpd_r_cut,       TYPE_DOUBLE, 1, "dpd_r_cut",         5 },         /* 4  from thermostat.cpp */
  {&langevin_gamma,  TYPE_DOUBLE, 1, "gamma",             1 },         /* 5  from thermostat.cpp */
  {&lees_edwards_offset, TYPE_DOUBLE, 1, "lees_edwards_offset", 2 },   /* 6  from lees_edwards.cpp */
  {&integ_switch,       TYPE_INT, 1, "integ_switch",      1 },         /* 7  from integrate.cpp */
  {local_box_l,      TYPE_DOUBLE, 3, "local_box_l",       2 },         /* 8  from global.cpp */
  {&max_cut,         TYPE_DOUBLE, 1, "max_cut",           7 },         /* 9  from interaction_data.cpp */
  {&max_num_cells,      TYPE_INT, 1, "max_num_cells",     5 },         /* 10 from cells.cpp */
  {&max_seen_particle,  TYPE_INT, 1, "max_part",          5 },         /* 11 from particle_data.cpp */
  {&max_range,       TYPE_DOUBLE, 1, "max_range",         5 },         /* 12 from integrate.cpp */
  {&max_skin,        TYPE_DOUBLE, 1, "max_skin",          5 },         /* 13 from integrate.cpp */
  {&min_num_cells,      TYPE_INT, 1, "min_num_cells",     5 },         /* 14  from cells.cpp */
  {&n_layers,           TYPE_INT, 1, "n_layers",          3 },         /* 15 from layered.cpp */
  {&n_nodes,            TYPE_INT, 1, "n_nodes",           3 },         /* 16 from communication.cpp */
  {&n_part,             TYPE_INT, 1, "n_part",            6 },         /* 17 from particle.cpp */
  {&n_particle_types,   TYPE_INT, 1, "n_part_types",      8 },         /* 18 from interaction_data.cpp */
  {&n_rigidbonds,       TYPE_INT, 1, "n_rigidbonds",      5 },         /* 19 from rattle.cpp */
  {node_grid,           TYPE_INT, 3, "node_grid",         2 },         /* 20 from grid.cpp */
  {&nptiso_gamma0,   TYPE_DOUBLE, 1, "nptiso_gamma0",    13 },         /* 21 from thermostat.cpp */
  {&nptiso_gammav,   TYPE_DOUBLE, 1, "nptiso_gammav",    13 },         /* 22 from thermostat.cpp */
  {&nptiso.p_ext,    TYPE_DOUBLE, 1, "npt_p_ext",         7 },         /* 23 from pressure.cpp */
  {&nptiso.p_inst,   TYPE_DOUBLE, 1, "npt_p_inst",       10 },         /* 24 from pressure.cpp */
  {&nptiso.p_inst_av,TYPE_DOUBLE, 1, "npt_p_inst_av",    10 },         /* 25 from pressure.cpp */
  {&nptiso.p_diff,   TYPE_DOUBLE, 1, "npt_p_diff",        7 },         /* 26 from pressure.cpp */
  {&nptiso.piston,   TYPE_DOUBLE, 1, "npt_piston",        6 },         /* 27 from pressure.cpp */
  {&periodic,          TYPE_BOOL, 3, "periodicity",       1 },         /* 28 from grid.cpp */
  {&skin,            TYPE_DOUBLE, 1, "skin",              2 },         /* 29 from integrate.cpp */
  {&temperature,     TYPE_DOUBLE, 1, "temperature",       2 },         /* 30 from thermostat.cpp */
  {&thermo_switch,      TYPE_INT, 1, "thermo_switch",     2 },         /* 31 from thermostat.cpp */
  {&sim_time,        TYPE_DOUBLE, 1, "time",              4 },         /* 32 from integrate.cpp */
  {&time_step,       TYPE_DOUBLE, 1, "time_step",         5 },         /* 33 from integrate.cpp */
  {&timing_samples,     TYPE_INT, 1, "timings",           4 },         /* 34 from tuning.cpp */
  {&transfer_rate,      TYPE_INT, 1, "transfer_rate",     2 },         /* 35 from imd.cpp */
  {&max_cut_nonbonded,TYPE_DOUBLE, 1, "max_cut_nonbonded",9 },         /* 36 from interaction_data.cpp */
  {&verlet_reuse,    TYPE_DOUBLE, 1, "verlet_reuse",      8 },         /* 37 from integrate.cpp */
  {&lattice_switch,     TYPE_INT, 1, "lattice_switch",    2 },         /* 38 from lattice.cpp */
  {&dpd_tgamma,      TYPE_DOUBLE, 1, "dpd_tgamma",        6 },         /* 39 from thermostat.cpp */
  {&dpd_tr_cut,      TYPE_DOUBLE, 1, "dpd_tr_cut",        6 },         /* 40 from thermostat.cpp */
  {&dpd_twf,            TYPE_INT, 1, "dpd_twf",           6 },         /* 41 from thermostat.cpp */
  {&dpd_wf,             TYPE_INT, 1, "dpd_wf",            5 },         /* 42 from thermostat.cpp */
  {&max_cut_bonded,  TYPE_DOUBLE, 1, "max_cut_bonded",    9 },         /* 43 from interaction_data.cpp */
  {&transfer_rate,      TYPE_INT, 1, "vmd_transfer_rate", 5 },         /* 44 from imd_tcl.cpp */
  {&min_global_cut,  TYPE_DOUBLE, 1, "min_global_cut",    5 },         /* 45 from interaction_data.cpp */
  {&ghmc_nmd,           TYPE_INT, 1, "ghmc_nmd",          6 },         /* 46 from thermostat.cpp */
  {&ghmc_phi,        TYPE_DOUBLE, 1, "ghmc_phi",          6 },         /* 47 from thermostat.cpp */
  {&ghmc_mc_res,        TYPE_INT, 1, "ghmc_mc_res",       7 },         /* 48 from ghmc.cpp */
  {&ghmc_mflip,         TYPE_INT, 1, "ghmc_mflip",        7 },         /* 49 from ghmc.cpp */
  {&ghmc_tscale,        TYPE_INT, 1, "ghmc_tscale",       6 },         /* 50 from ghmc.cpp */
  {&lb_components,      TYPE_INT, 1, "lb_components",     2 },         /* 51 from ghmc.cpp */
  {&warnings,           TYPE_INT, 1, "warnings",          1 },         /* 52 from global.cpp */
  {&dpd_ignore_fixed_particles, TYPE_INT, 1, "dpd_ignore_fixed_particles", 1 },         /* 53 from global.cpp */
  {&sd_viscosity,            TYPE_DOUBLE, 1, "sd_viscosoity",              4 },         /* 54 from integrate_sd.cpp */
  {&sd_radius,               TYPE_DOUBLE, 1, "sd_radius",                  4 },         /* 55 from integrate_sd.cpp */
  {&sd_seed,                    TYPE_INT, 2, "sd_seed",                    4 },         /* 56 from integrate_sd.cpp */
  {&sd_random_state,            TYPE_INT, 2, "sd_random_state",            4 },         /* 57 from integrate_sd.cpp */
  {&sd_random_precision,     TYPE_DOUBLE, 1, "sd_precision_random",        4 },         /* 58 from integrate_sd.cpp */
  {&smaller_time_step,TYPE_DOUBLE,1, "smaller_time_step", 5 },         /* 59 from integrate.cpp */
  {configtemp,       TYPE_DOUBLE, 2, "configtemp",        1 },         /* 60 from integrate.cpp */
  { NULL, 0, 0, NULL, 0 }
};

int warnings = 1;
