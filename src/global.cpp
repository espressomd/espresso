/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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
/** \file global.c
    Implementation of \ref global.h "global.h".
*/
#include "utils.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "domain_decomposition.hpp"
#include "dpd.hpp"
#include "layered.hpp"
#include "npt.hpp"
#include "tuning.hpp"
#include "adresso.hpp"
#include "rattle.hpp"
#include "imd.hpp"
#include "ghmc.hpp"
#include "lb.hpp"

/** This array contains the description of all global variables.

    Please declare where the variables come from.
*/
const Datafield fields[] = {
  {box_l,            TYPE_DOUBLE, 3, "box_l",             1 },         /* 0  from grid.c */
  {dd.cell_grid,        TYPE_INT, 3, "cell_grid",         6 },         /* 1  from cells.c */
  {dd.cell_size,     TYPE_DOUBLE, 3, "cell_size",         6 },         /* 2  from cells.c */
  {&dpd_gamma,       TYPE_DOUBLE, 1, "dpd_gamma",         5 },         /* 3  from thermostat.c */
  {&dpd_r_cut,       TYPE_DOUBLE, 1, "dpd_r_cut",         5 },         /* 4  from thermostat.c */
  {&langevin_gamma,  TYPE_DOUBLE, 1, "gamma",             1 },         /* 5  from thermostat.c */
  {&integ_switch,       TYPE_INT, 1, "integ_switch",      1 },         /* 6  from integrate.c */
  {local_box_l,      TYPE_DOUBLE, 3, "local_box_l",       2 },         /* 7  from global.c */
  {&max_cut,         TYPE_DOUBLE, 1, "max_cut",           7 },         /* 8  from interaction_data.c */
  {&max_num_cells,      TYPE_INT, 1, "max_num_cells",     5 },         /* 9 from cells.c */
  {&max_seen_particle,  TYPE_INT, 1, "max_part",          5 },         /* 10 from particle_data.c */
  {&max_range,       TYPE_DOUBLE, 1, "max_range",         5 },         /* 11 from integrate.c */
  {&max_skin,        TYPE_DOUBLE, 1, "max_skin",          5 },         /* 12 from integrate.c */
  {&min_num_cells,      TYPE_INT, 1, "min_num_cells",     5 },         /* 13  from cells.c */
  {&n_layers,           TYPE_INT, 1, "n_layers",          3 },         /* 14 from layered.c */
  {&n_nodes,            TYPE_INT, 1, "n_nodes",           3 },         /* 15 from communication.c */
  {&n_total_particles,  TYPE_INT, 1, "n_part",            6 },         /* 16 from particle.c */
  {&n_particle_types,   TYPE_INT, 1, "n_part_types",      8 },         /* 17 from interaction_data.c */
  {&n_rigidbonds,       TYPE_INT, 1, "n_rigidbonds",      5 },         /* 18 from rattle.c */
  {node_grid,           TYPE_INT, 3, "node_grid",         2 },         /* 19 from grid.c */
  {&nptiso_gamma0,   TYPE_DOUBLE, 1, "nptiso_gamma0",    13 },         /* 20 from thermostat.c */
  {&nptiso_gammav,   TYPE_DOUBLE, 1, "nptiso_gammav",    13 },         /* 21 from thermostat.c */
  {&nptiso.p_ext,    TYPE_DOUBLE, 1, "npt_p_ext",         7 },         /* 22 from pressure.c */
  {&nptiso.p_inst,   TYPE_DOUBLE, 1, "npt_p_inst",       10 },         /* 23 from pressure.c */
  {&nptiso.p_inst_av,TYPE_DOUBLE, 1, "npt_p_inst_av",    10 },         /* 24 from pressure.c */
  {&nptiso.p_diff,   TYPE_DOUBLE, 1, "npt_p_diff",        7 },         /* 25 from pressure.c */
  {&nptiso.piston,   TYPE_DOUBLE, 1, "npt_piston",        6 },         /* 26 from pressure.c */
  {&periodic,          TYPE_BOOL, 3, "periodicity",       1 },         /* 27 from grid.c */
  {&skin,            TYPE_DOUBLE, 1, "skin",              2 },         /* 28 from integrate.c */
  {&temperature,     TYPE_DOUBLE, 1, "temperature",       2 },         /* 29 from thermostat.c */
  {&thermo_switch,      TYPE_INT, 1, "thermo_switch",     2 },         /* 30 from thermostat.c */
  {&sim_time,        TYPE_DOUBLE, 1, "time",              4 },         /* 31 from integrate.c */
  {&time_step,       TYPE_DOUBLE, 1, "time_step",         5 },         /* 32 from integrate.c */
  {&timing_samples,     TYPE_INT, 1, "timings",           4 },         /* 33 from tuning.c */
  {&transfer_rate,      TYPE_INT, 1, "transfer_rate",     2 },         /* 34 from imd.c */
  {&max_cut_nonbonded,TYPE_DOUBLE, 1, "max_cut_nonbonded",9 },         /* 35 from interaction_data.c */
  {&verlet_reuse,    TYPE_DOUBLE, 1, "verlet_reuse",      8 },         /* 36 from integrate.c */
  {&lattice_switch,     TYPE_INT, 1, "lattice_switch",    2 },         /* 37 from lattice.c */
  {&dpd_tgamma,      TYPE_DOUBLE, 1, "dpd_tgamma",        6 },         /* 38 from thermostat.c */
  {&dpd_tr_cut,      TYPE_DOUBLE, 1, "dpd_tr_cut",        6 },         /* 39 from thermostat.c */
  {&dpd_twf,            TYPE_INT, 1, "dpd_twf",           6 },         /* 40 from thermostat.c */
  {&dpd_wf,             TYPE_INT, 1, "dpd_wf",            5 },         /* 41 from thermostat.c */
  {adress_vars,      TYPE_DOUBLE, 7, "adress_vars",       1 },         /* 42  from adresso.c */
  {&max_cut_bonded,  TYPE_DOUBLE, 1, "max_cut_bonded",    9 },         /* 43 from interaction_data.c */
  {&transfer_rate,      TYPE_INT, 1, "vmd_transfer_rate", 5 },         /* 44 from imd_tcl.c */
  {&min_global_cut,  TYPE_DOUBLE, 1, "min_global_cut",    5 },         /* 45 from interaction_data.c */
  {&ghmc_nmd,           TYPE_INT, 1, "ghmc_nmd",          6 },         /* 46 from thermostat.c */
  {&ghmc_phi,        TYPE_DOUBLE, 1, "ghmc_phi",          6 },         /* 47 from thermostat.c */
  {&ghmc_mc_res,     TYPE_INT,    1, "ghmc_mc_res",       7 },         /* 48 from ghmc.c */
  {&ghmc_mflip,      TYPE_INT,    1, "ghmc_mflip",        7 },         /* 49 from ghmc.c */
  {&ghmc_tscale,     TYPE_INT,    1, "ghmc_tscale",       6 },         /* 50 from ghmc.c */
  {&lb_components,   TYPE_INT,    1, "lb_components",     2 },         /* 51 from ghmc.c */
  { NULL, 0, 0, NULL, 0 }
};
