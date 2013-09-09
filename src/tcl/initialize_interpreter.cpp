/*
  Copyright (C) 2012,2013 The ESPResSo project
  
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
#include <tcl.h>
#include <unistd.h>
#include <cstdlib>
#include <cstring>
#include <cstdio>

#include "global.hpp"
#include "adresso_tcl.hpp"
#include "binary_file_tcl.hpp"
#include "constraint_tcl.hpp"
#include "domain_decomposition_tcl.hpp"
#include "dpd_tcl.hpp"
#include "galilei_tcl.hpp"
#include "global_tcl.hpp"
#include "grid_tcl.hpp"
#include "iccp3m_tcl.hpp"
#include "integrate_tcl.hpp"
#include "interaction_data_tcl.hpp"
#include "lb_tcl.hpp"
#include "lj_tcl.hpp"
#include "maggs_tcl.hpp"
#include "metadynamics_tcl.hpp"
#include "p3m-dipolar_tcl.hpp"
#include "p3m_tcl.hpp"
#include "polymer_tcl.hpp"
#include "pressure_tcl.hpp"
#include "random_tcl.hpp"
#include "reaction_tcl.hpp"
#include "statistics_chain_tcl.hpp"
#include "statistics_cluster_tcl.hpp"
#include "statistics_correlation_tcl.hpp"
#include "statistics_fluid_tcl.hpp"
#include "statistics_observable_tcl.hpp"
#include "statistics_tcl.hpp"
#include "thermostat_tcl.hpp"
#include "virtual_sites_com_tcl.hpp"
#include "ghmc_tcl.hpp"

#ifdef TK
#include <tk.hpp>
#endif

/****************************************
 * various forwards
 *****************************************/

/** Implementation of the tcl command bin, which can be used
    to bin data into arbitrary bins. See \ref bin_tcl.c
*/
int tclcommand_bin(ClientData data, Tcl_Interp *interp,
		   int argc, char **argv);
/** Implementation of the Tcl command blockfile. Allows to read and write
    blockfile comfortably from Tcl. See \ref blockfile_tcl.c */
int tclcommand_blockfile(ClientData data, Tcl_Interp *interp,
	      int argc, char **argv);
/** implementation of the Tcl command cellsystem. See \ref cells_tcl.c */
int tclcommand_cellsystem(ClientData data, Tcl_Interp *interp,
	       int argc, char **argv);
/** replaces one of TCLs standart channels with a named pipe. See \ref channels_tcl.c */
int tclcommand_replacestdchannel(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);
/** Implements the Tcl command code_info.  It provides information on the
    Version, Compilation status and the debug status of the used
    code. See \ref config_tcl.c */
int tclcommand_code_info(ClientData data, Tcl_Interp *interp,
	 int argc, char **argv);
/** Set the CUDA device to use or retrieve information
    available devices. See \ref cuda_init_tcl.c */
int tclcommand_cuda(ClientData data, Tcl_Interp *interp,
		    int argc, char **argv);
/** VMD connection. See \ref imd_tcl.c */
int tclcommand_imd(ClientData data, Tcl_Interp *interp,
		   int argc, char **argv);
/** tcl procedure for nemd steering.
    USAGE: nemd \<n_slabs\> \<n_exchange\>   
    see also \ref tclcommand_nemd. See \ref nemd_tcl.c */
int tclcommand_nemd(ClientData data, Tcl_Interp *interp,
		    int argc, char **argv);
/** Collision detection. See \ref collision_tcl.c */
int tclcommand_on_collision(ClientData data, Tcl_Interp *interp, int argc, char **argv);
/** Implementation of the tcl command "part". This command allows to
    modify particle data. See \ref particle_data_tcl.c */
int tclcommand_part(ClientData data, Tcl_Interp *interp,
		    int argc, char **argv);
/** The C implementation of the tcl function uwerr. See \ref uwerr_tcl.c */
int tclcommand_uwerr(ClientData data, Tcl_Interp *interp, int argc, char *argv[]);
/** callback for \ref timing_samples. See \ref tuning_tcl.c */
int tclcallback_timings(Tcl_Interp *interp, void *data);

/// from \ref scriptsdir.c
char *get_default_scriptsdir();

/****************************************
 * Registration functions
 *****************************************/

#define REGISTER_COMMAND(name, routine)					\
  Tcl_CreateCommand(interp, name, (Tcl_CmdProc *)routine, 0, NULL);

static void register_tcl_commands(Tcl_Interp* interp) {
  /* in cells.c */
  REGISTER_COMMAND("cellsystem", tclcommand_cellsystem);
  /* in integrate.c */
  REGISTER_COMMAND("invalidate_system", tclcommand_invalidate_system);
  REGISTER_COMMAND("integrate", tclcommand_integrate);
  /* in global.c */
  REGISTER_COMMAND("setmd", tclcommand_setmd);
  /* in grid.c */
  REGISTER_COMMAND("change_volume", tclcommand_change_volume);
  /* in config_tcl.c */
  REGISTER_COMMAND("code_info", tclcommand_code_info);
  /* in interaction_data.c */
  REGISTER_COMMAND("inter",tclcommand_inter);
  /* in particle_data.c */
  REGISTER_COMMAND("part",tclcommand_part);
  /* in file binaryfile.c */
  REGISTER_COMMAND("writemd", tclcommand_writemd);
  REGISTER_COMMAND("readmd", tclcommand_readmd);
  /* in file statistics.c */
  REGISTER_COMMAND("analyze", tclcommand_analyze);
  /* in file polymer.c */
  REGISTER_COMMAND("polymer", tclcommand_polymer);
  REGISTER_COMMAND("counterions", tclcommand_counterions);
  REGISTER_COMMAND("salt", tclcommand_salt);
  REGISTER_COMMAND("velocities", tclcommand_velocities);
  REGISTER_COMMAND("maxwell_velocities", tclcommand_maxwell_velocities);
  REGISTER_COMMAND("crosslink", tclcommand_crosslink);
  REGISTER_COMMAND("diamond", tclcommand_diamond);
  REGISTER_COMMAND("icosaeder", tclcommand_icosaeder);
  /* in file imd.c */
  REGISTER_COMMAND("imd", tclcommand_imd);
  /* in file random.c */
  REGISTER_COMMAND("t_random", tclcommand_t_random);
  REGISTER_COMMAND("bit_random", tclcommand_bit_random);
  /* in file blockfile_tcl.c */
  REGISTER_COMMAND("blockfile", tclcommand_blockfile);
  /* in constraint.c */
  REGISTER_COMMAND("constraint", tclcommand_constraint);
  /* in uwerr.c */
  REGISTER_COMMAND("uwerr", tclcommand_uwerr);
  /* in nemd.c */
  REGISTER_COMMAND("nemd", tclcommand_nemd);
  /* in thermostat.c */
  REGISTER_COMMAND("thermostat", tclcommand_thermostat);
  /* in bin.c */
  REGISTER_COMMAND("bin", tclcommand_bin);
  /* in ghmc.c */
  REGISTER_COMMAND("ghmc", tclcommand_ghmc);
  REGISTER_COMMAND("save_state", tclcommand_save_state);
  REGISTER_COMMAND("load_state", tclcommand_load_state);
  /* in lb.c */

  REGISTER_COMMAND("lbfluid", tclcommand_lbfluid);
  REGISTER_COMMAND("lbnode", tclcommand_lbnode);
  REGISTER_COMMAND("lbboundary", tclcommand_lbboundary);
  /* here */
  REGISTER_COMMAND("replacestdchannel", tclcommand_replacestdchannel);
  /* in iccp3m.h */
  REGISTER_COMMAND("observable", tclcommand_observable);
  /* in statistics_obsrvable.h */
  REGISTER_COMMAND("correlation", tclcommand_correlation);
  /* in statistics_correlation.h */
#ifdef ELECTROSTATICS
#ifdef P3M
  REGISTER_COMMAND("iccp3m", tclcommand_iccp3m);
#endif 
#endif 
  /* in adresso.h */
  REGISTER_COMMAND("adress", tclcommand_adress);
#ifdef ADRESS
  /* #ifdef THERMODYNAMIC_FORCE */
  REGISTER_COMMAND("thermodynamic_force", tclcommand_thermodynamic_force);
  /* #endif */
  REGISTER_COMMAND("update_adress_weights", tclcommand_update_adress_weights);
#endif
#ifdef METADYNAMICS
  /* in metadynamics.c */
  REGISTER_COMMAND("metadynamics", tclcommand_metadynamics);
#endif
#ifdef LB_GPU
  /* in lbgpu_cfile.c */
  REGISTER_COMMAND("lbnode_extforce", tclcommand_lbnode_extforce_gpu);
#endif
#ifdef CUDA
  REGISTER_COMMAND("cuda", tclcommand_cuda);
#endif
  /* from collision.c */
#ifdef COLLISION_DETECTION
  REGISTER_COMMAND("on_collision", tclcommand_on_collision);
#endif
#ifdef CATALYTIC_REACTIONS
  REGISTER_COMMAND("reaction", tclcommand_reaction);
#endif
  REGISTER_COMMAND("kill_particle_motion", tclcommand_kill_particle_motion);
  REGISTER_COMMAND("kill_particle_forces", tclcommand_kill_particle_forces);
  REGISTER_COMMAND("system_CMS", tclcommand_system_CMS);
  REGISTER_COMMAND("system_CMS_velocity", tclcommand_system_CMS_velocity);
  REGISTER_COMMAND("galilei_transform", tclcommand_galilei_transform);
}

static void register_global_variables(Tcl_Interp *interp)
{
  /* register all writeable TCL variables with their callback functions */
  register_global_callback(FIELD_BOXL, tclcallback_box_l);
  register_global_callback(FIELD_MAXNUMCELLS, tclcallback_max_num_cells);
  register_global_callback(FIELD_MINNUMCELLS, tclcallback_min_num_cells);
  register_global_callback(FIELD_NODEGRID, tclcallback_node_grid);
  register_global_callback(FIELD_NPTISO_PDIFF, tclcallback_npt_p_diff);
  register_global_callback(FIELD_NPTISO_PISTON, tclcallback_npt_piston);
  register_global_callback(FIELD_PERIODIC, tclcallback_periodicity);
  register_global_callback(FIELD_SKIN, tclcallback_skin);
  register_global_callback(FIELD_SIMTIME, tclcallback_time);
  register_global_callback(FIELD_TIMESTEP, tclcallback_time_step);
  register_global_callback(FIELD_TIMINGSAMP, tclcallback_timings);
  register_global_callback(FIELD_MIN_GLOBAL_CUT, tclcallback_min_global_cut);
}

int appinit(Tcl_Interp *interp)
{
  if (Tcl_Init(interp) == TCL_ERROR)
    return TCL_ERROR;

#ifdef TK
  if (Tk_Init(interp) == TCL_ERROR)
    return TCL_ERROR;
#endif

  /*
    installation of tcl commands
  */
  register_tcl_commands(interp);
  register_global_variables(interp);

  /* evaluate the Tcl initialization script */
  char *scriptdir = getenv("ESPRESSO_SCRIPTS");
  if (!scriptdir)
    scriptdir = get_default_scriptsdir();
  
  /*  fprintf(stderr,"Script directory: %s\n", scriptdir);*/

  char cwd[1024];
  if ((getcwd(cwd, 1024) == NULL) || (chdir(scriptdir) != 0)) {
    fprintf(stderr,
	    "\n\ncould not change to script dir %s, please check ESPRESSO_SCRIPTS.\n\n\n",
	    scriptdir);
    exit(1);
  }
  if (Tcl_EvalFile(interp, "init.tcl") == TCL_ERROR) {
    fprintf(stderr, "\n\nerror in initialization script: %s\n\n\n",
	    Tcl_GetStringResult(interp));
    exit(1);
  }
  if (chdir(cwd) != 0) {
    fprintf(stderr,
	    "\n\ncould not change back to execution dir %s ????\n\n\n",
	    cwd);
    exit(1);
  }

  return (TCL_OK);
}
