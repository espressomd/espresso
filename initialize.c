// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2009; all rights reserved unless otherwise stated.
/** \file initialize.c
    Implementation of \ref initialize.h "initialize.h"
*/
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "utils.h"
#include "initialize.h"
#include "global.h"
#include "particle_data.h"
#include "interaction_data.h"
#include "binary_file.h"
#include "integrate.h"
#include "statistics.h"
#include "energy.h"
#include "pressure.h"
#include "polymer.h"
#include "imd.h"
#include "random.h"
#include "communication.h"
#include "blockfile_tcl.h"
#include "cells.h"
#include "grid.h"
#include "thermostat.h"
#include "rotation.h"
#include "p3m.h"
#include "ewald.h"
#include "fft.h"
#include "mmm1d.h"
#include "mmm2d.h"
#include "maggs.h"
#include "elc.h"
#include "lb.h"
#include "ghosts.h"
#include "debye_hueckel.h"
#include "reaction_field.h"
#include "forces.h"
#include "uwerr.h"
#include "nsquare.h"
#include "nemd.h"
#include "domain_decomposition.h"
#include "errorhandling.h"
#include "rattle.h"
#include "bin.h"
#include "lattice.h"
#include "lb-boundaries.h"
#include "iccp3m.h" /* -iccp3m- */
#include "adresso.h"
#include "metadynamics.h"

/** whether before integration the thermostat has to be reinitialized */
static int reinit_thermo = 1;
static int reinit_electrostatics = 0;
static int reinit_magnetostatics = 0;

static void init_tcl(Tcl_Interp *interp);

int on_program_start(Tcl_Interp *interp)
{
  EVENT_TRACE(fprintf(stderr, "%d: on_program_start\n", this_node));
  /*
    call the initialization of the modules here
  */
  init_random();
  init_bit_random();

  setup_node_grid();
  /* calculate initial minimimal number of cells (see min_num_cells_callback) */
  min_num_cells = calc_processor_min_num_cells();

  cells_pre_init();
  ghost_init();
  /* Initialise force and energy tables */
  force_and_energy_tables_init();
#ifdef ADRESS
#ifdef INTERFACE_CORRECTION
  adress_force_and_energy_tables_init();
#endif
  /** #ifdef THERMODYNAMIC_FORCE */
  tf_tables_init();
  /** #endif */
#endif
#ifdef ELP3M
  fft_pre_init();
#endif

#ifdef LB
  lb_pre_init();
#endif


  /*
    call all initializations to do only on the master node here.
  */
  if (this_node == 0) {
    /* interaction_data.c: make sure 0<->0 ia always exists */
    make_particle_type_exist(0);
    
    init_tcl(interp);
  }
  return TCL_OK;
}

void on_integration_start()
{
  char *errtext;
  int i;

  EVENT_TRACE(fprintf(stderr, "%d: on_integration_start\n", this_node));
  INTEG_TRACE(fprintf(stderr,"%d: on_integration_start: reinit_thermo = %d, resort_particles=%d\n",
		      this_node,reinit_thermo,resort_particles));

  /********************************************/
  /* sanity checks                            */
  /********************************************/

  if ( time_step < 0.0 ) {
    errtext = runtime_error(128);
    ERROR_SPRINTF(errtext, "{010 time_step not set} ");
  }
  if ( skin < 0.0 ) {
    errtext = runtime_error(128);
    ERROR_SPRINTF(errtext,"{011 skin not set} ");
  }
  if ( temperature < 0.0 ) {
    errtext = runtime_error(128);
    ERROR_SPRINTF(errtext,"{012 thermostat not initialized} ");
  }

  for (i = 0; i < 3; i++)
    if (local_box_l[i] < max_range) {
      errtext = runtime_error(128 + TCL_INTEGER_SPACE);
      ERROR_SPRINTF(errtext,"{013 box_l in direction %d is still too small} ", i);
    }
  
#ifdef NPT
  if (integ_switch == INTEG_METHOD_NPT_ISO) {
    if (nptiso.piston <= 0.0) {
      char *errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"{014 npt on, but piston mass not set} ");
    }
    if(cell_structure.type!=CELL_STRUCTURE_DOMDEC) {
      char *errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"{014 npt requires domain decomposition cellsystem} ");
    }
  
#ifdef ELECTROSTATICS

    switch(coulomb.method) {
      case COULOMB_NONE:  break;
#ifdef ELP3M
      case COULOMB_P3M:   break;
#endif /*ELP3M*/
      case COULOMB_EWALD: break;
      default: {
        char *errtext = runtime_error(128);
        ERROR_SPRINTF(errtext,"{014 npt only works with Ewald sum or P3M} ");
      }
    }
#endif /*ELECTROSTATICS*/
  }
  
#endif /*NPT*/

  if (!check_obs_calc_initialized()) return;

#ifdef LB
  if(lattice_switch & LATTICE_LB) {
    if (lbpar.agrid < 0.0) {
      errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"{098 Lattice Boltzmann agrid not set} ");
    }
    if (lbpar.tau < 0.0) {
      errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"{099 Lattice Boltzmann time step not set} ");
    }
    if (lbpar.rho < 0.0) {
      errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"{100 Lattice Boltzmann fluid density not set} ");
    }
    if (lbpar.viscosity < 0.0) {
      errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"{101 Lattice Boltzmann fluid viscosity not set} ");
    }
  }
#endif

#ifdef METADYNAMICS
  meta_init();
#endif

  /********************************************/
  /* end sanity checks                        */
  /********************************************/


  /* Prepare the thermostat */
  if (reinit_thermo) {
    thermo_init();
    reinit_thermo = 0;
    recalc_forces = 1;
  }

  /* Ensemble preparation: NVT or NPT */
  integrate_ensemble_init();

  /* Update particle and observable information for routines in statistics.c */
  invalidate_obs();
  freePartCfg();

  on_observable_calc();
}

void on_observable_calc()
{
  /* Prepare particle structure: Communication step: number of ghosts and ghost information */

  if(resort_particles) {
    cells_resort_particles(CELL_GLOBAL_EXCHANGE);
    resort_particles = 0;
  }

#ifdef ELECTROSTATICS  
  if(reinit_electrostatics) {
    EVENT_TRACE(fprintf(stderr, "%d: reinit_electrostatics\n", this_node));
    switch (coulomb.method) {
#ifdef ELP3M
    case COULOMB_ELC_P3M:
    case COULOMB_P3M:
      P3M_count_charged_particles();
      break;
#endif
    case COULOMB_EWALD:
      EWALD_count_charged_particles();
      break;
    case COULOMB_MAGGS: 
      Maggs_init(); 
      break;
    default: break;
    }
    reinit_electrostatics = 0;
  }
#endif /*ifdef ELECTROSTATICS */

#ifdef MAGNETOSTATICS
  if(reinit_magnetostatics) {
    EVENT_TRACE(fprintf(stderr, "%d: reinit_magnetostatics\n", this_node));
    switch (coulomb.Dmethod) {
    #ifdef ELP3M
    case DIPOLAR_MDLC_P3M:
    case DIPOLAR_P3M:
      P3M_count_magnetic_particles();
      break;
#endif
    default: break;
    }
    reinit_magnetostatics = 0;
  }
#endif /*ifdef ELECTROSTATICS */
}

void on_particle_change()
{
  // EVENT_TRACE(fprintf(stderr, "%d: on_particle_change\n", this_node));
  resort_particles = 1;
  reinit_electrostatics = 1;
  reinit_magnetostatics = 1;
  rebuild_verletlist = 1;

  invalidate_obs();

  /* the particle information is no longer valid */
  freePartCfg();
}

void on_coulomb_change()
{
  EVENT_TRACE(fprintf(stderr, "%d: on_coulomb_change\n", this_node));
  invalidate_obs();

#ifdef ELECTROSTATICS
  if(temperature > 0.0)
    coulomb.prefactor = coulomb.bjerrum * temperature; 
  else
    coulomb.prefactor = coulomb.bjerrum;
  switch (coulomb.method) {
#ifdef ELP3M
  case COULOMB_ELC_P3M:
    ELC_init();
    // fall through
  case COULOMB_P3M:
    P3M_init_charges();
    integrate_vv_recalc_maxrange();
    on_parameter_change(FIELD_MAXRANGE);
    break;
#endif
  case COULOMB_EWALD:
    EWALD_init();
    integrate_vv_recalc_maxrange();
    on_parameter_change(FIELD_MAXRANGE);
    break;
  case COULOMB_MMM1D:
    MMM1D_init();
    break;
  case COULOMB_MMM2D:
    MMM2D_init();
    break;
  case COULOMB_MAGGS: 
    Maggs_init(); 
    break;
  default: break;
  }

  recalc_forces = 1;
#endif  /* ifdef ELECTROSTATICS */

#ifdef MAGNETOSTATICS
  if(temperature > 0.0)
    coulomb.Dprefactor = coulomb.Dbjerrum * temperature; 
  else
    coulomb.Dprefactor = coulomb.Dbjerrum;
  
  switch (coulomb.Dmethod) {
#ifdef ELP3M
    case DIPOLAR_MDLC_P3M:
       // fall through
  case DIPOLAR_P3M:
    P3M_init_dipoles();
    integrate_vv_recalc_maxrange();
    on_parameter_change(FIELD_MAXRANGE);
    break;
#endif
  default: break;
  }

  recalc_forces = 1;
#endif  /* ifdef MAGNETOSTATICS */

}

void on_short_range_ia_change()
{
  EVENT_TRACE(fprintf(stderr, "%d: on_short_range_ia_changes\n", this_node));
  invalidate_obs();

  integrate_vv_recalc_maxrange();
  on_parameter_change(FIELD_MAXRANGE);

  recalc_forces = 1;
}

void on_constraint_change()
{
  EVENT_TRACE(fprintf(stderr, "%d: on_constraint_change\n", this_node));
  invalidate_obs();

#ifdef LB
#ifdef CONSTRAINTS
  if(lattice_switch & LATTICE_LB) {
    lb_init_boundaries();
  }
#endif
#endif

  recalc_forces = 1;
}

void on_lb_boundary_change()
{
  EVENT_TRACE(fprintf(stderr, "%d: on_lb_boundary_change\n", this_node));
  invalidate_obs();

#ifdef LB_BOUNDARIES
  printf("executing on_lb_boundary_change on node %d\n", this_node);
  
  if(lattice_switch & LATTICE_LB) {
    lb_init_boundaries();
  }
#endif

  recalc_forces = 1;
}

void on_cell_structure_change()
{
  EVENT_TRACE(fprintf(stderr, "%d: on_cell_structure_change\n", this_node));
  on_coulomb_change();
}

void on_resort_particles()
{
  EVENT_TRACE(fprintf(stderr, "%d: on_resort_particles\n", this_node));
#ifdef ELECTROSTATICS
  switch (coulomb.method) {
#ifdef ELP3M
  case COULOMB_ELC_P3M:
    ELC_on_resort_particles();
    break;
#endif
  case COULOMB_MMM2D:
    MMM2D_on_resort_particles();
    break;
  case COULOMB_EWALD:
    EWALD_on_resort_particles();
    break;
  default: break;
  }
#endif /* ifdef ELECTROSTATICS */


#ifdef MAGNETOSTATICS
  switch (coulomb.Dmethod) {
#ifdef ELP3M
  case DIPOLAR_MDLC_P3M:
    /* dlc_on_resort_particles();   NOT NECESSARY DUE TO HOW WE COMPUTE THINGS*/
    break;
#endif
 default: break;
  }
#endif /* ifdef MAGNETOSTATICS*/

}

#ifdef NPT
void on_NpT_boxl_change(double scal1) {
  grid_changed_box_l();
  
#ifdef ELECTROSTATICS
  switch(coulomb.method) {
#ifdef ELP3M
  case COULOMB_P3M:
    P3M_scaleby_box_l_charges();
    integrate_vv_recalc_maxrange();
    break;
#endif
  case COULOMB_EWALD:
    EWALD_scaleby_box_l();
    integrate_vv_recalc_maxrange();
    break;
  default: break;
  }
#endif

#ifdef MAGNETOSTATICS
  switch(coulomb.Dmethod) {
#ifdef ELP3M
  case DIPOLAR_P3M:
    P3M_scaleby_box_l_dipoles();
    integrate_vv_recalc_maxrange();
    break;
#endif
  default: break;
  }
#endif

  if(cell_structure.type==CELL_STRUCTURE_DOMDEC)
    dd_NpT_update_cell_grid(scal1);
}
#endif

void on_parameter_change(int field)
{
  /* to prevent two on_coulomb_change */
#if defined(ELECTROSTATICS) || defined(MAGNETOSTATICS)
  int cc = 0;
#endif

  EVENT_TRACE(fprintf(stderr, "%d: on_parameter_change %s\n", this_node, fields[field].name));

  if (field == FIELD_SKIN) {
    integrate_vv_recalc_maxrange();
    on_parameter_change(FIELD_MAXRANGE);
  }

  if (field == FIELD_NODEGRID)
    grid_changed_n_nodes();
  if (field == FIELD_BOXL || field == FIELD_NODEGRID)
    grid_changed_box_l();
  if (field == FIELD_TIMESTEP || field == FIELD_TEMPERATURE || field == FIELD_LANGEVIN_GAMMA || field == FIELD_DPD_TGAMMA
      || field == FIELD_DPD_GAMMA || field == FIELD_NPTISO_G0 || field == FIELD_NPTISO_GV || field == FIELD_NPTISO_PISTON )
    reinit_thermo = 1;

#ifdef NPT
  if ((field == FIELD_INTEG_SWITCH) && (integ_switch != INTEG_METHOD_NPT_ISO))
    nptiso.invalidate_p_vel = 1;  
#endif

#ifdef ADRESS
//   if (field == FIELD_BOXL)
//    adress_changed_box_l();
#endif

#ifdef ELECTROSTATICS
  switch (coulomb.method) {
#ifdef ELP3M
  case COULOMB_ELC_P3M:
    if (field == FIELD_TEMPERATURE || field == FIELD_BOXL)
      cc = 1;
    // fall through
  case COULOMB_P3M:
    if (field == FIELD_TEMPERATURE || field == FIELD_NODEGRID || field == FIELD_SKIN)
      cc = 1;
    else if (field == FIELD_BOXL) {
      P3M_scaleby_box_l_charges();
      integrate_vv_recalc_maxrange(); 
    }
    break;
#endif
  case COULOMB_EWALD:
    if (field == FIELD_TEMPERATURE || field == FIELD_SKIN)
      cc = 1;
    else if (field == FIELD_BOXL) {
      EWALD_scaleby_box_l();
      integrate_vv_recalc_maxrange(); 
    }
    break;
  case COULOMB_DH:
    if (field == FIELD_TEMPERATURE)
      cc = 1;
    break;
  case COULOMB_RF:
  case COULOMB_INTER_RF:
    if (field == FIELD_TEMPERATURE)
      cc = 1;
    break;
  case COULOMB_MMM1D:
    if (field == FIELD_TEMPERATURE || field == FIELD_BOXL)
      cc = 1;
    break;
  case COULOMB_MMM2D:
    if (field == FIELD_TEMPERATURE || field == FIELD_BOXL || field == FIELD_NLAYERS)
      cc = 1;
    break;
  case COULOMB_MAGGS:
    /* Maggs electrostatics needs ghost velocities */
    on_ghost_flags_change();
    cells_re_init(CELL_STRUCTURE_CURRENT);    
    break;
  default: break;
  }
#endif /*ifdef ELECTROSTATICS */

#ifdef MAGNETOSTATICS
  switch (coulomb.Dmethod) {
   #ifdef ELP3M
    case DIPOLAR_MDLC_P3M:
     if (field == FIELD_TEMPERATURE || field == FIELD_BOXL)
       cc = 1;
      // fall through
    case DIPOLAR_P3M:
      if (field == FIELD_TEMPERATURE || field == FIELD_NODEGRID || field == FIELD_SKIN)
        cc = 1;
      else if (field == FIELD_BOXL) {
        P3M_scaleby_box_l_dipoles();
        integrate_vv_recalc_maxrange(); 
      }
      break;
#endif
  default: break;
  }
#endif /*ifdef MAGNETOSTATICS */

#if defined(ELECTROSTATICS) || defined(MAGNETOSTATICS)
  if (cc)
    on_coulomb_change();
#endif

  /* DPD needs ghost velocities, other thermostats not */
  if (field == FIELD_THERMO_SWITCH) {
    on_ghost_flags_change();
    cells_re_init(CELL_STRUCTURE_CURRENT);
  }

  if (field == FIELD_MAXRANGE)
    rebuild_verletlist = 1;

  switch (cell_structure.type) {
  case CELL_STRUCTURE_LAYERED:
    if (field == FIELD_NODEGRID) {
      if (node_grid[0] != 1 || node_grid[1] != 1) {
	char *errtext = runtime_error(128);
	ERROR_SPRINTF(errtext, "{091 layered cellsystem requires 1 1 n node grid} ");
      }
    }
    if (field == FIELD_BOXL || field == FIELD_MAXRANGE || field == FIELD_THERMO_SWITCH)
      cells_re_init(CELL_STRUCTURE_LAYERED);
    break;
  case CELL_STRUCTURE_DOMDEC:
    if (field == FIELD_BOXL || field == FIELD_NODEGRID || field == FIELD_MAXRANGE ||
	field == FIELD_MINNUMCELLS || field == FIELD_MAXNUMCELLS || field == FIELD_THERMO_SWITCH)
      cells_re_init(CELL_STRUCTURE_DOMDEC);
    break;
  }

#ifdef LB
  /* LB needs ghost velocities */
  if (field == FIELD_LATTICE_SWITCH) {
    on_ghost_flags_change();
    cells_re_init(CELL_STRUCTURE_CURRENT);
  }

  if (lattice_switch & LATTICE_LB) {
    if (field == FIELD_TEMPERATURE) {
      lb_reinit_parameters();
    }

    if (field == FIELD_BOXL || field == FIELD_CELLGRID || field == FIELD_NNODES || field == FIELD_NODEGRID) {
      lb_init();
    }
  }
#endif
}

#ifdef LB
void on_lb_params_change(int field) {

  if (field == LBPAR_AGRID) {
    lb_init();
  }
  if (field == LBPAR_DENSITY) {
    lb_reinit_fluid();
  }

  lb_reinit_parameters();

}
#endif

void on_ghost_flags_change()
{
  /* that's all we change here */
  extern int ghosts_have_v;
  ghosts_have_v = 0;
  
  /* DPD and LB need also ghost velocities */
  if (thermo_switch & THERMO_DPD)
    ghosts_have_v = 1;
#ifdef LB
  if (lattice_switch & LATTICE_LB)
    ghosts_have_v = 1;
#endif
#ifdef BOND_CONSTRAINT
  else if (n_rigidbonds)
    ghosts_have_v = 1;
#endif
#ifdef ELECTROSTATICS
  /* Maggs electrostatics needs ghost velocities too */
  else if(coulomb.method == COULOMB_MAGGS)
      ghosts_have_v = 1;
#endif
#ifdef INTER_DPD
  //maybe we have to add a new global to differ between compile in and acctual use.
  ghosts_have_v = 1;
#endif
#ifdef VIRTUAL_SITES 
  //VIRUTAL_SITES need v to update v of virtual sites
  ghosts_have_v = 1;
#endif
}

static void init_tcl(Tcl_Interp *interp)
{
  char cwd[1024];
  char *scriptdir;

  /*
    installation of tcl commands
  */

#define REGISTER_COMMAND(name, routine)					\
  Tcl_CreateCommand(interp, name, (Tcl_CmdProc *)routine, 0, NULL);

  /* in cells.c */
  REGISTER_COMMAND("cellsystem", cellsystem);
  /* in integrate.c */
  REGISTER_COMMAND("invalidate_system", invalidate_system);
  REGISTER_COMMAND("integrate", integrate);
  /* in global.c */
  REGISTER_COMMAND("setmd", setmd);
  /* in grid.c */
  REGISTER_COMMAND("change_volume", change_volume);
  /* in interaction_data.c */
  REGISTER_COMMAND("code_info", code_info);
  /* in interaction_data.c */
  REGISTER_COMMAND("inter", inter);
  /* in particle_data.c */
  REGISTER_COMMAND("part", part);
  /* in file binaryfile.c */
  REGISTER_COMMAND("writemd", writemd);
  REGISTER_COMMAND("readmd", readmd);
  /* in file statistics.c */
  REGISTER_COMMAND("analyze", analyze);
  /* in file polymer.c */
  REGISTER_COMMAND("polymer", polymer);
  REGISTER_COMMAND("counterions", counterions);
  REGISTER_COMMAND("salt", salt);
  REGISTER_COMMAND("velocities", velocities);
  REGISTER_COMMAND("maxwell_velocities", maxwell_velocities);
  REGISTER_COMMAND("crosslink", crosslink);
  REGISTER_COMMAND("diamond", diamond);
  REGISTER_COMMAND("icosaeder", icosaeder);
  /* in file imd.c */
  REGISTER_COMMAND("imd", imd);
  /* in file random.c */
  REGISTER_COMMAND("t_random", t_random);
  REGISTER_COMMAND("bit_random", bit_random);
  /* in file blockfile_tcl.c */
  REGISTER_COMMAND("blockfile", blockfile);
  /* in constraint.c */
  REGISTER_COMMAND("constraint", constraint);
  /* in lb-boundaries.c */
  REGISTER_COMMAND("lb_boundary", lb_boundary);
  /* in uwerr.c */
  REGISTER_COMMAND("uwerr", uwerr);
  /* in nemd.c */
  REGISTER_COMMAND("nemd", nemd);
  /* in thermostat.c */
  REGISTER_COMMAND("thermostat", thermostat);
  /* in bin.c */
  REGISTER_COMMAND("bin", bin);
  /* in lb.c */
  REGISTER_COMMAND("lbfluid", lbfluid_cmd);
  /* in utils.h */
  REGISTER_COMMAND("replacestdchannel", replacestdchannel);
  /* in iccp3m.h */
#ifdef ELECTROSTATICS
#ifdef ELP3M
  REGISTER_COMMAND("iccp3m", iccp3m);
#endif 
#endif 
  /* in adresso.h */
  REGISTER_COMMAND("adress", adress_tcl);
#ifdef ADRESS
  /** #ifdef THERMODYNAMIC_FORCE */
  REGISTER_COMMAND("thermodynamic_force", tf_tcl);
  /** #endif */
  REGISTER_COMMAND("update_adress_weights", manual_update_weights);
#endif
#ifdef METADYNAMICS
  /* in metadynamics.c */
  REGISTER_COMMAND("metadynamics", metadynamics);
#endif
  /* evaluate the Tcl initialization script */
  scriptdir = getenv("ESPRESSO_SCRIPTS");
  fprintf(stderr,"%d: Script directory: %s\n",this_node,scriptdir);
  if (!scriptdir)
    scriptdir= "scripts";


  if ((getcwd(cwd, 1024) == NULL) || (chdir(scriptdir) != 0)) {
    fprintf(stderr,
	    "\n\ncould not change to script dir %s, please check ESPRESSO_SCRIPTS.\n\n\n",
	    scriptdir);
    exit(-1);
  }
  if (Tcl_EvalFile(interp, "init.tcl") == TCL_ERROR) {
    fprintf(stderr, "\n\nerror in initialization script: %s\n\n\n",
	    Tcl_GetStringResult(interp));
    exit(-1);
  }
  if (chdir(cwd) != 0) {
    fprintf(stderr,
	    "\n\ncould not change back to execution dir %s ????\n\n\n",
	    cwd);
    exit(-1);
  }
}
