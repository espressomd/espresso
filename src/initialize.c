/*
  Copyright (C) 2010,2011,2012 The ESPResSo project
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
/** \file initialize.c
    Implementation of \ref initialize.h "initialize.h"
*/
#include "utils.h"
#include "initialize.h"
#include "global.h"
#include "particle_data.h"
#include "interaction_data.h"
#include "reaction.h"
#include "statistics.h"
#include "energy.h"
#include "pressure.h"
#include "polymer.h"
#include "imd.h"
#include "random.h"
#include "communication.h"
#include "cells.h"
#include "grid.h"
#include "thermostat.h"
#include "rotation.h"
#include "p3m.h"
#include "p3m-dipolar.h"
#include "mmm1d.h"
#include "mmm2d.h"
#include "maggs.h"
#include "elc.h"
#include "lb.h"
#include "ghosts.h"
#include "debye_hueckel.h"
#include "reaction_field.h"
#include "forces.h"
#include "nsquare.h"
#include "nemd.h"
#include "domain_decomposition.h"
#include "errorhandling.h"
#include "rattle.h"
#include "lattice.h"
#include "iccp3m.h" /* -iccp3m- */
#include "adresso.h"
#include "metadynamics.h"
#include "statistics_observable.h"
#include "statistics_correlation.h"
#include "lb-boundaries.h"
#include "ghmc.h"
#include "domain_decomposition.h"

/** whether the thermostat has to be reinitialized before integration */
static int reinit_thermo = 1;
static int reinit_electrostatics = 0;
static int reinit_magnetostatics = 0;
#ifdef LB_GPU
static int lb_reinit_particles_gpu = 1;
#endif

void on_program_start()
{
  EVENT_TRACE(fprintf(stderr, "%d: on_program_start\n", this_node));

  /* tell Electric fence that we do realloc(0) on purpose. */
#ifdef EFENCE
  extern int EF_ALLOW_MALLOC_0;
  EF_ALLOW_MALLOC_0 = 1;
#endif

  register_sigint_handler();

  if (this_node == 0) {
    /* master node */
#ifdef FORCE_CORE
    /* core should be the last exit handler (process dies) */
    atexit(core);
#endif
    atexit(mpi_stop);
  }

  /*
    call the initialization of the modules here
  */
  init_random();
  init_bit_random();

  setup_node_grid();
  /* calculate initial minimimal number of cells (see tclcallback_min_num_cells) */
  min_num_cells = calc_processor_min_num_cells();

  /* initially go for domain decomposition */
  dd_topology_init(&local_cells);

  ghost_init();
  /* Initialise force and energy tables */
  force_and_energy_tables_init();
#ifdef ADRESS
#ifdef INTERFACE_CORRECTION
  adress_force_and_energy_tables_init();
#endif
  /* #ifdef THERMODYNAMIC_FORCE */
  tf_tables_init();
  /* #endif */
#endif
#ifdef P3M
  p3m_pre_init();
#endif
#ifdef DP3M
  dp3m_pre_init();
#endif
#ifdef LB
  lb_pre_init();
#endif
#ifdef REACTIONS
  reaction.back_rate=-1.0;
#endif

  /*
    call all initializations to do only on the master node here.
  */
  if (this_node == 0) {
    /* interaction_data.c: make sure 0<->0 ia always exists */
    make_particle_type_exist(0);
  }
}


void on_integration_start()
{
  char *errtext;

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
  
#ifdef NPT
  if (integ_switch == INTEG_METHOD_NPT_ISO) {
    if (nptiso.piston <= 0.0) {
      char *errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"{014 npt on, but piston mass not set} ");
    }
  
#ifdef ELECTROSTATICS

    switch(coulomb.method) {
    case COULOMB_NONE:  break;
#ifdef P3M
    case COULOMB_P3M:   break;
#endif /*P3M*/
    default: {
      char *errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"{014 npt only works with P3M} ");
    }
    }
#endif /*ELECTROSTATICS*/

#ifdef DIPOLES

    switch (coulomb.Dmethod) {
    case DIPOLAR_NONE: break;
#ifdef DP3M
    case DIPOLAR_P3M: break;
#endif
    default: {
      char *errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"NpT does not work with your dipolar method, please use P3M.");
    }
    }
#endif  /* ifdef DIPOLES */
  }
#endif /*NPT*/

  if (!check_obs_calc_initialized()) return;

#ifdef LB
  if(lattice_switch & LATTICE_LB) {
    if (lbpar.agrid <= 0.0) {
      errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"{098 Lattice Boltzmann agrid not set} ");
    }
    if (lbpar.tau <= 0.0) {
      errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"{099 Lattice Boltzmann time step not set} ");
    }
    if (lbpar.rho <= 0.0) {
      errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"{100 Lattice Boltzmann fluid density not set} ");
    }
    if (lbpar.viscosity <= 0.0) {
      errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"{101 Lattice Boltzmann fluid viscosity not set} ");
    }
    if (dd.use_vList && skin>=lbpar.agrid/2.0) {
      errtext = runtime_error(128);
      ERROR_SPRINTF(errtext, "{104 LB requires either no Verlet lists or that the skin of the verlet list to be less than half of lattice-Boltzmann grid spacing.} ");
    }
  }
#endif
#ifdef LB_GPU
if(this_node == 0){
  if(lattice_switch & LATTICE_LB_GPU) {
    if (lbpar_gpu.agrid < 0.0) {
      errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"{098 Lattice Boltzmann agrid not set} ");
    }
    if (lbpar_gpu.tau < 0.0) {
      errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"{099 Lattice Boltzmann time step not set} ");
    }
    if (lbpar_gpu.rho < 0.0) {
      errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"{100 Lattice Boltzmann fluid density not set} ");
    }
    if (lbpar_gpu.viscosity < 0.0) {
      errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"{101 Lattice Boltzmann fluid viscosity not set} ");
    }
    if (lb_reinit_particles_gpu) {
	    lb_realloc_particles_gpu();
	    lb_reinit_particles_gpu = 0;
    }
  }
}

#endif

#ifdef METADYNAMICS
  meta_init();
#endif

#ifdef REACTIONS
if(reaction.rate != 0.0) {
  if(max_cut < reaction.range) {
    errtext = runtime_error(128);
    ERROR_SPRINTF(errtext,"{105 Reaction range of %f exceeds maximum cutoff of %f} ", reaction.range, max_cut);
  }
}
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
  EVENT_TRACE(fprintf(stderr, "%d: on_observable_calc\n", this_node));
  /* Prepare particle structure: Communication step: number of ghosts and ghost information */

  if(resort_particles)
    cells_resort_particles(CELL_GLOBAL_EXCHANGE);

#ifdef ELECTROSTATICS  
  if(reinit_electrostatics) {
    EVENT_TRACE(fprintf(stderr, "%d: reinit_electrostatics\n", this_node));
    switch (coulomb.method) {
#ifdef P3M
    case COULOMB_ELC_P3M:
    case COULOMB_P3M:
      p3m_count_charged_particles();
      break;
#endif
    case COULOMB_MAGGS: 
      maggs_init(); 
      break;
    default: break;
    }
    reinit_electrostatics = 0;
  }
#endif /*ifdef ELECTROSTATICS */

#ifdef DIPOLES
  if(reinit_magnetostatics) {
    EVENT_TRACE(fprintf(stderr, "%d: reinit_magnetostatics\n", this_node));
    switch (coulomb.Dmethod) {
#ifdef DP3M
    case DIPOLAR_MDLC_P3M:
      // fall through
    case DIPOLAR_P3M:
      dp3m_count_magnetic_particles();
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
  EVENT_TRACE(fprintf(stderr, "%d: on_particle_change\n", this_node));
  resort_particles = 1;
  reinit_electrostatics = 1;
  reinit_magnetostatics = 1;

#ifdef LB_GPU
  lb_reinit_particles_gpu = 1;
#endif

  invalidate_obs();

  /* the particle information is no longer valid */
  freePartCfg();
}

void on_coulomb_change()
{
  EVENT_TRACE(fprintf(stderr, "%d: on_coulomb_change\n", this_node));
  invalidate_obs();

  recalc_coulomb_prefactor();

#ifdef ELECTROSTATICS
  switch (coulomb.method) {
  case COULOMB_DH:
    break;    
#ifdef P3M
  case COULOMB_ELC_P3M:
    ELC_init();
    // fall through
  case COULOMB_P3M:
    p3m_init();
    break;
#endif
  case COULOMB_MMM1D:
    MMM1D_init();
    break;
  case COULOMB_MMM2D:
    MMM2D_init();
    break;
  case COULOMB_MAGGS: 
    maggs_init();
    /* Maggs electrostatics needs ghost velocities */
    on_ghost_flags_change();
    break;
  default: break;
  }
#endif  /* ifdef ELECTROSTATICS */

#ifdef DIPOLES
  switch (coulomb.Dmethod) {
#ifdef DP3M
  case DIPOLAR_MDLC_P3M:
    // fall through
  case DIPOLAR_P3M:
    dp3m_init();
    break;
#endif
  default: break;
  }
#endif  /* ifdef DIPOLES */

  /* all Coulomb methods have a short range part, aka near field
     correction. Even in case of switching off, we should call this,
     since the required cutoff might have reduced. */
  on_short_range_ia_change();

  recalc_forces = 1;
}

void on_short_range_ia_change()
{
  EVENT_TRACE(fprintf(stderr, "%d: on_short_range_ia_changes\n", this_node));
  invalidate_obs();

  recalc_maximal_cutoff();
  cells_on_geometry_change(0);

  recalc_forces = 1;
}

void on_constraint_change()
{
  EVENT_TRACE(fprintf(stderr, "%d: on_constraint_change\n", this_node));
  invalidate_obs();
  recalc_forces = 1;
}

void on_lbboundary_change()
{
  EVENT_TRACE(fprintf(stderr, "%d: on_lbboundary_change\n", this_node));
  invalidate_obs();

#ifdef LB_BOUNDARIES
  if(lattice_switch & LATTICE_LB) {
    lb_init_boundaries();
  }
#endif
#ifdef LB_BOUNDARIES_GPU
  if(this_node == 0){
    if(lattice_switch & LATTICE_LB_GPU) {
      lb_init_boundaries();
    }
  }
#endif

  recalc_forces = 1;
}

void on_resort_particles()
{
  EVENT_TRACE(fprintf(stderr, "%d: on_resort_particles\n", this_node));
#ifdef ELECTROSTATICS
  switch (coulomb.method) {
#ifdef P3M
  case COULOMB_ELC_P3M:
    ELC_on_resort_particles();
    break;
#endif
  case COULOMB_MMM2D:
    MMM2D_on_resort_particles();
    break;
  default: break;
  }
#endif /* ifdef ELECTROSTATICS */
  
  /* DIPOLAR interactions so far don't need this */

  recalc_forces = 1;
}

void on_boxl_change() {
  EVENT_TRACE(fprintf(stderr, "%d: on_boxl_change\n", this_node));

  /* Now give methods a chance to react to the change in box length */
#ifdef ELECTROSTATICS
  switch(coulomb.method) {
#ifdef P3M
  case COULOMB_ELC_P3M:
    ELC_init();
    // fall through
  case COULOMB_P3M:
    p3m_scaleby_box_l();
    break;
#endif
  case COULOMB_MMM1D:
    MMM1D_init();
    break;
  case COULOMB_MMM2D:
    MMM2D_init();
    break;
  case COULOMB_MAGGS: 
    maggs_init();
    break;
  }
#endif

#ifdef DIPOLES
  switch(coulomb.Dmethod) {
#ifdef DP3M
  case DIPOLAR_MDLC_P3M:
    // fall through
  case DIPOLAR_P3M:
    dp3m_scaleby_box_l();
    break;
#endif
  }
#endif

#ifdef LB
  if(lattice_switch & LATTICE_LB) {
    lb_init();
#ifdef LB_BOUNDARIES
    lb_init_boundaries();
#endif
  }
#endif
#ifdef LB_GPU
  if(this_node == 0) {
    if(lattice_switch & LATTICE_LB_GPU) {
      lb_init_gpu();
#ifdef LB_BOUNDARIES_GPU
    lb_init_boundaries();
#endif
    }
  }
#endif
}

void on_cell_structure_change()
{
  EVENT_TRACE(fprintf(stderr, "%d: on_cell_structure_change\n", this_node));

  /* Now give methods a chance to react to the change in cell
     structure.  Most ES methods need to reinitialize, as they depend
     on skin, node grid and so on. Only for a change in box length we
     have separate, faster methods, as this might happend frequently
     in a NpT simulation. */
#ifdef ELECTROSTATICS
  switch (coulomb.method) {
  case COULOMB_DH:
    break;    
#ifdef P3M
  case COULOMB_ELC_P3M:
    ELC_init();
    // fall through
  case COULOMB_P3M:
    p3m_init();
    break;
#endif
  case COULOMB_MMM1D:
    MMM1D_init();
    break;
  case COULOMB_MMM2D:
    MMM2D_init();
    break;
  case COULOMB_MAGGS: 
    maggs_init();
    /* Maggs electrostatics needs ghost velocities */
    on_ghost_flags_change();
    break;
  default: break;
  }
#endif  /* ifdef ELECTROSTATICS */

#ifdef DIPOLES
  switch (coulomb.Dmethod) {
#ifdef DP3M
  case DIPOLAR_MDLC_P3M:
    // fall through
  case DIPOLAR_P3M:
    dp3m_init();
    break;
#endif
  default: break;
  }
#endif  /* ifdef DIPOLES */

#ifdef LB
  if(lattice_switch & LATTICE_LB) {
    lb_init();
  }
#endif
}

void on_temperature_change()
{
  EVENT_TRACE(fprintf(stderr, "%d: on_temperature_change\n", this_node));

#ifdef ELECTROSTATICS

#endif
#ifdef LB
  if (lattice_switch & LATTICE_LB) {
    lb_reinit_parameters();
  }
#endif
#ifdef LB_GPU
  if(this_node == 0) {
    if (lattice_switch & LATTICE_LB_GPU) {
      lb_reinit_parameters_gpu();
    }
  }
#endif
}

void on_parameter_change(int field)
{
  EVENT_TRACE(fprintf(stderr, "%d: on_parameter_change %s\n", this_node, fields[field].name));

  switch (field) {
  case FIELD_BOXL:
    grid_changed_box_l();
    /* Electrostatics cutoffs mostly depend on the system size,
       therefore recalculate them. */
    recalc_maximal_cutoff();
    cells_on_geometry_change(0);
    break;
  case FIELD_MIN_GLOBAL_CUT:
    recalc_maximal_cutoff();
    cells_on_geometry_change(0);
    break;
  case FIELD_SKIN:
    cells_on_geometry_change(0);
  case FIELD_PERIODIC:
    cells_on_geometry_change(CELL_FLAG_GRIDCHANGED);
    break;
  case FIELD_NODEGRID:
    grid_changed_n_nodes();
    cells_on_geometry_change(CELL_FLAG_GRIDCHANGED);
    break;
  case FIELD_MINNUMCELLS:
  case FIELD_MAXNUMCELLS:
    cells_re_init(CELL_STRUCTURE_CURRENT);
  case FIELD_TEMPERATURE:
    on_temperature_change();
    reinit_thermo = 1;
    break;
  case FIELD_TIMESTEP:
#ifdef LB_GPU
    if(this_node == 0) {
      if (lattice_switch & LATTICE_LB_GPU) {
        lb_reinit_parameters_gpu();
      }
    }  
#endif    
#ifdef LB
    if (lattice_switch & LATTICE_LB) {
      lb_reinit_parameters();
    }
#endif
  case FIELD_LANGEVIN_GAMMA:
  case FIELD_DPD_TGAMMA:
  case FIELD_DPD_GAMMA:
  case FIELD_NPTISO_G0:
  case FIELD_NPTISO_GV:
  case FIELD_NPTISO_PISTON:
    reinit_thermo = 1;
    break;
#ifdef NPT
  case FIELD_INTEG_SWITCH:
    if (integ_switch != INTEG_METHOD_NPT_ISO)
      nptiso.invalidate_p_vel = 1;
    break;
#endif
  case FIELD_THERMO_SWITCH:
    /* DPD needs ghost velocities, other thermostats not */
    on_ghost_flags_change();
    break;
#ifdef LB
  case FIELD_LATTICE_SWITCH:
    /* LB needs ghost velocities */
    on_ghost_flags_change();
    break;
#endif
  }
}

#ifdef LB
void on_lb_params_change(int field) {
  EVENT_TRACE(fprintf(stderr, "%d: on_lb_params_change\n", this_node));

  if (field == LBPAR_AGRID) {
    lb_init();
  }
  if (field == LBPAR_DENSITY) {
    lb_reinit_fluid();
  }

  lb_reinit_parameters();

}
#endif

#if defined (LB) || defined (LB_GPU)
void on_lb_params_change_gpu(int field) {
  EVENT_TRACE(fprintf(stderr, "%d: on_lb_params_change_gpu\n", this_node));

#ifdef LB_GPU
  if (field == LBPAR_AGRID) {
    lb_init_gpu();
#ifdef LB_BOUNDARIES_GPU
    lb_init_boundaries();
#endif
  }
  if (field == LBPAR_DENSITY) {
    lb_reinit_fluid_gpu();
  }

  lb_reinit_parameters_gpu();
#endif
}
#endif

void on_ghost_flags_change()
{
  EVENT_TRACE(fprintf(stderr, "%d: on_ghost_flags_change\n", this_node));
  /* that's all we change here */
  extern int ghosts_have_v;

  int old_have_v = ghosts_have_v;

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

  if (old_have_v != ghosts_have_v)
    cells_re_init(CELL_STRUCTURE_CURRENT);    
}

