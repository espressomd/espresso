// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2005; all rights reserved unless otherwise stated.
/** \file communication.c
    Implementation of \ref communication.h "communication.h".
*/
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "communication.h"
#include "interaction_data.h"
#include "particle_data.h"
#include "integrate.h"
#include "cells.h"
#include "global.h"
#include "grid.h"
#include "initialize.h"
#include "forces.h"
#include "rotation.h"
#include "p3m.h"
#include "statistics.h"
#include "energy.h"
#include "pressure.h"
#include "random.h"
#include "lj.h"
#include "lb.h"
#include "morse.h"
#include "buckingham.h"
#include "tab.h"
#include "ljcos.h"
#include "gb.h"
#include "mmm1d.h"
#include "mmm2d.h"
#include "maggs.h"
#include "elc.h"
#include "statistics_chain.h"
#include "topology.h"
#include "errorhandling.h"

int this_node = -1;
int n_nodes = -1;

/**********************************************
 * slave callbacks.
 **********************************************/

/** slave callback procedure */
typedef void (SlaveCallback)(int node, int param);

/** \name request codes in asynchronous mode
    <em> Keep this in sync with \ref communication::callbacks and \ref communication::names. </em> */
/*@{*/
/** Action number for \ref mpi_stop. */
#define REQ_TERM      0
/** Action number for \ref mpi_bcast_parameter. */
#define REQ_BCAST_PAR 1
/** Action number for \ref mpi_who_has. */
#define REQ_WHO_HAS   2
/** Action number for \ref mpi_bcast_event. */
#define REQ_EVENT     3
/** Action number for \ref mpi_place_particle. */
#define REQ_PLACE     4

/** Action number for \ref mpi_send_v. */
#define REQ_SET_V     5
/** Action number for \ref mpi_send_f. */
#define REQ_SET_F     6
/** Action number for \ref mpi_send_q. */
#define REQ_SET_Q     7
/** Action number for \ref mpi_send_type. */
#define REQ_SET_TYPE  8
/** Action number for \ref mpi_send_bond. */
#define REQ_SET_BOND  9

/** Action number for \ref mpi_recv_part. */
#define REQ_GET_PART  10
/** Action number for \ref mpi_integrate. */
#define REQ_INTEGRATE 11
/** Action number for \ref mpi_bcast_ia_params. */
#define REQ_BCAST_IA  12
/** Action number for \ref mpi_bcast_n_particle_types. */
#define REQ_BCAST_IA_SIZE  13
/** Action number for \ref mpi_gather_stats. */
#define REQ_GATHER    14

/** Action number for \ref mpi_set_time_step. */
#define REQ_SET_TIME_STEP  15
/** Action number for \ref mpi_get_particles. */
#define REQ_GETPARTS  16
/** Action number for \ref mpi_bcast_coulomb_params. */
#define REQ_BCAST_COULOMB 17
/** Action number for \ref mpi_send_ext. */
#define REQ_SET_EXT 18
/** Action number for \ref mpi_place_particle. */
#define REQ_PLACE_NEW   19

/** Action number for \ref mpi_remove_particle */
#define REQ_REM_PART   20
/** Action number for \ref mpi_bcast_constraint */
#define REQ_BCAST_CONSTR 21
/** Action number for \ref mpi_random_seed */
#define REQ_RANDOM_SEED 22
/** Action number for \ref mpi_random_stat */
#define REQ_RANDOM_STAT 23
/** Action number for \ref mpi_lj_cap_forces. */
#define REQ_BCAST_LFC 24

/** Action number for \ref mpi_tab_cap_forces. */
#define REQ_BCAST_TFC 25
/** Action number for \ref mpi_random_seed */
#define REQ_BIT_RANDOM_SEED 26
/** Action number for \ref mpi_random_stat */
#define REQ_BIT_RANDOM_STAT 27
/** Action number for \ref mpi_get_constraint_force */
#define REQ_GET_CONSFOR 28
/** Action number for \ref mpi_rescale_particles */
#define REQ_RESCALE_PART 29

/** Action number for \ref mpi_bcast_cell_structure */
#define REQ_BCAST_CS  30
/** Action number for \ref mpi_send_quat. */
#define REQ_SET_QUAT  31
/** Action number for \ref mpi_send_omega. */
#define REQ_SET_OMEGA 32
/** Action number for \ref mpi_send_torque. */
#define REQ_SET_TORQUE 33
/** Action number for \ref mpi_send_mol_id. */
#define REQ_SET_MOLID 34

/** Action number for \ref mpi_bcast_nptiso_geom */
#define REQ_BCAST_NPTISO_GEOM 35
/** Action number for \ref mpi_update_mol_ids  */
#define REQ_UPDATE_MOL_IDS 36
/** Action number for \ref mpi_sync_topo_part_info  */
#define REQ_SYNC_TOPO 37
/** Action number for \ref mpi_send_mass. */
#define REQ_SET_MASS  38
/** Action number for \ref mpi_buck_cap_forces. */
#define REQ_BCAST_BFC 39

/** Action number for \ref mpi_gather_runtime_errors  */
#define REQ_GET_ERRS  40
/** Action number for \ref mpi_send_exclusion. */
#define REQ_SET_EXCL  41
/** Action number for \ref mpi_morse_cap_forces. */
#define REQ_BCAST_MFC 42
/** Action number for \ref lb_set_params. */
#define REQ_LB_BCAST 43
/** Total number of action numbers. */
#define REQ_MAXIMUM   44

/*@}*/

/** \name Slave Callbacks
    These functions are the slave node counterparts for the
    commands the master node issues. The master node function *
    corresponds to the slave node function *_slave.
*/
/*@{*/
void mpi_stop_slave(int node, int parm);
void mpi_bcast_parameter_slave(int node, int parm);
void mpi_who_has_slave(int node, int parm);
void mpi_bcast_event_slave(int node, int parm);
void mpi_place_particle_slave(int node, int parm);
void mpi_send_v_slave(int node, int parm);
void mpi_send_f_slave(int node, int parm);
void mpi_send_q_slave(int node, int parm);
void mpi_send_type_slave(int node, int parm);
void mpi_send_bond_slave(int node, int parm);
void mpi_recv_part_slave(int node, int parm);
void mpi_integrate_slave(int node, int parm);
void mpi_bcast_ia_params_slave(int node, int parm);
void mpi_bcast_n_particle_types_slave(int node, int parm);
void mpi_gather_stats_slave(int node, int parm);
void mpi_set_time_step_slave(int node, int parm);
void mpi_get_particles_slave(int node, int parm);
void mpi_bcast_coulomb_params_slave(int node, int parm);
void mpi_send_ext_slave(int node, int parm);
void mpi_remove_particle_slave(int node, int parm);
void mpi_bcast_constraint_slave(int node, int parm);
void mpi_random_seed_slave(int node, int parm);
void mpi_random_stat_slave(int node, int parm);
void mpi_lj_cap_forces_slave(int node, int parm);
void mpi_tab_cap_forces_slave(int node, int parm);
void mpi_bit_random_seed_slave(int node, int parm);
void mpi_bit_random_stat_slave(int node, int parm);
void mpi_get_constraint_force_slave(int node, int parm);
void mpi_rescale_particles_slave(int node, int parm);
void mpi_bcast_cell_structure_slave(int node, int parm);
void mpi_send_quat_slave(int node, int parm);
void mpi_send_omega_slave(int node, int parm);
void mpi_send_torque_slave(int node, int parm);
void mpi_send_mol_id_slave(int node, int parm);
void mpi_bcast_nptiso_geom_slave(int node, int parm);
void mpi_update_mol_ids_slave(int node, int parm);
void mpi_sync_topo_part_info_slave(int node, int parm);
void mpi_send_mass_slave(int node, int parm);
void mpi_buck_cap_forces_slave(int node, int parm);
void mpi_gather_runtime_errors_slave(int node, int parm);
void mpi_send_exclusion_slave(int node, int parm);
void mpi_morse_cap_forces_slave(int node, int parm);
void mpi_bcast_lb_params_slave(int node, int parm);

/*@}*/

/** A list of which function has to be called for
    the issued command. */
static SlaveCallback *slave_callbacks[] = {
  mpi_stop_slave,                   /*  0: REQ_TERM */
  mpi_bcast_parameter_slave,        /*  1: REQ_BCAST_PAR */
  mpi_who_has_slave,                /*  2: REQ_WHO_HAS */
  mpi_bcast_event_slave,            /*  3: REQ_EVENT */
  mpi_place_particle_slave,         /*  4: REQ_SET_POS */
  mpi_send_v_slave,                 /*  5: REQ_SET_V */
  mpi_send_f_slave,                 /*  6: REQ_SET_F */
  mpi_send_q_slave,                 /*  7: REQ_SET_Q */
  mpi_send_type_slave,              /*  8: REQ_SET_TYPE */
  mpi_send_bond_slave,              /*  9: REQ_SET_BOND */
  mpi_recv_part_slave,              /* 10: REQ_GET_PART */
  mpi_integrate_slave,              /* 11: REQ_INTEGRATE */
  mpi_bcast_ia_params_slave,        /* 12: REQ_BCAST_IA */
  mpi_bcast_n_particle_types_slave, /* 13: REQ_BCAST_IA_SIZE */
  mpi_gather_stats_slave,           /* 14: REQ_GATHER */
  mpi_set_time_step_slave,          /* 15: REQ_SET_TIME_STEP */
  mpi_get_particles_slave,          /* 16: REQ_GETPARTS */
  mpi_bcast_coulomb_params_slave,   /* 17: REQ_BCAST_COULOMB */
  mpi_send_ext_slave,               /* 18: REQ_SEND_EXT */
  mpi_place_particle_slave,         /* 19: REQ_PLACE_NEW */
  mpi_remove_particle_slave,        /* 20: REQ_REM_PART */
  mpi_bcast_constraint_slave,       /* 21: REQ_BCAST_CONSTR */
  mpi_random_seed_slave,            /* 22: REQ_RANDOM_SEED */
  mpi_random_stat_slave,            /* 23: REQ_RANDOM_STAT */
  mpi_lj_cap_forces_slave,          /* 24: REQ_BCAST_LFC */
  mpi_tab_cap_forces_slave,         /* 25: REQ_BCAST_TFC */
  mpi_bit_random_seed_slave,        /* 26: REQ_RANDOM_SEED */
  mpi_bit_random_stat_slave,        /* 27: REQ_RANDOM_STAT */
  mpi_get_constraint_force_slave,   /* 28: REQ_GET_CONSFOR */
  mpi_rescale_particles_slave,      /* 29: REQ_RESCALE_PART */
  mpi_bcast_cell_structure_slave,   /* 30: REQ_BCAST_CS */
  mpi_send_quat_slave,              /* 31: REQ_SET_QUAT */
  mpi_send_omega_slave,             /* 32: REQ_SET_OMEGA */
  mpi_send_torque_slave,            /* 33: REQ_SET_TORQUE */
  mpi_send_mol_id_slave,            /* 34: REQ_SET_MOLID */
  mpi_bcast_nptiso_geom_slave,      /* 35: REQ_BCAST_NPTISO_GEOM */
  mpi_update_mol_ids_slave,         /* 36: REQ_UPDATE_MOL_IDS */
  mpi_sync_topo_part_info_slave,    /* 37: REQ_SYNC_TOPO */
  mpi_send_mass_slave,              /* 38: REQ_SET_MASS */
  mpi_buck_cap_forces_slave,        /* 39: REQ_BCAST_LFC */
  mpi_gather_runtime_errors_slave,  /* 40: REQ_GET_ERRS */
  mpi_send_exclusion_slave,         /* 41: REQ_SET_EXCL */
  mpi_morse_cap_forces_slave,       /* 42: REQ_BCAST_MFC */
  mpi_bcast_lb_params_slave,        /* 43: REQ_LB_BCAST */
};

/** Names to be printed when communication debugging is on. */
char *names[] = {
  "TERM"      ,     /*  0 */
  "BCAST_PAR" ,     /*  1 */
  "WHO_HAS"   ,     /*  2 */
  "EVENT"     ,     /*  3 */
  "SET_POS"   ,     /*  4 */

  "SET_V"     ,     /*  5 */
  "SET_F"     ,     /*  6 */
  "SET_Q"     ,     /*  7 */
  "SET_TYPE"  ,     /*  8 */
  "SET_BOND"  ,     /*  9 */

  "GET_PART"  ,     /* 10 */
  "INTEGRATE" ,     /* 11 */
  "BCAST_IA"  ,     /* 12 */
  "BCAST_IAS" ,     /* 13 */
  "GATHER"    ,     /* 14 */

  "TIME_STEP" ,     /* 15 */
  "GET_PARTS" ,     /* 16 */
  "BCAST_CIA" ,     /* 17 */
  "SEND_EXT"  ,     /* 18 */
  "PLACE_NEW" ,     /* 19 */

  "REM_PART"  ,     /* 20 */
  "BCAST_CON" ,     /* 21 */
  "RAND_SEED" ,     /* 22 */
  "RAND_STAT" ,     /* 23 */
  "BCAST_LFC" ,     /* 24 */

  "BCAST_TFC" ,     /* 25 */
  "BIT_RAND_SEED",  /* 26 */
  "BIT_RAND_STAT",  /* 27 */
  "GET_CONSTR",     /* 28 */
  "RESCALE_PART",   /* 29 */

  "BCAST_CS",       /* 30 */
  "SET_QUAT"  ,     /* 31 */
  "SET_OMEGA",      /* 32 */
  "SET_TORQUE",     /* 33 */
  "SET_MOLID",      /* 34 */

  "BCAST_NPT_GEOM", /* 35 */
  "UPDATE_MOL_IDS", /* 36 */
  "SYNC_TOPO",      /* 37 */
  "SET_MASS",       /* 38 */
  "BCAST_BFC" ,     /* 39 */

  "GET_ERRS",       /* 40 */
  "SET_EXCL",       /* 41 */
  "BCAST_MFC" ,     /* 42 */
  "BCAST_LB",       /* 43 */
};

/** the requests are compiled here. So after a crash you get the last issued request */
static int request[3];

/**********************************************
 * procedures
 **********************************************/

#ifdef MPI_CORE
void mpi_core(MPI_Comm *comm, int *errcode,...) {
  fprintf(stderr, "%d: Aborting due to MPI error %d, forcing core dump\n", this_node, *errcode);
  fflush(stderr);
  core();
}
#endif

void mpi_init(int *argc, char ***argv)
{
#ifdef MPI_CORE
  MPI_Errhandler mpi_errh;
#endif

  MPI_Init(argc, argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &this_node);
  MPI_Comm_size(MPI_COMM_WORLD, &n_nodes);

#ifdef MPI_CORE
  MPI_Errhandler_create((MPI_Handler_function *)mpi_core, &mpi_errh);
  MPI_Errhandler_set(MPI_COMM_WORLD, mpi_errh);
#endif

}

static void mpi_issue(int reqcode, int node, int param)
{
  request[0] = reqcode;
  request[1] = node;
  request[2] = param;

  COMM_TRACE(fprintf(stderr, "0: issuing %s(%d), assigned to node %d\n",
		     names[reqcode], param, node));
#ifdef ASYNC_BARRIER
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  MPI_Bcast(request, 3, MPI_INT, 0, MPI_COMM_WORLD);
}

/**************** REQ_TERM ************/

static int terminated = 0;

void mpi_stop()
{
  if (terminated)
    return;

  mpi_issue(REQ_TERM, -1, 0);

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  regular_exit = 1;
  terminated = 1;
}

void mpi_stop_slave(int node, int param)
{
  COMM_TRACE(fprintf(stderr, "%d: exiting\n", this_node));

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  regular_exit = 1;
  exit(0);
}

/*************** REQ_BCAST_PAR ************/
static void common_bcast_parameter(int i)
{
  switch (fields[i].type) {
  case TYPE_INT:
    MPI_Bcast((int *)fields[i].data, fields[i].dimension,
	      MPI_INT, 0, MPI_COMM_WORLD);
    break;
  case TYPE_BOOL:
    MPI_Bcast((int *)fields[i].data, 1,
	      MPI_INT, 0, MPI_COMM_WORLD);
    break;
  case TYPE_DOUBLE:
    MPI_Bcast((double *)fields[i].data, fields[i].dimension,
	      MPI_DOUBLE, 0, MPI_COMM_WORLD);
    break;
  default: break;
  }

  on_parameter_change(i);
}

int mpi_bcast_parameter(int i)
{
  mpi_issue(REQ_BCAST_PAR, -1, i);

  common_bcast_parameter(i);

  return check_runtime_errors();
}

void mpi_bcast_parameter_slave(int node, int i)
{
  common_bcast_parameter(i);
  check_runtime_errors();
}

/*************** REQ_WHO_HAS ****************/

void mpi_who_has()
{
  Cell *cell;
  int *sizes = malloc(sizeof(int)*n_nodes);
  int *pdata = NULL;
  int pdata_s = 0, i, c;
  int pnode;
  int n_part;
  MPI_Status status;

  mpi_issue(REQ_WHO_HAS, -1, 0);

  n_part = cells_get_n_particles();
  /* first collect number of particles on each node */
  MPI_Gather(&n_part, 1, MPI_INT, sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);

  for (i = 0; i <= max_seen_particle; i++)
    particle_node[i] = -1;

  /* then fetch particle locations */
  for (pnode = 0; pnode < n_nodes; pnode++) {
    COMM_TRACE(fprintf(stderr, "node %d reports %d particles\n",
		       pnode, sizes[pnode]));
    if (pnode == this_node) {
      for (c = 0; c < local_cells.n; c++) {
	cell = local_cells.cell[c];
	for (i = 0; i < cell->n; i++)
	  particle_node[cell->part[i].p.identity] = this_node;
      }
    }
    else if (sizes[pnode] > 0) {
      if (pdata_s < sizes[pnode]) {
	pdata_s = sizes[pnode];
	pdata = (int *)realloc(pdata, sizeof(int)*pdata_s);
      }
      MPI_Recv(pdata, sizes[pnode], MPI_INT, pnode, REQ_WHO_HAS,
	       MPI_COMM_WORLD, &status);
      for (i = 0; i < sizes[pnode]; i++)
	particle_node[pdata[i]] = pnode;
    }
  }

  free(pdata);
  free(sizes);
}

void mpi_who_has_slave(int node, int param)
{
  Cell *cell;
  int npart, i, c;
  int *sendbuf;
  int n_part;

  n_part = cells_get_n_particles();
  MPI_Gather(&n_part, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
  if (n_part == 0)
    return;

  sendbuf = malloc(sizeof(int)*n_part);
  npart = 0;
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    for (i = 0; i < cell->n; i++)
      sendbuf[npart++] = cell->part[i].p.identity;
  }
  MPI_Send(sendbuf, npart, MPI_INT, 0, REQ_WHO_HAS, MPI_COMM_WORLD);
  free(sendbuf);
}

/**************** REQ_CHTOPL ***********/
void mpi_bcast_event(int event)
{
  mpi_issue(REQ_EVENT, -1, event);
  mpi_bcast_event_slave(-1, event);
}

void mpi_bcast_event_slave(int node, int event)
{
  switch (event) {
#ifdef ELECTROSTATICS
  case P3M_COUNT_CHARGES:
    P3M_count_charged_particles();
    break; 
  case MAGGS_COUNT_CHARGES:
    maggs_count_charged_particles();
    break; 
#endif
  case INVALIDATE_SYSTEM:
    local_invalidate_system();
    break;
#ifdef ADDITIONAL_CHECKS
  case CHECK_PARTICLES:
    check_particles();
    break;
#endif
  default:;
  }
}

/****************** REQ_PLACE/REQ_PLACE_NEW ************/

void mpi_place_particle(int pnode, int part, int new, double p[3])
{
  if (new) {
    mpi_issue(REQ_PLACE_NEW, pnode, part);
    added_particle(part);
  }
  else
    mpi_issue(REQ_PLACE, pnode, part);

  if (pnode == this_node)
    local_place_particle(part, p, new);
  else
    MPI_Send(p, 3, MPI_DOUBLE, pnode, REQ_PLACE, MPI_COMM_WORLD);

  on_particle_change();
}

void mpi_place_particle_slave(int pnode, int part)
{
  double p[3];
  MPI_Status status;
  int new = (request[0] == REQ_PLACE_NEW); 

  if (new)
    added_particle(part);

  if (pnode == this_node) {
    MPI_Recv(p, 3, MPI_DOUBLE, 0, REQ_PLACE, MPI_COMM_WORLD, &status);
    local_place_particle(part, p, new);
  }

  on_particle_change();
}

/****************** REQ_SET_V ************/
void mpi_send_v(int pnode, int part, double v[3])
{
  mpi_issue(REQ_SET_V, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    memcpy(p->m.v, v, 3*sizeof(double));
  }
  else
    MPI_Send(v, 3, MPI_DOUBLE, pnode, REQ_SET_V, MPI_COMM_WORLD);

  on_particle_change();
}

void mpi_send_v_slave(int pnode, int part)
{
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Status status;
    MPI_Recv(p->m.v, 3, MPI_DOUBLE, 0, REQ_SET_V,
	     MPI_COMM_WORLD, &status);
  }

  on_particle_change();
}

/****************** REQ_SET_F ************/
void mpi_send_f(int pnode, int part, double F[3])
{
  mpi_issue(REQ_SET_F, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    memcpy(p->f.f, F, 3*sizeof(double));
  }
  else
    MPI_Send(F, 3, MPI_DOUBLE, pnode, REQ_SET_F, MPI_COMM_WORLD);

  on_particle_change();
}

void mpi_send_f_slave(int pnode, int part)
{
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Status status;
    MPI_Recv(p->f.f, 3, MPI_DOUBLE, 0, REQ_SET_F,
	     MPI_COMM_WORLD, &status);
  }

  on_particle_change();
}

/********************* REQ_SET_Q ********/
void mpi_send_q(int pnode, int part, double q)
{
#ifdef ELECTROSTATICS
  mpi_issue(REQ_SET_Q, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->p.q = q;
  }
  else {
    MPI_Send(&q, 1, MPI_DOUBLE, pnode, REQ_SET_Q, MPI_COMM_WORLD);
  }

  on_particle_change();
#endif
}

void mpi_send_q_slave(int pnode, int part)
{
#ifdef ELECTROSTATICS
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Status status;
    MPI_Recv(&p->p.q, 1, MPI_DOUBLE, 0, REQ_SET_Q,
	     MPI_COMM_WORLD, &status);
  }

  on_particle_change();
#endif
}


/********************* REQ_SET_M ********/
void mpi_send_mass(int pnode, int part, double mass)
{
#ifdef MASS
  mpi_issue(REQ_SET_MASS, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->p.mass = mass;
  }
  else {
    MPI_Send(&mass, 1, MPI_DOUBLE, pnode, REQ_SET_MASS, MPI_COMM_WORLD);
  }

  on_particle_change();
#endif
}

void mpi_send_mass_slave(int pnode, int part)
{
#ifdef MASS
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Status status;
    MPI_Recv(&p->p.mass, 1, MPI_DOUBLE, 0, REQ_SET_MASS,
	     MPI_COMM_WORLD, &status);
  }

  on_particle_change();
#endif
}

/********************* REQ_SET_TYPE ********/
void mpi_send_type(int pnode, int part, int type)
{
  mpi_issue(REQ_SET_TYPE, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->p.type = type;
  }
  else
    MPI_Send(&type, 1, MPI_INT, pnode, REQ_SET_TYPE, MPI_COMM_WORLD);

  on_particle_change();
}

void mpi_send_type_slave(int pnode, int part)
{
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Status status;
    MPI_Recv(&p->p.type, 1, MPI_INT, 0, REQ_SET_TYPE,
	     MPI_COMM_WORLD, &status);
  }

  on_particle_change();
}

/********************* REQ_SET_MOLID ********/
void mpi_send_mol_id(int pnode, int part, int mid)
{
  mpi_issue(REQ_SET_MOLID, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->p.mol_id = mid;
  }
  else
    MPI_Send(&mid, 1, MPI_INT, pnode, REQ_SET_MOLID, MPI_COMM_WORLD);

  on_particle_change();
}

void mpi_send_mol_id_slave(int pnode, int part)
{
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Status status;
    MPI_Recv(&p->p.mol_id, 1, MPI_INT, 0, REQ_SET_MOLID,
	     MPI_COMM_WORLD, &status);
  }

  on_particle_change();
}

/********************* REQ_SET_QUAT ********/

void mpi_send_quat(int pnode, int part, double quat[4])
{
#ifdef ROTATION
  mpi_issue(REQ_SET_QUAT, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->r.quat[0] = quat[0];
    p->r.quat[1] = quat[1];
    p->r.quat[2] = quat[2];
    p->r.quat[3] = quat[3];    
  }
  else {
    MPI_Send(quat, 4, MPI_DOUBLE, pnode, REQ_SET_QUAT, MPI_COMM_WORLD);
  }

  on_particle_change();
#endif
}

void mpi_send_quat_slave(int pnode, int part)
{
#ifdef ROTATION
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Status status;
    MPI_Recv(p->r.quat, 4, MPI_DOUBLE, 0, REQ_SET_QUAT,
	     MPI_COMM_WORLD, &status);
  }

  on_particle_change();
#endif
}

/********************* REQ_SET_OMEGA ********/

void mpi_send_omega(int pnode, int part, double omega[3])
{
#ifdef ROTATION
  mpi_issue(REQ_SET_OMEGA, pnode, part);

  if (pnode == this_node) {
   Particle *p = local_particles[part];
/*  memcpy(p->omega, omega, 3*sizeof(double));*/
    p->m.omega[0] = omega[0];
    p->m.omega[1] = omega[1];
    p->m.omega[2] = omega[2];
  }
  else {
    MPI_Send(omega, 3, MPI_DOUBLE, pnode, REQ_SET_OMEGA, MPI_COMM_WORLD);
  }

  on_particle_change();
#endif
}

void mpi_send_omega_slave(int pnode, int part)
{
#ifdef ROTATION
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Status status;
    MPI_Recv(p->m.omega, 3, MPI_DOUBLE, 0, REQ_SET_OMEGA,
	     MPI_COMM_WORLD, &status);
  }

  on_particle_change();
#endif
}

/********************* REQ_SET_TORQUE ********/

void mpi_send_torque(int pnode, int part, double torque[3])
{
#ifdef ROTATION
  mpi_issue(REQ_SET_TORQUE, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->f.torque[0] = torque[0];
    p->f.torque[1] = torque[1];
    p->f.torque[2] = torque[2];
  }
  else {
    MPI_Send(torque, 3, MPI_DOUBLE, pnode, REQ_SET_TORQUE, MPI_COMM_WORLD);
  }

  on_particle_change();
#endif
}

void mpi_send_torque_slave(int pnode, int part)
{
#ifdef ROTATION
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Status status;
    MPI_Recv(p->f.torque, 3, MPI_DOUBLE, 0, REQ_SET_TORQUE,
	     MPI_COMM_WORLD, &status);
  }

  on_particle_change();
#endif
}


/********************* REQ_SET_BOND ********/
int mpi_send_bond(int pnode, int part, int *bond, int delete)
{
  int bond_size, stat;
  MPI_Status status;

  mpi_issue(REQ_SET_BOND, pnode, part);

  bond_size = (bond) ? bonded_ia_params[bond[0]].num + 1 : 0;

  if (pnode == this_node) {
    stat = local_change_bond(part, bond, delete);
    on_particle_change();
    return stat;
  }
  /* else */
  MPI_Send(&bond_size, 1, MPI_INT, pnode, REQ_SET_BOND, MPI_COMM_WORLD);
  if (bond_size)
    MPI_Send(bond, bond_size, MPI_INT, pnode, REQ_SET_BOND, MPI_COMM_WORLD);
  MPI_Send(&delete, 1, MPI_INT, pnode, REQ_SET_BOND, MPI_COMM_WORLD);
  MPI_Recv(&stat, 1, MPI_INT, pnode, REQ_SET_BOND, MPI_COMM_WORLD, &status);
  on_particle_change();
  return stat;
}

void mpi_send_bond_slave(int pnode, int part)
{
  int bond_size, *bond, delete, stat;
  MPI_Status status;

  if (pnode == this_node) {
    MPI_Recv(&bond_size, 1, MPI_INT, 0, REQ_SET_BOND, MPI_COMM_WORLD, &status);
    if (bond_size) {
      bond = (int *)malloc(bond_size*sizeof(int));
      MPI_Recv(bond, bond_size, MPI_INT, 0, REQ_SET_BOND, MPI_COMM_WORLD, &status);
    }
    else
      bond = NULL;
    MPI_Recv(&delete, 1, MPI_INT, 0, REQ_SET_BOND, MPI_COMM_WORLD, &status);
    stat = local_change_bond(part, bond, delete);
    if (bond)
      free(bond);
    MPI_Send(&stat, 1, MPI_INT, 0, REQ_SET_BOND, MPI_COMM_WORLD);
  }

  on_particle_change();
}

/****************** REQ_GET_PART ************/
void mpi_recv_part(int pnode, int part, Particle *pdata)
{
  IntList *bl = &(pdata->bl);
#ifdef EXCLUSIONS
  IntList *el = &(pdata->el);
#endif
  MPI_Status status;

  /* fetch fixed data */
  if (pnode == this_node)
    memcpy(pdata, local_particles[part], sizeof(Particle));
  else {
    mpi_issue(REQ_GET_PART, pnode, part);
    MPI_Recv(pdata, sizeof(Particle), MPI_BYTE, pnode,
	     REQ_GET_PART, MPI_COMM_WORLD, &status);
  }

  /* copy dynamic data */
  /* bonds */
  bl->max = bl->n;
  if (bl->n > 0) {
    alloc_intlist(bl, bl->n);
    if (pnode == this_node)
      memcpy(bl->e, local_particles[part]->bl.e, sizeof(int)*bl->n);
    else
      MPI_Recv(bl->e, bl->n, MPI_INT, pnode,
               REQ_GET_PART, MPI_COMM_WORLD, &status);
  }
  else
    bl->e = NULL;

#ifdef EXCLUSIONS
  /* exclusions */
  el->max = el->n;
  if (el->n > 0) {
    alloc_intlist(el, el->n);
    if (pnode == this_node)
      memcpy(el->e, local_particles[part]->el.e, sizeof(int)*el->n);
    else
      MPI_Recv(el->e, el->n, MPI_INT, pnode,
               REQ_GET_PART, MPI_COMM_WORLD, &status);
  }
  else
    el->e = NULL;
#endif
}

void mpi_recv_part_slave(int pnode, int part)
{
  Particle *p;
  if (pnode != this_node)
    return;

  p = local_particles[part];

  MPI_Send(p, sizeof(Particle), MPI_BYTE, 0, REQ_GET_PART,
	   MPI_COMM_WORLD);
  if (p->bl.n > 0)
    MPI_Send(p->bl.e, p->bl.n, MPI_INT, 0, REQ_GET_PART,
	     MPI_COMM_WORLD);
#ifdef EXCLUSIONS
  if (p->el.n > 0)
    MPI_Send(p->el.e, p->el.n, MPI_INT, 0, REQ_GET_PART,
	     MPI_COMM_WORLD);
#endif
}

/****************** REQ_REM_PART ************/
void mpi_remove_particle(int pnode, int part)
{
  mpi_issue(REQ_REM_PART, pnode, part);
  mpi_remove_particle_slave(pnode, part);
}

void mpi_remove_particle_slave(int pnode, int part)
{
  if (part != -1) {
    n_total_particles--;

    if (pnode == this_node)
      local_remove_particle(part);
    
    remove_all_bonds_to(part);
  }
  else
    local_remove_all_particles();

  on_particle_change();
}

/********************* REQ_INTEGRATE ********/
int mpi_integrate(int n_steps)
{
  mpi_issue(REQ_INTEGRATE, -1, n_steps);

  integrate_vv(n_steps);
  COMM_TRACE(fprintf(stderr, "%d: integration task %d done.\n", this_node, n_steps));

  return check_runtime_errors();
}

void mpi_integrate_slave(int pnode, int task)
{
  integrate_vv(task);
  COMM_TRACE(fprintf(stderr, "%d: integration task %d done.\n", this_node, task));

  check_runtime_errors();
}

/*************** REQ_BCAST_IA ************/
void mpi_bcast_ia_params(int i, int j)
{
  int tablesize;

  mpi_issue(REQ_BCAST_IA, i, j);
  tablesize = tabulated_forces.max;

  if (j>=0) {
    /* non-bonded interaction parameters */
    /* INCOMPATIBLE WHEN NODES USE DIFFERENT ARCHITECTURES */
    MPI_Bcast(get_ia_param(i, j), sizeof(IA_parameters), MPI_BYTE,
	      0, MPI_COMM_WORLD);

#ifdef TABULATED
    /* If there are tabulated forces broadcast those as well */
    if ( get_ia_param(i,j)->TAB_maxval > 0) {
      /* First let all nodes know the new size for force and energy tables */
      MPI_Bcast(&tablesize, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD); // Don't do anything until all nodes have this information

      /* Communicate the data */
      MPI_Bcast(tabulated_forces.e,tablesize, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
      MPI_Bcast(tabulated_energies.e,tablesize, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
    }    
#endif
  }
  else {
    /* bonded interaction parameters */
    /* INCOMPATIBLE WHEN NODES USE DIFFERENT ARCHITECTURES */
    MPI_Bcast(&(bonded_ia_params[i]),sizeof(Bonded_ia_parameters), MPI_BYTE,
	      0, MPI_COMM_WORLD);
#ifdef TABULATED
    /* For tabulated potentials we have to send the tables extra */
    if(bonded_ia_params[i].type == BONDED_IA_TABULATED) {
      int size = bonded_ia_params[i].p.tab.npoints;
      MPI_Bcast(bonded_ia_params[i].p.tab.f, size, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
      MPI_Bcast(bonded_ia_params[i].p.tab.e, size, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
    }
#endif
  }

  on_short_range_ia_change();
}

void mpi_bcast_ia_params_slave(int i, int j)
{
  if(j >= 0) { /* non-bonded interaction parameters */
    /* INCOMPATIBLE WHEN NODES USE DIFFERENT ARCHITECTURES */
    MPI_Bcast(get_ia_param(i, j), sizeof(IA_parameters), MPI_BYTE,
	      0, MPI_COMM_WORLD);
#ifdef TABULATED
    {
      int tablesize;
      /* If there are tabulated forces broadcast those as well */
      if ( get_ia_param(i,j)->TAB_maxval > 0) {
	/* Determine the new size for force and energy tables */
	MPI_Bcast(&tablesize,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	/* Allocate sizes accordingly */
	realloc_doublelist(&tabulated_forces, tablesize);
	realloc_doublelist(&tabulated_energies, tablesize);
	/* Now communicate the data */
	MPI_Bcast(tabulated_forces.e,tablesize, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
	MPI_Bcast(tabulated_energies.e,tablesize, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
      }
    }
#endif
  }
  else { /* bonded interaction parameters */
    make_bond_type_exist(i); /* realloc bonded_ia_params on slave nodes! */
    MPI_Bcast(&(bonded_ia_params[i]),sizeof(Bonded_ia_parameters), MPI_BYTE,
	      0, MPI_COMM_WORLD);
#ifdef TABULATED
    /* For tabulated potentials we have to send the tables extra */
    if(bonded_ia_params[i].type == BONDED_IA_TABULATED) {
      int size = bonded_ia_params[i].p.tab.npoints;
      /* alloc force and energy tables on slave nodes! */
      bonded_ia_params[i].p.tab.f = (double*)malloc(size*sizeof(double));
      bonded_ia_params[i].p.tab.e = (double*)malloc(size*sizeof(double));
      MPI_Bcast(bonded_ia_params[i].p.tab.f, size, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
      MPI_Bcast(bonded_ia_params[i].p.tab.e, size, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
    }
#endif
  }

  on_short_range_ia_change();
}

/*************** REQ_BCAST_IA_SIZE ************/

void mpi_bcast_n_particle_types(int ns)
{
  mpi_issue(REQ_BCAST_IA_SIZE, -1, ns);
  mpi_bcast_n_particle_types_slave(-1, ns);
}

void mpi_bcast_n_particle_types_slave(int pnode, int ns)
{
  realloc_ia_params(ns);
}

/*************** REQ_GATHER ************/
void mpi_gather_stats(int job, void *result)
{
  switch (job) {
  case 1:
    mpi_issue(REQ_GATHER, -1, 1);
    energy_calc(result);
    break;
  case 2:
    mpi_issue(REQ_GATHER, -1, 2);
    pressure_calc(result,0);
    break;
  case 3:
    mpi_issue(REQ_GATHER, -1, 3);
    pressure_calc(result,1);
    break;
  default:
    fprintf(stderr, "%d: INTERNAL ERROR: illegal request %d for REQ_GATHER\n", this_node, job);
    errexit();
  }
}

void mpi_gather_stats_slave(int ana_num, int job)
{
  switch (job) {
  case 1:
    /* calculate and reduce (sum up) energies */
    energy_calc(NULL);
    break;
  case 2:
    /* calculate and reduce (sum up) virials for 'analyze pressure' */
    pressure_calc(NULL,0);
    break;
  case 3:
    /* calculate and reduce (sum up) virials, revert velocities half a timestep for 'analyze p_inst' */
    pressure_calc(NULL,1);
    break;
  default:
    fprintf(stderr, "%d: INTERNAL ERROR: illegal request %d for REQ_GATHER\n", this_node, job);
    errexit();
  }
}

/*************** REQ_GETPARTS ************/
void mpi_get_particles(Particle *result, IntList *bi)
{
  IntList local_bi;
  int n_part;
  int tot_size, i, g, pnode;
  int *sizes;
  Cell *cell;
  int c;
  MPI_Status status;

  mpi_issue(REQ_GETPARTS, -1, bi != NULL);

  sizes = malloc(sizeof(int)*n_nodes);
  n_part = cells_get_n_particles();

  /* first collect number of particles on each node */
  MPI_Gather(&n_part, 1, MPI_INT, sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
  tot_size = 0;
  for (i = 0; i < n_nodes; i++)
    tot_size += sizes[i];

  if (tot_size!=n_total_particles) {
    fprintf(stderr,"%d: ERROR: mpi_get_particles: n_total_particles %d, but I counted %d. Exiting...\n",
	    this_node, n_total_particles, tot_size);
    errexit();
  }

  /* fetch particle informations into 'result' */
  init_intlist(&local_bi);
  g = 0;
  for (pnode = 0; pnode < n_nodes; pnode++) {
    if (sizes[pnode] > 0) {
      if (pnode == this_node) {
	for (c = 0; c < local_cells.n; c++) {
	  Particle *part;
	  int npart;
	  cell = local_cells.cell[c];
	  part = cell->part;
	  npart = cell->n;
	  memcpy(&result[g], part, npart*sizeof(Particle));
	  g += npart;
	  if (bi) {
	    int pc;
	    for (pc = 0; pc < npart; pc++) {
	      Particle *p = &part[pc];
	      realloc_intlist(&local_bi, local_bi.n + p->bl.n);
	      memcpy(&local_bi.e[local_bi.n], p->bl.e, p->bl.n*sizeof(int));
	      local_bi.n += p->bl.n;
	    }
	  }
	}
      }
      else {
	MPI_Recv(&result[g], sizes[pnode]*sizeof(Particle), MPI_BYTE, pnode, REQ_GETPARTS,
		 MPI_COMM_WORLD, &status);
	g += sizes[pnode];
      }
    }
  }

  /* perhaps add some debugging output */
#ifdef ELECTROSTATICS
  COMM_TRACE(for(i = 0; i < tot_size; i++) {
    printf("%d: %d -> %d %d %f (%f, %f, %f)\n", this_node, i, result[i].p.identity, result[i].p.type,
	   result[i].p.q, result[i].r.p[0], result[i].r.p[1], result[i].r.p[2]);
  });
#endif

  /* gather bonding information */
  if (bi) {
    int *bonds;

    init_intlist(bi);
    MPI_Gather(&local_bi.n, 1, MPI_INT, sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
    for (pnode = 0; pnode < n_nodes; pnode++) {
      if (sizes[pnode] > 0) {
	realloc_intlist(bi, bi->n + sizes[pnode]);

	if (pnode == this_node)
	  memcpy(&bi->e[bi->n], local_bi.e, sizes[pnode]*sizeof(int));
	else
	  MPI_Recv(&bi->e[bi->n], sizes[pnode], MPI_INT, pnode, REQ_GETPARTS,
		   MPI_COMM_WORLD, &status);

	bi->n += sizes[pnode];
      }
    }

    /* setup particle bond pointers into bi */
    bonds = bi->e;
    for (i = 0; i < tot_size; i++) {
      result[i].bl.e = bonds;
      bonds += result[i].bl.n;
      COMM_TRACE(if (result[i].bl.n > 0) {
	printf("(%d) part %d: bonds ", i, result[i].p.identity);
	for(g = 0; g < result[i].bl.n; g++) printf("%d ", result[i].bl.e[g]);
	printf("\n");
      });
    }
    realloc_intlist(&local_bi, 0);
  }

  COMM_TRACE(fprintf(stderr, "%d: finished\n", this_node));

  free(sizes);
}

void mpi_get_particles_slave(int pnode, int bi)
{
  int n_part;
  int g;
  Particle *result;
  Cell *cell;
  int c;

  n_part = cells_get_n_particles();

  COMM_TRACE(fprintf(stderr, "%d: get_particles_slave, %d particles\n", this_node, n_part));

  /* first collect number of particles on each node */
  MPI_Gather(&n_part, 1, MPI_INT, NULL, 1, MPI_INT,
	     0, MPI_COMM_WORLD);

  if (n_part > 0) {
    IntList local_bi;

    /* get (unsorted) particle informations as an array of type 'particle' */
    /* then get the particle information */
    result = malloc(n_part*sizeof(Particle));

    init_intlist(&local_bi);
    
    g = 0;
    for (c = 0; c < local_cells.n; c++) {
      Particle *part;
      int npart;
      cell = local_cells.cell[c];
      part = cell->part;
      npart = cell->n;
      memcpy(&result[g],part,npart*sizeof(Particle));
      g+=cell->n;
      if (bi) {
	int pc;
	for (pc = 0; pc < npart; pc++) {
	  Particle *p = &part[pc];
	  realloc_intlist(&local_bi, local_bi.n + p->bl.n);
	  memcpy(&local_bi.e[local_bi.n], p->bl.e, p->bl.n*sizeof(int));
	  local_bi.n += p->bl.n;
	}
      }
    }
    /* and send it back to the master node */
    MPI_Send(result, n_part*sizeof(Particle), MPI_BYTE, 0, REQ_GETPARTS, MPI_COMM_WORLD);
    free(result);

    if (bi) {
      COMM_TRACE(fprintf(stderr, "%d: sending bonds\n", this_node));

      MPI_Gather(&local_bi.n, 1, MPI_INT, NULL, 1, MPI_INT, 0, MPI_COMM_WORLD);      
      if (local_bi.n > 0)
	MPI_Send(local_bi.e, local_bi.n, MPI_INT, 0, REQ_GETPARTS, MPI_COMM_WORLD);
      realloc_intlist(&local_bi, 0);
    }
  }
  else {
    if (bi) {
      /* inform master node that we do not have bonds (as we don't have particles) */
      g = 0;
      MPI_Gather(&g, 1, MPI_INT, NULL, 1, MPI_INT, 0, MPI_COMM_WORLD);      
    }
  }
}

/*************** REQ_SET_TIME_STEP ************/
void mpi_set_time_step(double time_s)
{
  double old_ts = time_step;

  mpi_issue(REQ_SET_TIME_STEP, -1, 0);

  time_step = time_s;
  MPI_Bcast(&time_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  rescale_velocities(time_step / old_ts);
  on_parameter_change(FIELD_TIMESTEP);
}

void mpi_set_time_step_slave(int node, int i)
{
  double old_ts = time_step;

  MPI_Bcast(&time_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  rescale_velocities(time_step / old_ts);
  on_parameter_change(FIELD_TIMESTEP);
}

/*************** REQ_BCAST_COULOMB ************/
void mpi_bcast_coulomb_params()
{
#ifdef ELECTROSTATICS
  mpi_issue(REQ_BCAST_COULOMB, 1, 0);
  mpi_bcast_coulomb_params_slave(-1, 0);
#endif
}

void mpi_bcast_coulomb_params_slave(int node, int parm)
{   
#ifdef ELECTROSTATICS
  MPI_Bcast(&coulomb, sizeof(Coulomb_parameters), MPI_BYTE, 0, MPI_COMM_WORLD);
  switch (coulomb.method) {
  case COULOMB_NONE:
    break;
  case COULOMB_P3M:
    MPI_Bcast(&p3m, sizeof(p3m_struct), MPI_BYTE, 0, MPI_COMM_WORLD);
    break;
  case COULOMB_DH:
  case COULOMB_DH_PW:
    MPI_Bcast(&dh_params, sizeof(Debye_hueckel_params), MPI_BYTE, 0, MPI_COMM_WORLD);
    break;
  case COULOMB_MMM1D:
    MPI_Bcast(&mmm1d_params, sizeof(MMM1D_struct), MPI_BYTE, 0, MPI_COMM_WORLD);
    break;
  case COULOMB_MMM2D:
    MPI_Bcast(&mmm2d_params, sizeof(MMM2D_struct), MPI_BYTE, 0, MPI_COMM_WORLD);
    break;  
  case COULOMB_MAGGS:
    MPI_Bcast(&maggs, sizeof(MAGGS_struct), MPI_BYTE, 0, MPI_COMM_WORLD); 
    break;
  default:
    fprintf(stderr, "%d: INTERNAL ERROR: cannot bcast coulomb params for unknown method %d\n", this_node, coulomb.method);
    errexit();
  }
  if (coulomb.use_elc) {
    MPI_Bcast(&elc_params, sizeof(ELC_struct), MPI_BYTE, 0, MPI_COMM_WORLD);
  }
  on_coulomb_change();
  on_short_range_ia_change();
#endif
}

/****************** REQ_SET_EXT ************/
void mpi_send_ext(int pnode, int part, int flag, int mask, double force[3])
{
#ifdef EXTERNAL_FORCES
  int s_buf[2];
  mpi_issue(REQ_SET_EXT, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    /* mask out old flags */
    p->l.ext_flag &= ~mask;
    /* set new values */
    p->l.ext_flag |= flag;
    if (mask & PARTICLE_EXT_FORCE)
      memcpy(p->l.ext_force, force, 3*sizeof(double));
  }
  else {
    s_buf[0] = flag; s_buf[1] = mask;
    MPI_Send(s_buf, 2, MPI_INT, pnode, REQ_SET_EXT, MPI_COMM_WORLD);
    if (mask & PARTICLE_EXT_FORCE)
      MPI_Send(force, 3, MPI_DOUBLE, pnode, REQ_SET_EXT, MPI_COMM_WORLD);
  }

  on_particle_change();
#endif
}

void mpi_send_ext_slave(int pnode, int part)
{
#ifdef EXTERNAL_FORCES
  int s_buf[2];
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Status status;
    MPI_Recv(s_buf, 2, MPI_INT, 0, REQ_SET_EXT, MPI_COMM_WORLD, &status);
    /* mask out old flags */
    p->l.ext_flag &= ~s_buf[1];
    /* set new values */
    p->l.ext_flag |= s_buf[0];
    
    if (s_buf[1] & PARTICLE_EXT_FORCE)
      MPI_Recv(p->l.ext_force, 3, MPI_DOUBLE, 0, REQ_SET_EXT, MPI_COMM_WORLD, &status);
  }

  on_particle_change();
#endif
}

/*************** REQ_BCAST_CONSTR ************/
void mpi_bcast_constraint(int del_num)
{
#ifdef CONSTRAINTS
  mpi_issue(REQ_BCAST_CONSTR, 0, del_num);

  if (del_num == -1) {
    /* bcast new constraint */
    MPI_Bcast(&constraints[n_constraints-1], sizeof(Constraint), MPI_BYTE, 0, MPI_COMM_WORLD);
  }
  else if (del_num == -2) {
    /* delete all constraints */
    n_constraints = 0;
    constraints = realloc(constraints,n_constraints*sizeof(Constraint));
  }
  else {
    memcpy(&constraints[del_num],&constraints[n_constraints-1],sizeof(Constraint));
    n_constraints--;
    constraints = realloc(constraints,n_constraints*sizeof(Constraint));
  }

  on_constraint_change();
#endif
}

void mpi_bcast_constraint_slave(int node, int parm)
{   
#ifdef CONSTRAINTS
  if(parm == -1) {
    n_constraints++;
    constraints = realloc(constraints,n_constraints*sizeof(Constraint));
    MPI_Bcast(&constraints[n_constraints-1], sizeof(Constraint), MPI_BYTE, 0, MPI_COMM_WORLD);
  }
  else if (parm == -2) {
    /* delete all constraints */
    n_constraints = 0;
    constraints = realloc(constraints,n_constraints*sizeof(Constraint));
  }
  else {
    memcpy(&constraints[parm],&constraints[n_constraints-1],sizeof(Constraint));
    n_constraints--;
    constraints = realloc(constraints,n_constraints*sizeof(Constraint));    
  }

  on_constraint_change();
#endif
}

/*************** REQ_RANDOM_SEED ************/
void mpi_random_seed(int cnt, long *seed) {
  long this_idum = print_random_seed();

  mpi_issue(REQ_RANDOM_SEED, -1, cnt);

  if (cnt==0) {
    RANDOM_TRACE(printf("%d: Have seed %ld\n",this_node,this_idum));
    MPI_Gather(&this_idum,1,MPI_LONG,seed,1,MPI_LONG,0,MPI_COMM_WORLD); }
  else {
    MPI_Scatter(seed,1,MPI_LONG,&this_idum,1,MPI_LONG,0,MPI_COMM_WORLD);
    RANDOM_TRACE(printf("%d: Received seed %ld\n",this_node,this_idum));
    init_random_seed(this_idum);
  }
}

void mpi_random_seed_slave(int pnode, int cnt) {
  long this_idum = print_random_seed();

  if (cnt==0) {
    RANDOM_TRACE(printf("%d: Have seed %ld\n",this_node,this_idum));
    MPI_Gather(&this_idum,1,MPI_LONG,NULL,0,MPI_LONG,0,MPI_COMM_WORLD); }
  else {
    MPI_Scatter(NULL,1,MPI_LONG,&this_idum,1,MPI_LONG,0,MPI_COMM_WORLD);
    RANDOM_TRACE(printf("%d: Received seed %ld\n",this_node,this_idum));
    init_random_seed(this_idum);
  }
}

/*************** REQ_RANDOM_STAT ************/
void mpi_random_stat(int cnt, RandomStatus *stat) {
  RandomStatus this_stat = print_random_stat();

  mpi_issue(REQ_RANDOM_STAT, -1, cnt);

  if (cnt==0) {
    RANDOM_TRACE(printf("%d: Have status %ld/%ld/...\n",this_node,this_stat.idum,this_stat.iy));
    MPI_Gather(&this_stat,1*sizeof(RandomStatus),MPI_BYTE,stat,1*sizeof(RandomStatus),MPI_BYTE,0,MPI_COMM_WORLD); }
  else {
    MPI_Scatter(stat,1*sizeof(RandomStatus),MPI_BYTE,&this_stat,1*sizeof(RandomStatus),MPI_BYTE,0,MPI_COMM_WORLD);
    RANDOM_TRACE(printf("%d: Received status %ld/%ld/...\n",this_node,this_stat.idum,this_stat.iy));
    init_random_stat(this_stat);
  }
}

void mpi_random_stat_slave(int pnode, int cnt) {
  RandomStatus this_stat = print_random_stat();

  if (cnt==0) {
    RANDOM_TRACE(printf("%d: Have status %ld/%ld/...\n",this_node,this_stat.idum,this_stat.iy));
    MPI_Gather(&this_stat,1*sizeof(RandomStatus),MPI_BYTE,NULL,0,MPI_BYTE,0,MPI_COMM_WORLD); }
  else {
    MPI_Scatter(NULL,0,MPI_BYTE,&this_stat,1*sizeof(RandomStatus),MPI_BYTE,0,MPI_COMM_WORLD);
    RANDOM_TRACE(printf("%d: Received status %ld/%ld/...\n",this_node,this_stat.idum,this_stat.iy));
    init_random_stat(this_stat);
  }
}

/*************** REQ_BCAST_LJFORCECAP ************/
void mpi_lj_cap_forces(double fc)
{
  lj_force_cap = fc;
  mpi_issue(REQ_BCAST_LFC, 1, 0);
  mpi_lj_cap_forces_slave(1, 0);
}

void mpi_lj_cap_forces_slave(int node, int parm)
{
#ifdef LENNARD_JONES
  MPI_Bcast(&lj_force_cap, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  calc_lj_cap_radii(lj_force_cap);
  on_short_range_ia_change();
#endif
}

/*************** REQ_BCAST_MORSEFORCECAP ************/
void mpi_morse_cap_forces(double fc)
{
  morse_force_cap = fc;
  mpi_issue(REQ_BCAST_MFC, 1, 0);
  mpi_morse_cap_forces_slave(1, 0);
}

void mpi_morse_cap_forces_slave(int node, int parm)
{
#ifdef MORSE
  MPI_Bcast(&morse_force_cap, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  calc_morse_cap_radii(morse_force_cap);
  on_short_range_ia_change();
#endif
}

/*************** REQ_BCAST_BUCKFORCECAP ************/
void mpi_buck_cap_forces(double fc)
{
#ifdef BUCKINGHAM
  buck_force_cap = fc;
  mpi_issue(REQ_BCAST_BFC, 1, 0);
  mpi_buck_cap_forces_slave(1, 0);
#endif
}

void mpi_buck_cap_forces_slave(int node, int parm)
{
#ifdef BUCKINGHAM
  MPI_Bcast(&buck_force_cap, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  calc_buck_cap_radii(buck_force_cap);
  on_short_range_ia_change();
#endif
}

/*************** REQ_BCAST_TABFORCECAP ************/
void mpi_tab_cap_forces(double fc)
{
#ifdef TABULATED
  tab_force_cap = fc;
  mpi_issue(REQ_BCAST_TFC, 1, 0);
  mpi_tab_cap_forces_slave(1, 0);
#endif
}

void mpi_tab_cap_forces_slave(int node, int parm)
{
#ifdef TABULATED
  MPI_Bcast(&tab_force_cap, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  check_tab_forcecap(tab_force_cap);
  on_short_range_ia_change();
#endif
}

/*************** REQ_GET_CONSFOR ************/
void mpi_get_constraint_force(int cons, double force[3])
{
#ifdef CONSTRAINTS
  mpi_issue(REQ_GET_CONSFOR, -1, cons);
  MPI_Reduce(constraints[cons].part_rep.f.f, force, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
}

void mpi_get_constraint_force_slave(int node, int parm)
{
#ifdef CONSTRAINTS
  MPI_Reduce(constraints[parm].part_rep.f.f, NULL, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
}

/*************** REQ_BIT_RANDOM_SEED ************/
void mpi_bit_random_seed(int cnt, int *seed) {
  int this_idum = print_bit_random_seed();

  mpi_issue(REQ_BIT_RANDOM_SEED, -1, cnt);

  if (cnt==0) {
    RANDOM_TRACE(printf("%d: Have seed %d\n",this_node,this_idum));
    MPI_Gather(&this_idum,1,MPI_INT,seed,1,MPI_INT,0,MPI_COMM_WORLD); }
  else {
    MPI_Scatter(seed,1,MPI_INT,&this_idum,1,MPI_INT,0,MPI_COMM_WORLD);
    RANDOM_TRACE(printf("%d: Received seed %d\n",this_node,this_idum));
    init_bit_random_generator(this_idum);
  }
}

void mpi_bit_random_seed_slave(int pnode, int cnt) {
  int this_idum = print_bit_random_seed();

  if (cnt==0) {
    RANDOM_TRACE(printf("%d: Have seed %d\n",this_node,this_idum));
    MPI_Gather(&this_idum,1,MPI_INT,NULL,0,MPI_INT,0,MPI_COMM_WORLD); }
  else {
    MPI_Scatter(NULL,1,MPI_INT,&this_idum,1,MPI_INT,0,MPI_COMM_WORLD);
    RANDOM_TRACE(printf("%d: Received seed %d\n",this_node,this_idum));
    init_bit_random_generator(this_idum);
  }
}

/*************** REQ_BIT_RANDOM_STAT ************/
void mpi_bit_random_stat(int cnt, BitRandomStatus *stat) {
  BitRandomStatus this_stat = print_bit_random_stat();

  mpi_issue(REQ_BIT_RANDOM_STAT, -1, cnt);

  if (cnt==0) {
    RANDOM_TRACE(printf("%d: Have status %d/%d/...\n",this_node,this_stat.random_pointer_1,this_stat.random_pointer_2));
    MPI_Gather(&this_stat,1*sizeof(BitRandomStatus),MPI_BYTE,stat,1*sizeof(BitRandomStatus),MPI_BYTE,0,MPI_COMM_WORLD); }
  else {
    MPI_Scatter(stat,1*sizeof(BitRandomStatus),MPI_BYTE,&this_stat,1*sizeof(BitRandomStatus),MPI_BYTE,0,MPI_COMM_WORLD);
    RANDOM_TRACE(printf("%d: Received status %d/%d/...\n",this_node,this_stat.random_pointer_1,this_stat.random_pointer_2));
    init_bit_random_stat(this_stat);
  }
}

void mpi_bit_random_stat_slave(int pnode, int cnt) {
  BitRandomStatus this_stat = print_bit_random_stat();

  if (cnt==0) {
    RANDOM_TRACE(printf("%d: Have status %d/%d/...\n",this_node,this_stat.random_pointer_1,this_stat.random_pointer_2));
    MPI_Gather(&this_stat,1*sizeof(BitRandomStatus),MPI_BYTE,NULL,0,MPI_BYTE,0,MPI_COMM_WORLD); }
  else {
    MPI_Scatter(NULL,0,MPI_BYTE,&this_stat,1*sizeof(BitRandomStatus),MPI_BYTE,0,MPI_COMM_WORLD);
    RANDOM_TRACE(printf("%d: Received status %d/%d/...\n",this_node,this_stat.random_pointer_1,this_stat.random_pointer_2));
    init_bit_random_stat(this_stat);
  }
}

/****************** REQ_RESCALE_PART ************/

void mpi_rescale_particles(int dir, double scale) {
  int pnode;

  mpi_issue(REQ_RESCALE_PART, -1, dir);
  for (pnode = 0; pnode < n_nodes; pnode++) {
    if (pnode == this_node) {
      local_rescale_particles(dir, scale); }
    else {
      MPI_Send(&scale, 1, MPI_DOUBLE, pnode, REQ_PLACE, MPI_COMM_WORLD); }
  }
  on_particle_change();
}

void mpi_rescale_particles_slave(int pnode, int dir) {
  double scale; MPI_Status status;

  MPI_Recv(&scale, 1, MPI_DOUBLE, 0, REQ_PLACE, MPI_COMM_WORLD, &status);
  local_rescale_particles(dir, scale);
  on_particle_change();
}

/*************** REQ_BCAST_CS *****************/

void mpi_bcast_cell_structure(int cs)
{
  mpi_issue(REQ_BCAST_CS, -1, cs);
  mpi_bcast_cell_structure_slave(-1, cs);
}

void mpi_bcast_cell_structure_slave(int pnode, int cs)
{
  cells_re_init(cs);
  on_cell_structure_change();
}

/*************** REQ_BCAST_NPTISO_GEOM *****************/

void mpi_bcast_nptiso_geom()
{
  mpi_issue(REQ_BCAST_NPTISO_GEOM, -1 , 0);
  mpi_bcast_nptiso_geom_slave(-1,0);

}

void mpi_bcast_nptiso_geom_slave(int node,int parm)
{
  MPI_Bcast(&nptiso.geometry, 1,MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nptiso.dimension, 1,MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nptiso.cubic_box, 1,MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nptiso.non_const_dim, 1,MPI_INT, 0, MPI_COMM_WORLD);

}

/***************REQ_UPDATE_MOL_IDS *********************/

void mpi_update_mol_ids()
{
  mpi_issue(REQ_UPDATE_MOL_IDS, -1, 0);
  mpi_update_mol_ids_slave(-1, 0);
}

void mpi_update_mol_ids_slave(int node,int parm)
{
  update_mol_ids_setchains();
}

/******************* REQ_SYNC_TOPO ********************/
int mpi_sync_topo_part_info() {
  int i;
  int molsize;
  int moltype;
  int n_mols;
  mpi_issue(REQ_SYNC_TOPO,-1,0);
  n_mols = n_molecules;
  MPI_Bcast(&n_mols,1,MPI_INT,0,MPI_COMM_WORLD);

  for ( i = 0 ; i < n_molecules ; i++) {
    molsize = topology[i].part.n;
    moltype = topology[i].type;
    MPI_Bcast(&molsize,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&moltype,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(topology[i].part.e,topology[i].part.n,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&topology[i].type,1,MPI_INT,0,MPI_COMM_WORLD);
    
  }
  
  sync_topo_part_info();
  return 1;
}

void mpi_sync_topo_part_info_slave(int node,int parm ) {
  int i;
  int molsize;
  int moltype;
  int n_mols;

  MPI_Bcast(&n_mols,1,MPI_INT,0,MPI_COMM_WORLD);
  realloc_topology(n_mols);
  for ( i = 0 ; i < n_molecules ; i++) {

    MPI_Bcast(&molsize,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&moltype,1,MPI_INT,0,MPI_COMM_WORLD);
    topology[i].type = moltype;
    realloc_intlist(&topology[i].part,topology[i].part.n = molsize);

    MPI_Bcast(topology[i].part.e,topology[i].part.n,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&topology[i].type,1,MPI_INT,0,MPI_COMM_WORLD);

  }
  

  sync_topo_part_info();

}

int lb_set_params(double temp, double lbfriction, double viscosity, double tgrid,
		  double density, double agrid)
{
#ifdef LB
  temperature      = temp;
  lbpar.lbfriction = lbfriction;
  lbpar.viscos     = viscosity;
  lbpar.agrid      = agrid;
  lbpar.tau        = tgrid;
  lbpar.rho        = density;  
  lbpar.gridpoints = box_l[2]/agrid;  

  mpi_issue(REQ_LB_BCAST, 1, 0);
  mpi_bcast_lb_params_slave(-1, 0);
#endif  
  return TCL_OK;
}


void mpi_bcast_lb_params_slave(int node, int parm)
{
  /* broadcast lb parameters */
#ifdef LB  
  MPI_Bcast(&lbpar, sizeof(LB_structure), MPI_BYTE, 0, MPI_COMM_WORLD);
  thermo_init_lb();
#endif  
}

/******************* REQ_GET_ERRS ********************/

int mpi_gather_runtime_errors(Tcl_Interp *interp, int error_code)
{
  char nr_buf[TCL_INTEGER_SPACE + 3];
  char *other_error_msg;
  int *errcnt;
  int node, n_other_error_msg;
  MPI_Status status;

  mpi_issue(REQ_GET_ERRS, -1, 0);

  if (!check_runtime_errors())
    return error_code;

  if (error_code != TCL_ERROR)
    Tcl_ResetResult(interp);
  else
    Tcl_AppendResult(interp, " ", (char *) NULL);

  Tcl_AppendResult(interp, "background_errors ", (char *) NULL);

  /* gather the maximum length of the error messages, and allocate transfer space */
  errcnt = malloc(n_nodes*sizeof(int));
  MPI_Gather(&n_error_msg, 1, MPI_INT, errcnt, 1, MPI_INT, 0, MPI_COMM_WORLD);

  /* allocate transfer buffer for maximal error message length */
  n_other_error_msg = n_error_msg;
  for (node = 1; node < n_nodes; node++)
    if (errcnt[node] > n_other_error_msg)
      n_other_error_msg = errcnt[node];
  other_error_msg = malloc(n_other_error_msg);

  /* first handle node master errors. */
  if (n_error_msg > 0)
    Tcl_AppendResult(interp, "0 ", error_msg, (char *) NULL);
    
  for (node = 1; node < n_nodes; node++) {
    if (errcnt[node] > 0) {
      MPI_Recv(other_error_msg, errcnt[node], MPI_CHAR, node, 0, MPI_COMM_WORLD, &status);

      sprintf(nr_buf, "%d ", node);

      /* check wether it's the same message as from the master, then just consent */
      if (strcmp(other_error_msg, error_msg) == 0) {
	Tcl_AppendResult(interp, nr_buf, "<consent> ", (char *) NULL);
      }
      else {
	Tcl_AppendResult(interp, nr_buf, other_error_msg, (char *) NULL);
      }
    }
  }

  /* reset error message on master node */
  error_msg = realloc(error_msg, n_error_msg = 0);
  free(other_error_msg);
  free(errcnt);

  return TCL_ERROR;
}

void mpi_gather_runtime_errors_slave(int node, int parm)
{
  if (!check_runtime_errors())
    return;

  MPI_Gather(&n_error_msg, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);

  if (n_error_msg > 0) {
    MPI_Send(error_msg, n_error_msg, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    /* reset error message on slave node */
    error_msg = realloc(error_msg, n_error_msg = 0);
  }
}

/********************* REQ_SET_EXCL ********/
void mpi_send_exclusion(int part1, int part2, int delete)
{
#ifdef EXCLUSIONS
  mpi_issue(REQ_SET_EXCL, part1, part2);

  MPI_Bcast(&delete, 1, MPI_INT, 0, MPI_COMM_WORLD);
  local_change_exclusion(part1, part2, delete);
  on_particle_change();
#endif
}

void mpi_send_exclusion_slave(int part1, int part2)
{
#ifdef EXCLUSIONS
  int delete;
  MPI_Bcast(&delete, 1, MPI_INT, 0, MPI_COMM_WORLD);  
  local_change_exclusion(part1, part2, delete);
  on_particle_change();
#endif
}

/*********************** MAIN LOOP for slaves ****************/

void mpi_loop()
{
  for (;;) {
#ifdef ASYNC_BARRIER
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    MPI_Bcast(request, 3, MPI_INT, 0, MPI_COMM_WORLD);
    COMM_TRACE(fprintf(stderr, "%d: processing %s %d...\n", this_node,
		       names[request[0]], request[1]));

    if ((request[0] < 0) || (request[0] >= REQ_MAXIMUM)) {
      fprintf(stderr, "%d: INTERNAL ERROR: unknown request %d\n", this_node, request[0]);
      errexit();
    }
    slave_callbacks[request[0]](request[1], request[2]);
    COMM_TRACE(fprintf(stderr, "%d: finished %s %d %d\n", this_node,
		       names[request[0]], request[1], request[2]));
  }
}

