/** \file communication.c
    Implementation of \ref communication.h "communication.h".
*/
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"
#include "communication.h"
#include "interaction_data.h"
#include "integrate.h"
#include "cells.h"
#include "global.h"
#include "grid.h"
#include "debug.h"
#include "utils.h"
#include "initialize.h"
#include "forces.h"
#include "p3m.h"
#include "statistics.h"
#include "random.h"
#include "lj.h"
#include "ljcos.h"

int this_node = -1;
int n_nodes = -1;

/**********************************************
 * slave callbacks.
 **********************************************/

/*************************************
 * requests in random access mode
 * KEEP ORDER IN SYNC WITH callbacks[]
 * ALSO ADD FOR EVERY REQUEST
 * THE CALLBACK IN callbacks[]
 *************************************/
/* slave callback procedure */
typedef void (SlaveCallback)(int node, int param);

/** Action number for \ref mpi_stop. */
#define REQ_TERM      0
/** Action number for \ref mpi_bcast_parameter. */
#define REQ_BCAST_PAR 1
/** Action number for \ref mpi_who_has. */
#define REQ_WHO_HAS   2
/** Action number for \ref mpi_bcast_event. */
#define REQ_EVENT    3
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
/** Total number of action numbers. */
#define REQ_MAXIMUM 25

#ifdef DIPOLAR_INTERACTION
/** Action number for \ref mpi_send_quat. */
#define REQ_SET_QUAT 25
/** Action number for \ref mpi_send_lambda. */
#define REQ_SET_LAMBDA 26
/** Action number for \ref mpi_send_torque. */
#define REQ_SET_TORQUE 27
/** Total number of action numbers. */
#undef REQ_MAXIMUM
#define REQ_MAXIMUM 28
#endif

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
#ifdef DIPOLAR_INTERACTION
void mpi_send_quat_slave(int node, int parm);
void mpi_send_lambda_slave(int node, int parm);
void mpi_send_torque_slave(int node, int parm);
#endif
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
/*@}*/

/** A list of which function has to be called for
    the issued command. */
SlaveCallback *callbacks[] = {
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
#ifdef DIPOLAR_INTERACTION
  mpi_send_quat_slave,              /* 25: REQ_SET_QUAT */
  mpi_send_lambda_slave,            /* 26: REQ_SET_LAMBDA */
  mpi_send_torque_slave,            /* 27: REQ_SET_TORQUE */
#endif  
};

/** Names to be printed when communication debugging is on. */
char *names[] = {
  "TERM"      , /*  0 */
  "BCAST_PAR" , /*  1 */
  "WHO_HAS"   , /*  2 */
  "EVENT"     , /*  3 */
  "SET_POS"   , /*  4 */
  "SET_V"     , /*  5 */
  "SET_F"     , /*  6 */
  "SET_Q"     , /*  7 */
  "SET_TYPE"  , /*  8 */
  "SET_BOND"  , /*  9 */
  "GET_PART"  , /* 10 */
  "INTEGRATE" , /* 11 */
  "BCAST_IA"  , /* 12 */
  "BCAST_IAS" , /* 13 */
  "GATHER"    , /* 14 */
  "TIME_STEP" , /* 15 */
  "GET_PARTS" , /* 16 */
  "BCAST_CIA" , /* 17 */
  "SEND_EXT"  , /* 18 */
  "PLACE_NEW" , /* 19 */
  "REM_PART"  , /* 20 */
  "BCAST_CON" , /* 21 */
  "RAND_SEED" , /* 22 */
  "RAND_STAT" , /* 23 */
  "BCAST_LFC" , /* 24 */
#ifdef DIPOLAR_INTERACTION  
  "SET_QUAT"  , /* 25 */
  "SET_LAMBDA", /* 26 */
  "SET_TORQUE", /* 27 */
#endif  
};

/** the requests are compiled here. So after a crash you get the last issued request */
static int request[3];

/**********************************************
 * procedures
 **********************************************/

#ifdef MPI_CORE
void mpi_core(MPI_Comm *comm, int *errcode,...) {
  fprintf(stderr, "Aborting due to MPI error %d, forcing core dump\n", *errcode);
  fflush(stderr);
  sleep(10);
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
void mpi_bcast_parameter(int i)
{
  mpi_issue(REQ_BCAST_PAR, -1, i);

  mpi_bcast_parameter_slave(-1, i);
}

void mpi_bcast_parameter_slave(int node, int i)
{
  switch (fields[i].type) {
  case TYPE_INT:
    MPI_Bcast((int *)fields[i].data, fields[i].dimension,
	      MPI_INT, 0, MPI_COMM_WORLD);
    break;
  case TYPE_DOUBLE:
    MPI_Bcast((double *)fields[i].data, fields[i].dimension,
	      MPI_DOUBLE, 0, MPI_COMM_WORLD);
    break;
  default: break;
  }
}

/*************** REQ_WHO_HAS ****************/

void mpi_who_has()
{
  Cell *c;
  int *sizes = malloc(sizeof(int)*n_nodes);
  int *pdata = NULL;
  int pdata_s = 0, i, m, n, o;
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
      INNER_CELLS_LOOP(m,n,o) {
	c = CELL_PTR(m,n,o);
	for (i = 0; i < c->pList.n; i++)
	  particle_node[c->pList.part[i].r.identity] = this_node;
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
  Cell *c;
  int npart, i, m, n, o;
  int *sendbuf;
  int n_part;

  n_part = cells_get_n_particles();
  MPI_Gather(&n_part, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
  if (n_part == 0)
    return;

  sendbuf = malloc(sizeof(int)*n_part);
  npart = 0;
  INNER_CELLS_LOOP(m,n,o) {
    c = CELL_PTR(m,n,o);
    for (i = 0; i < c->pList.n; i++)
      sendbuf[npart++] = c->pList.part[i].r.identity;
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
  case PARTICLE_CHANGED:
    on_particle_change();
    break;
  case INTERACTION_CHANGED:
    on_ia_change();
    break;
  case PARAMETER_CHANGED:
    on_parameter_change();  
    break;
  case TOPOLOGY_CHANGED:
    on_topology_change();  
    break;
  case P3M_COUNT_CHARGES:
    P3M_count_charged_particles();
    break;
  default:;
  }
}

/****************** REQ_PLACE/RE_PLACE_NEW ************/

void mpi_place_particle(int pnode, int part, int new, double p[3])
{
  if (new) {
    mpi_issue(REQ_PLACE_NEW, pnode, part);
    added_particle(part);
  }
  else
    mpi_issue(REQ_PLACE, pnode, part);

  if (pnode == this_node)
    local_place_particle(part, p);
  else
    MPI_Send(p, 3, MPI_DOUBLE, pnode, REQ_PLACE, MPI_COMM_WORLD);

  on_particle_change();
}

void mpi_place_particle_slave(int pnode, int part)
{
  double p[3];
  MPI_Status status;

  if (request[0] == REQ_PLACE_NEW)
    added_particle(part);

  if (pnode == this_node) {
    MPI_Recv(p, 3, MPI_DOUBLE, 0, REQ_PLACE, MPI_COMM_WORLD, &status);
    local_place_particle(part, p);
  }

  on_particle_change();
}

/****************** REQ_SET_V ************/
void mpi_send_v(int pnode, int part, double v[3])
{
  mpi_issue(REQ_SET_V, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    memcpy(p->v, v, 3*sizeof(double));
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
    MPI_Recv(p->v, 3, MPI_DOUBLE, 0, REQ_SET_V,
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
    memcpy(p->f, F, 3*sizeof(double));
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
    MPI_Recv(p->f, 3, MPI_DOUBLE, 0, REQ_SET_F,
	     MPI_COMM_WORLD, &status);
  }

  on_particle_change();
}

/********************* REQ_SET_Q ********/
void mpi_send_q(int pnode, int part, double q)
{
  mpi_issue(REQ_SET_Q, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->r.q = q;
  }
  else {
    MPI_Send(&q, 1, MPI_DOUBLE, pnode, REQ_SET_Q, MPI_COMM_WORLD);
  }

  on_particle_change();
}

void mpi_send_q_slave(int pnode, int part)
{
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Status status;
    MPI_Recv(&p->r.q, 1, MPI_DOUBLE, 0, REQ_SET_Q,
	     MPI_COMM_WORLD, &status);
  }

  on_particle_change();
}

/********************* REQ_SET_TYPE ********/
void mpi_send_type(int pnode, int part, int type)
{
  mpi_issue(REQ_SET_TYPE, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->r.type = type;
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
    MPI_Recv(&p->r.type, 1, MPI_INT, 0, REQ_SET_TYPE,
	     MPI_COMM_WORLD, &status);
  }

  on_particle_change();
}

#ifdef DIPOLAR_INTERACTION

/********************* REQ_SET_QUAT ********/

void mpi_send_quat(int pnode, int part, double quat[4])
{
  mpi_issue(REQ_SET_QUAT, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    memcpy(p->quat, quat, 4*sizeof(double));
  }
  else {
    MPI_Send(quat, 4, MPI_DOUBLE, pnode, REQ_SET_QUAT, MPI_COMM_WORLD);
  }

  on_particle_change();
}

void mpi_send_quat_slave(int pnode, int part)
{
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Status status;
    MPI_Recv(p->quat, 4, MPI_DOUBLE, 0, REQ_SET_QUAT,
	     MPI_COMM_WORLD, &status);
  }

  on_particle_change();
}


/********************* REQ_SET_LAMBDA ********/

void mpi_send_lambda(int pnode, int part, double lambda)
{
  mpi_issue(REQ_SET_LAMBDA, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    memcpy(&p->lambda, &lambda, 1*sizeof(double));
  }
  else {
    MPI_Send(&lambda, 1, MPI_DOUBLE, pnode, REQ_SET_LAMBDA, MPI_COMM_WORLD);
  }

  on_particle_change();
}

void mpi_send_lambda_slave(int pnode, int part)
{
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Status status;
    MPI_Recv(&p->lambda, 1, MPI_DOUBLE, 0, REQ_SET_LAMBDA,
	     MPI_COMM_WORLD, &status);
  }

  on_particle_change();
}

/********************* REQ_SET_TORQUE ********/

void mpi_send_torque(int pnode, int part, double torque[3])
{
  mpi_issue(REQ_SET_TORQUE, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    memcpy(p->torque, torque, 3*sizeof(double));
  }
  else {
    MPI_Send(torque, 3, MPI_DOUBLE, pnode, REQ_SET_TORQUE, MPI_COMM_WORLD);
  }

  on_particle_change();
}

void mpi_send_torque_slave(int pnode, int part)
{
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Status status;
    MPI_Recv(p->torque, 3, MPI_DOUBLE, 0, REQ_SET_TORQUE,
	     MPI_COMM_WORLD, &status);
  }

  on_particle_change();
}

#endif


/********************* REQ_SET_BOND ********/
int mpi_send_bond(int pnode, int part, int *bond, int delete)
{
  int bond_size, stat;
  MPI_Status status;

  mpi_issue(REQ_SET_BOND, pnode, part);

  bond_size = (bond) ? bonded_ia_params[bond[0]].num + 1 : 0;

  if (pnode == this_node)
    return local_change_bond(part, bond, delete);
  else {
    MPI_Send(&bond_size, 1, MPI_INT, pnode, REQ_SET_BOND, MPI_COMM_WORLD);
    if (bond_size)
      MPI_Send(bond, bond_size, MPI_INT, pnode, REQ_SET_BOND, MPI_COMM_WORLD);
    MPI_Send(&delete, 1, MPI_INT, pnode, REQ_SET_BOND, MPI_COMM_WORLD);
    MPI_Recv(&stat, 1, MPI_INT, pnode, REQ_SET_BOND, MPI_COMM_WORLD, &status);
    return stat;
  }

  on_particle_change();

  return 0;
}

void mpi_send_bond_slave(int pnode, int part)
{
  int bond_size, *bond, delete, stat;
  MPI_Status status;

  if (pnode != this_node)
    return;

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

/****************** REQ_GET_PART ************/
void mpi_recv_part(int pnode, int part, Particle *pdata)
{
  IntList *bl = &(pdata->bl);
  MPI_Status status;

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    memcpy(pdata, p, sizeof(Particle));
    bl->max = bl->n;
    if (bl->n > 0) {
      alloc_intlist(bl, bl->n);
      memcpy(bl->e, p->bl.e, sizeof(int)*bl->n);
    }
    else
      bl->e = NULL;
  }
  else {
    mpi_issue(REQ_GET_PART, pnode, part);

    pdata->r.identity = part;
    MPI_Recv(&pdata->r.type, 1, MPI_INT, pnode,
	     REQ_GET_PART, MPI_COMM_WORLD, &status);
    MPI_Recv(pdata->r.p, 3, MPI_DOUBLE, pnode,
	     REQ_GET_PART, MPI_COMM_WORLD, &status);
    MPI_Recv(&pdata->r.q, 1, MPI_DOUBLE, pnode,
	     REQ_GET_PART, MPI_COMM_WORLD, &status);
    MPI_Recv(pdata->f, 3, MPI_DOUBLE, pnode,
	     REQ_GET_PART, MPI_COMM_WORLD, &status);
#ifdef DIPOLAR_INTERACTION	      
    MPI_Recv(pdata->quat, 4, MPI_DOUBLE, pnode,
	     REQ_GET_PART, MPI_COMM_WORLD, &status);     
    MPI_Recv(&pdata->lambda, 1, MPI_DOUBLE, pnode,
	     REQ_GET_PART, MPI_COMM_WORLD, &status);     
    MPI_Recv(pdata->torque, 3, MPI_DOUBLE, pnode,
	     REQ_GET_PART, MPI_COMM_WORLD, &status);
#endif   
    MPI_Recv(pdata->i, 3, MPI_INT, pnode,
	     REQ_GET_PART, MPI_COMM_WORLD, &status);
    MPI_Recv(pdata->v, 3, MPI_DOUBLE, pnode,
	     REQ_GET_PART, MPI_COMM_WORLD, &status);
    MPI_Recv(&(bl->n), 1, MPI_INT, pnode,
	     REQ_GET_PART, MPI_COMM_WORLD, &status);
    bl->max = bl->n;
    if (bl->n > 0) {
      alloc_intlist(bl, bl->n);
      MPI_Recv(bl->e, bl->n, MPI_INT, pnode,
	       REQ_GET_PART, MPI_COMM_WORLD, &status);
    }
    else
      pdata->bl.e = NULL;
  }
}

void mpi_recv_part_slave(int pnode, int part)
{
  Particle *p;
  if (pnode != this_node)
    return;

  p = local_particles[part];

  MPI_Send(&p->r.type, 1, MPI_INT, 0, REQ_GET_PART,
	   MPI_COMM_WORLD);
  MPI_Send(p->r.p, 3, MPI_DOUBLE, 0, REQ_GET_PART,
	   MPI_COMM_WORLD);
  MPI_Send(&p->r.q, 1, MPI_DOUBLE, 0, REQ_GET_PART,
	   MPI_COMM_WORLD);
  MPI_Send(p->f, 3, MPI_DOUBLE, 0, REQ_GET_PART,
	   MPI_COMM_WORLD);
#ifdef DIPOLAR_INTERACTION
  MPI_Send(p->quat, 4, MPI_DOUBLE, 0, REQ_GET_PART,
	   MPI_COMM_WORLD);
  MPI_Send(&p->lambda, 1, MPI_DOUBLE, 0, REQ_GET_PART,
	   MPI_COMM_WORLD);	   
  MPI_Send(p->torque, 3, MPI_DOUBLE, 0, REQ_GET_PART,
	   MPI_COMM_WORLD);
#endif	 
  MPI_Send(p->i, 3, MPI_INT, 0, REQ_GET_PART,
	   MPI_COMM_WORLD);
  MPI_Send(p->v, 3, MPI_DOUBLE, 0, REQ_GET_PART,
	   MPI_COMM_WORLD);
  MPI_Send(&p->bl.n, 1, MPI_INT, 0, REQ_GET_PART,
	   MPI_COMM_WORLD);
  if (p->bl.n > 0)
    MPI_Send(p->bl.e, p->bl.n, MPI_INT, 0, REQ_GET_PART,
	     MPI_COMM_WORLD);
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
}

/********************* REQ_INTEGRATE ********/
void mpi_integrate(int n_steps)
{
  mpi_issue(REQ_INTEGRATE, -1, n_steps);

  mpi_integrate_slave(-1, n_steps);
}

void mpi_integrate_slave(int pnode, int task)
{
  integrate_vv(task);

  COMM_TRACE(fprintf(stderr, "%d: integration task %d done.\n",
		     this_node, task));
}

/*************** REQ_BCAST_IA ************/
void mpi_bcast_ia_params(int i, int j)
{
  mpi_issue(REQ_BCAST_IA, i, j);

  if (j>=0) {
    /* non-bonded interaction parameters */
    /* INCOMPATIBLE WHEN NODES USE DIFFERENT ARCHITECTURES */
    MPI_Bcast(get_ia_param(i, j), sizeof(IA_parameters), MPI_BYTE,
	      0, MPI_COMM_WORLD);
  }
  else {
    /* bonded interaction parameters */
    /* INCOMPATIBLE WHEN NODES USE DIFFERENT ARCHITECTURES */
    MPI_Bcast(&(bonded_ia_params[i]),sizeof(Bonded_ia_parameters), MPI_BYTE,
	      0, MPI_COMM_WORLD);
  }

  on_ia_change();
}

void mpi_bcast_ia_params_slave(int i, int j)
{
  if(j >= 0) { /* non-bonded interaction parameters */
    /* INCOMPATIBLE WHEN NODES USE DIFFERENT ARCHITECTURES */
    MPI_Bcast(get_ia_param(i, j), sizeof(IA_parameters), MPI_BYTE,
	      0, MPI_COMM_WORLD);
  }
  else { /* bonded interaction parameters */
    make_bond_type_exist(i); /* realloc bonded_ia_params on slave nodes! */
    MPI_Bcast(&(bonded_ia_params[i]),sizeof(Bonded_ia_parameters), MPI_BYTE,
	      0, MPI_COMM_WORLD);
  }

  on_ia_change();
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
  int task[2];

  mpi_issue(REQ_GATHER, -1, job);

  switch (job) {
  case 0:
    /* gather minimum distance */
    MPI_Gather(&minimum_part_dist, 1, MPI_DOUBLE, (double *)result,
	       1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    break;
  case 1:
    /* calculate and reduce (sum up) energies */
    task[0] = energy.init_status;
    task[1] = energy.ana_num;
    MPI_Bcast(task, 2, MPI_INT, 0, MPI_COMM_WORLD);
    calc_energy();
    break;
  case 2:
    /* calculate and reduce (sum up) virials */
    task[0] = virials.init_status;
    task[1] = virials.ana_num;
    MPI_Bcast(task, 2, MPI_INT, 0, MPI_COMM_WORLD);
    calc_virials();
    break;
  default:
    fprintf(stderr, "%d: illegal request %d for REQ_GATHER\n", this_node, job);
    errexit();
  }
}

void mpi_gather_stats_slave(int pnode, int job)
{
  int task[2];

  switch (job) {
  case 0:
    MPI_Gather(&minimum_part_dist, 1, MPI_DOUBLE, NULL,
               1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    break;
  case 1:
    /* calculate and reduce (sum up) energies */
    MPI_Bcast(task, 2, MPI_INT, 0, MPI_COMM_WORLD);
    if( task[0]==0 ) init_energies();
    energy.ana_num = task[1];
    calc_energy();
    break;
  case 2:
    /* calculate and reduce (sum up) virials */
    MPI_Bcast(task, 2, MPI_INT, 0, MPI_COMM_WORLD);
    if( task[0]==0 ) init_virials();
    virials.ana_num = task[1];
    calc_virials();
    break;
  default:
    fprintf(stderr, "%d: illegal request %d for REQ_GATHER\n", this_node, job);
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
  int m,n,o;
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
    fprintf(stderr,"%d: mpi_get_particles: n_total_particles %d, but I counted %d. Exiting...\n",
	    this_node, n_total_particles, tot_size);
    errexit();
  }

  /* fetch particle informations into 'result' */
  init_intlist(&local_bi);
  g = 0;
  for (pnode = 0; pnode < n_nodes; pnode++) {
    if (sizes[pnode] > 0) {
      if (pnode == this_node) {
	INNER_CELLS_LOOP(m, n, o) {
	  Particle *part = CELL_PTR(m,n,o)->pList.part;
	  int npart = CELL_PTR(m,n,o)->pList.n;
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
  COMM_TRACE(for(i = 0; i < tot_size; i++) {
	       printf("%d: %d %d %f (%f, %f, %f)\n", i, result[i].r.identity, result[i].r.type,
		      result[i].r.q, result[i].r.p[0], result[i].r.p[1], result[i].r.p[2]);
	     });

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
	printf("(%d) part %d: bonds ", i, result[i].r.identity);
	for(g = 0; g < result[i].bl.n; g++) printf("%d ", result[i].bl.e[g]);
	printf("\n");
      });
    }
    realloc_intlist(&local_bi, 0);
  }

  free(sizes); 
}

void mpi_get_particles_slave(int pnode, int bi)
{
  IntList local_bi;
  int n_part;
  int g;
  int m,n,o;
  Particle *result;

  n_part = cells_get_n_particles();
  /* first collect number of particles on each node */
  MPI_Gather(&n_part, 1, MPI_INT, NULL, 1, MPI_INT,
	     0, MPI_COMM_WORLD);

  if (n_part > 0) {
    /* get (unsorted) particle informations as an array of type 'particle' */
    /* then get the particle information */
    result = malloc(n_part*sizeof(Particle));

    init_intlist(&local_bi);
    
    g = 0;
    INNER_CELLS_LOOP(m, n, o) {
      Particle *part = CELL_PTR(m,n,o)->pList.part;
      int npart = CELL_PTR(m,n,o)->pList.n;
      memcpy(&result[g],part,npart*sizeof(Particle));
      g+=CELL_PTR(m,n,o)->pList.n;
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
      MPI_Gather(&local_bi.n, 1, MPI_INT, NULL, 1, MPI_INT, 0, MPI_COMM_WORLD);      
      if (local_bi.n > 0)
	MPI_Send(local_bi.e, local_bi.n, MPI_INT, 0, REQ_GETPARTS, MPI_COMM_WORLD);
      realloc_intlist(&local_bi, 0);
    }
  }
}

/*************** REQ_SET_TIME_STEP ************/
void mpi_set_time_step()
{

  mpi_issue(REQ_SET_TIME_STEP, -1, 0);
  MPI_Bcast(&time_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  rescale_velocities();
}

void mpi_set_time_step_slave(int node, int i)
{

  old_time_step = time_step;
  MPI_Bcast(&time_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  rescale_velocities();
}

/*************** REQ_BCAST_COULOMB ************/
void mpi_bcast_coulomb_params()
{
  mpi_issue(REQ_BCAST_COULOMB, 1, 0);

  MPI_Bcast(&coulomb, sizeof(Coulomb_parameters), MPI_BYTE, 0, MPI_COMM_WORLD);
  if(coulomb.method == COULOMB_P3M) {
    MPI_Bcast(&p3m, sizeof(p3m_struct), MPI_BYTE, 0, MPI_COMM_WORLD);
  }
  else if(coulomb.method == COULOMB_DH) {
    MPI_Bcast(&dh_params, sizeof(Debye_hueckel_params), MPI_BYTE, 0, MPI_COMM_WORLD);
  }
  on_ia_change();
}

void mpi_bcast_coulomb_params_slave(int node, int parm)
{   
  MPI_Bcast(&coulomb, sizeof(Coulomb_parameters), MPI_BYTE, 0, MPI_COMM_WORLD);
  if(coulomb.method == COULOMB_P3M) {
    MPI_Bcast(&p3m, sizeof(p3m_struct), MPI_BYTE, 0, MPI_COMM_WORLD);
  }
  else if(coulomb.method == COULOMB_DH) {
    MPI_Bcast(&dh_params, sizeof(Debye_hueckel_params), MPI_BYTE, 0, MPI_COMM_WORLD);
  }
  on_ia_change();
}

/****************** REQ_SET_EXT ************/
void mpi_send_ext(int pnode, int part, int flag, double force[3])
{
#ifdef EXTERNAL_FORCES
  mpi_issue(REQ_SET_EXT, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->ext_flag = flag;
    memcpy(p->ext_force, force, 3*sizeof(double));
  }
  else {
    MPI_Send(&flag, 1, MPI_INT, pnode, REQ_SET_EXT, MPI_COMM_WORLD);
    MPI_Send(force, 3, MPI_DOUBLE, pnode, REQ_SET_EXT, MPI_COMM_WORLD);
  }

  on_particle_change();
#endif
}

void mpi_send_ext_slave(int pnode, int part)
{
#ifdef EXTERNAL_FORCES
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Status status;
    MPI_Recv(&(p->ext_flag), 1, MPI_INT, 0, REQ_SET_EXT, MPI_COMM_WORLD, &status);
    MPI_Recv(p->ext_force, 3, MPI_DOUBLE, 0, REQ_SET_EXT, MPI_COMM_WORLD, &status);
  }

  on_particle_change();
#endif
}

/*************** REQ_BCAST_CONSTR ************/
void mpi_bcast_constraint(int del_num)
{
#ifdef CONSTRAINTS
  mpi_issue(REQ_BCAST_CONSTR, 0, del_num);

  if(del_num < 0) {
    MPI_Bcast(&constraints[n_constraints-1], sizeof(Constraint), MPI_BYTE, 0, MPI_COMM_WORLD);
  }
  else {
    memcpy(&constraints[del_num],&constraints[n_constraints-1],sizeof(Constraint));
    n_constraints--;
    constraints = realloc(constraints,n_constraints*sizeof(Constraint));
  }

  on_ia_change();
#endif
}

void mpi_bcast_constraint_slave(int node, int parm)
{   
#ifdef CONSTRAINTS
  if(parm < 0) {
    n_constraints++;
    constraints = realloc(constraints,n_constraints*sizeof(Constraint));
    MPI_Bcast(&constraints[n_constraints-1], sizeof(Constraint), MPI_BYTE, 0, MPI_COMM_WORLD);
  }
  else {
    memcpy(&constraints[parm],&constraints[n_constraints-1],sizeof(Constraint));
    n_constraints--;
    constraints = realloc(constraints,n_constraints*sizeof(Constraint));    
  }

  on_ia_change();
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
  MPI_Bcast(&lj_force_cap, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  calc_lj_cap_radii(lj_force_cap);
  on_ia_change();
}

/*********************** MAIN LOOP for slaves ****************/

void mpi_loop()
{
  for (;;) {
    MPI_Bcast(request, 3, MPI_INT, 0, MPI_COMM_WORLD);
    COMM_TRACE(fprintf(stderr, "%d: processing %s %d...\n", this_node,
		       names[request[0]], request[1]));
 
    if ((request[0] < 0) || (request[0] >= REQ_MAXIMUM)) {
      fprintf(stderr, "%d: FATAL internal error: unknown request %d\n", this_node, request[0]);
      errexit();
    }
    callbacks[request[0]](request[1], request[2]);
    COMM_TRACE(fprintf(stderr, "%d: finished %s %d %d\n", this_node,
		       names[request[0]], request[1], request[2]));
  }
}
