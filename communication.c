/** \file communication.c
    Implementation of \ref communication.h "communication.h".
*/
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "communication.h"
#include "interaction_data.h"
#include "integrate.h"
#include "cells.h"
#include "forces.h"
#include "global.h"
#include "debug.h"
#include "utils.h"
#include "grid.h"

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
typedef void (SlaveCallback)(int param);

/** Action number for \ref mpi_stop. */
#define REQ_TERM      0
/** Action number for \ref mpi_bcast_parameter. */
#define REQ_BCAST_PAR 1
/** Action number for \ref mpi_who_has. */
#define REQ_WHO_HAS   2
/** Action number for \ref mpi_attach_particle. */
#define REQ_CHTOPL    3
/** Action number for \ref mpi_send_pos. */
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
#define REQ_GATHER  14
/** Total number of action numbers. */
#define REQ_MAXIMUM   15

/** \name Slave Callbacks
    These functions are the slave node counterparts for the
    commands the master node issues. The master node function *
    corresponds to the slave node function *_slave.
*/
/*@{*/
void mpi_stop_slave(int parm);
void mpi_bcast_parameter_slave(int parm);
void mpi_who_has_slave(int parm);
void mpi_changed_topology_slave(int parm);
void mpi_place_particle_slave(int parm);
void mpi_send_v_slave(int parm);
void mpi_send_f_slave(int parm);
void mpi_send_q_slave(int parm);
void mpi_send_type_slave(int parm);
void mpi_send_bond_slave(int parm);
void mpi_recv_part_slave(int parm);
void mpi_integrate_slave(int parm);
void mpi_bcast_ia_params_slave(int parm);
void mpi_bcast_n_particle_types_slave(int parm);
void mpi_gather_stats_slave(int parm);
/*@}*/

/** A list of wich function has to be called for
    the issued command. */
SlaveCallback *callbacks[] = {
  mpi_stop_slave,                /*  0: REQ_TERM */
  mpi_bcast_parameter_slave,     /*  1: REQ_BCAST_PAR */
  mpi_who_has_slave,             /*  2: REQ_WHO_HAS */
  mpi_changed_topology_slave,    /*  3: REQ_CHTOPL */
  mpi_place_particle_slave,      /*  4: REQ_SET_POS */
  mpi_send_v_slave,              /*  5: REQ_SET_V */
  mpi_send_f_slave,              /*  6: REQ_SET_F */
  mpi_send_q_slave,              /*  7: REQ_SET_Q */
  mpi_send_type_slave,           /*  8: REQ_SET_TYPE */
  mpi_send_bond_slave,           /*  9: REQ_SET_BOND */
  mpi_recv_part_slave,           /* 10: REQ_GET_PART */
  mpi_integrate_slave,           /* 11: REQ_INTEGRATE */
  mpi_bcast_ia_params_slave,     /* 12: REQ_BCAST_IA */ 
  mpi_bcast_n_particle_types_slave, /* 13: REQ_BCAST_IA_SIZE */ 
  mpi_gather_stats_slave         /* 14: REQ_GATHER */ 
};

/** Names to be printed when communication debugging is on. */
char *names[] = {
  "TERM",       /*  0 */
  "BCAST_PAR" , /*  1 */
  "WHO_HAS"   , /*  2 */
  "CHANGETOP" , /*  3 */
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
};

/** the requests are compiled here. So after a crash you get the last issued request */
static int request[2];

/**********************************************
 * procedures
 **********************************************/

#ifdef MPI_CORE
void mpi_core(MPI_Comm *comm, int *errcode,...) {
  fprintf(stderr, "Aborting due to MPI error %d, forcing core dump\n", *errcode);
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

static void mpi_issue(int reqcode, int param)
{
  request[0] = reqcode;
  request[1] = param;

  COMM_TRACE(fprintf(stderr, "0: issuing %s\n", names[reqcode]));
  MPI_Bcast(request, 2, MPI_INT, 0, MPI_COMM_WORLD);
}

/**************** REQ_TERM ************/

static int terminated = 0;

void mpi_stop()
{
  if (terminated)
    return;

  mpi_issue(REQ_TERM, 0);

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  regular_exit = 1;
  terminated = 1;
}

void mpi_stop_slave(int param)
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
  mpi_issue(REQ_BCAST_PAR, i);

  mpi_bcast_parameter_slave(i);
}

void mpi_bcast_parameter_slave(int i)
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

  mpi_issue(REQ_WHO_HAS, 0);

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
	  particle_node[c->pList.part[i].r.identity] = pnode;
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

void mpi_who_has_slave(int ident)
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
void mpi_changed_topology()
{
  mpi_issue(REQ_CHTOPL, 0);
  mpi_changed_topology_slave(0);
}

void mpi_changed_topology_slave(int dummy)
{
  grid_changed_topology();  
}

/****************** REQ_SET_POS ************/

static void __place_particle(int part, double p[3])
{
  int i[3];
  Particle *pt = cells_got_particle(part);

  i[0] = 0;
  i[1] = 0;
  i[2] = 0;
  fold_particle(p, i);

  if (!pt)
    pt = cells_alloc_particle(part, p);
  memcpy(pt->r.p, p, 3*sizeof(double));
  memcpy(pt->i, i, 3*sizeof(int));
}

void mpi_place_particle(int pnode, int part, double p[3])
{
  if (pnode == this_node)
    __place_particle(part, p);
  else {
    mpi_issue(REQ_PLACE, pnode);
    MPI_Send(&part, 1, MPI_INT, pnode, REQ_PLACE, MPI_COMM_WORLD);
    MPI_Send(p, 3, MPI_DOUBLE, pnode, REQ_PLACE, MPI_COMM_WORLD);
  }
}

void mpi_place_particle_slave(int pnode)
{
  int part;
  double p[3];
  MPI_Status status;

  if (pnode != this_node)
    return;

  MPI_Recv(&part, 1, MPI_INT, 0, REQ_PLACE, MPI_COMM_WORLD, &status);
  MPI_Recv(p, 3, MPI_DOUBLE, 0, REQ_PLACE, MPI_COMM_WORLD, &status);
  __place_particle(part, p);
}

/****************** REQ_SET_V ************/
void mpi_send_v(int pnode, int part, double v[3])
{
  COMM_TRACE(fprintf(stderr, "0: issuing SET_V %d %d\n", pnode, part));
  if (pnode == this_node) {
    Particle *p = cells_got_particle(part);
    memcpy(p->v, v, 3*sizeof(double));
  }
  else {
    mpi_issue(REQ_SET_V, part);
    MPI_Send(v, 3, MPI_DOUBLE, pnode, REQ_SET_V, MPI_COMM_WORLD);
  }
}

void mpi_send_v_slave(int part)
{
  Particle *p = cells_got_particle(part);
  MPI_Status status;
  if (p)
    MPI_Recv(p->v, 3, MPI_DOUBLE, 0, REQ_SET_V,
	     MPI_COMM_WORLD, &status);
}

/****************** REQ_SET_F ************/
void mpi_send_f(int pnode, int part, double F[3])
{
  COMM_TRACE(fprintf(stderr, "0: issuing SET_F %d %d\n", pnode, part));
  if (pnode == this_node) {
    Particle *p = cells_got_particle(part);
    memcpy(p->f, F, 3*sizeof(double));
  }
  else {
    mpi_issue(REQ_SET_F, part);
    MPI_Send(F, 3, MPI_DOUBLE, pnode, REQ_SET_F, MPI_COMM_WORLD);
  }
}

void mpi_send_f_slave(int part)
{
  Particle *p = cells_got_particle(part);
  MPI_Status status;
  if (p)
    MPI_Recv(p->f, 3, MPI_DOUBLE, 0, REQ_SET_F,
	     MPI_COMM_WORLD, &status);
}

/********************* REQ_SET_Q ********/
void mpi_send_q(int pnode, int part, double q)
{
  COMM_TRACE(fprintf(stderr, "0: issuing SET_Q %d %d\n", pnode, part));
  if (pnode == this_node) {
    Particle *p = cells_got_particle(part);
    p->r.q = q;
  }
  else {
    mpi_issue(REQ_SET_Q, part);
    MPI_Send(&q, 1, MPI_DOUBLE, pnode, REQ_SET_Q, MPI_COMM_WORLD);
  }
}

void mpi_send_q_slave(int part)
{
  Particle *p = cells_got_particle(part);
  MPI_Status status;
  if (p)
    MPI_Recv(&p->r.q, 1, MPI_DOUBLE, 0, REQ_SET_Q,
	     MPI_COMM_WORLD, &status);
}

/********************* REQ_SET_TYPE ********/
void mpi_send_type(int pnode, int part, int type)
{
  COMM_TRACE(fprintf(stderr, "0: issuing SET_TYPE %d %d\n", pnode, part));
  if (pnode == this_node) {
    Particle *p = cells_got_particle(part);
    p->r.type = type;
  }
  else {
    mpi_issue(REQ_SET_TYPE, part);
    MPI_Send(&type, 1, MPI_INT, pnode, REQ_SET_TYPE, MPI_COMM_WORLD);
  }
}

void mpi_send_type_slave(int part)
{
  Particle *p = cells_got_particle(part);
  MPI_Status status;
  if (p)
    MPI_Recv(&p->r.type, 1, MPI_INT, 0, REQ_SET_TYPE,
	     MPI_COMM_WORLD, &status);
}

/********************* REQ_SET_BOND ********/
int mpi_send_bond(int pnode, int part, int *bond, int delete)
{
  IntList *bl;
  int i,bond_size, stat;
  MPI_Status status;

  bond_size = bonded_ia_params[bond[0]].num + 1;
  if (pnode == this_node) {
    Particle *p = cells_got_particle(part);
    if (delete)
      return try_delete_bond(p, bond);
    bl = &(p->bl);
    realloc_intlist(bl, bl->n + bond_size);
    for(i = 0; i < bond_size; i++)
      bl->e[bl->n++] = bond[i];
    return 1;
  }
  else {
    mpi_issue(REQ_SET_BOND, part);
    MPI_Send(&delete, 1, MPI_INT, pnode, REQ_SET_BOND, MPI_COMM_WORLD);
    MPI_Send(&bond_size, 1, MPI_INT, pnode, REQ_SET_BOND, MPI_COMM_WORLD);
    MPI_Send(bond, bond_size, MPI_INT, pnode, REQ_SET_BOND, MPI_COMM_WORLD);
    MPI_Recv(&stat, 1, MPI_INT, pnode, REQ_SET_BOND, MPI_COMM_WORLD, &status);
    return stat;
  }
  return 0;
}

void mpi_send_bond_slave(int part)
{
  IntList *bl = NULL;
  Particle *p = cells_got_particle(part);
  int start, bond_size, delete, stat, *insert;
  MPI_Status status;

  if (!p)
    return;

  MPI_Recv(&delete, 1, MPI_INT, 0, REQ_SET_BOND,
	   MPI_COMM_WORLD, &status);
  MPI_Recv(&bond_size, 1, MPI_INT, 0, REQ_SET_BOND,
	   MPI_COMM_WORLD, &status);
  if (delete)
    insert = (int *)malloc(bond_size*sizeof(int));
  else {
    bl = &(p->bl);
    start = bl->n;
    bl->n += bond_size;
    realloc_intlist(bl, bl->n);
    insert = bl->e + start;
  }
  MPI_Recv(insert, bond_size, MPI_INT, 0, REQ_SET_BOND, MPI_COMM_WORLD, &status);
  if (delete) {
    stat = try_delete_bond(p, insert);
    free(insert);
  }
  else
    stat = 1;
  MPI_Send(&stat, 1, MPI_INT, 0, REQ_SET_BOND, MPI_COMM_WORLD);
}

/****************** REQ_GET_PART ************/
void mpi_recv_part(int pnode, int part, Particle *pdata)
{
  IntList *bl = &(pdata->bl);
  MPI_Status status;

  if (pnode == this_node) {
    Particle *p = cells_got_particle(part);
    
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
    mpi_issue(REQ_GET_PART, part);

    pdata->r.identity = part;
    MPI_Recv(&pdata->r.type, 1, MPI_INT, pnode,
	     REQ_GET_PART, MPI_COMM_WORLD, &status);
    MPI_Recv(pdata->r.p, 3, MPI_DOUBLE, pnode,
	     REQ_GET_PART, MPI_COMM_WORLD, &status);
    MPI_Recv(&pdata->r.q, 1, MPI_DOUBLE, pnode,
	     REQ_GET_PART, MPI_COMM_WORLD, &status);
    MPI_Recv(pdata->f, 3, MPI_DOUBLE, pnode,
	     REQ_GET_PART, MPI_COMM_WORLD, &status);
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

void mpi_recv_part_slave(int part)
{
  Particle *p;
  if (!(p = cells_got_particle(part)))
    return;
  
  MPI_Send(&p->r.type, 1, MPI_INT, 0, REQ_GET_PART,
	   MPI_COMM_WORLD);
  MPI_Send(p->r.p, 3, MPI_DOUBLE, 0, REQ_GET_PART,
	   MPI_COMM_WORLD);
  MPI_Send(&p->r.q, 1, MPI_DOUBLE, 0, REQ_GET_PART,
	   MPI_COMM_WORLD);
  MPI_Send(p->f, 3, MPI_DOUBLE, 0, REQ_GET_PART,
	   MPI_COMM_WORLD);
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

/********************* REQ_INTEGRATE ********/
void mpi_integrate(int n_steps)
{
  mpi_issue(REQ_INTEGRATE, n_steps);

  mpi_integrate_slave(n_steps);
}

void mpi_integrate_slave(int task)
{
  if(task == -1) {
    integrate_vv_init();
  }
  else if(task > 0)
    integrate_vv(task);
  else if(task < -1)
    integrate_vv_exit();
  
  COMM_TRACE(fprintf(stderr, "%d: integration task %d done.\n",
		     this_node, task));
}

/*************** REQ_BCAST_IA ************/
void mpi_bcast_ia_params(int i, int j)
{
  mpi_issue(REQ_BCAST_IA, i);

  MPI_Bcast(&j,  1, MPI_INT, 0, MPI_COMM_WORLD);

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
}

void mpi_bcast_ia_params_slave(int i)
{
  int j;
  MPI_Bcast(&j,  1, MPI_INT, 0, MPI_COMM_WORLD);

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
}

/*************** REQ_BCAST_IA_SIZE ************/

void mpi_bcast_n_particle_types(int ns)
{
  mpi_issue(REQ_BCAST_IA_SIZE, ns);
  mpi_bcast_n_particle_types_slave(ns);
}

void mpi_bcast_n_particle_types_slave(int ns)
{
  realloc_ia_params(ns);
}

/*************** REQ_GATHER ************/
void mpi_gather_stats(int job, void *result)
{
  int n_part;
  int tot_size, i, c, g, pnode;
  int *sizes;
  float *cdata;
  float_packed_particle_data *fdata;
  MPI_Status status;

  mpi_issue(REQ_GATHER, job);

  switch (job) {
  case 0:
    /* gather minimum distance */
    MPI_Gather(&minimum_part_dist, 1, MPI_DOUBLE, (double *)result,
	       1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    break;
  case 1:
    /* get packed particle coords - similar to REQ_WHO_HAS */
    sizes = malloc(sizeof(int)*n_nodes);
    fdata = (float_packed_particle_data *)result;
    n_part = cells_get_n_particles();

    /* first collect number of particles on each node */
    MPI_Gather(&n_part, 1, MPI_INT, sizes, 1, MPI_INT,
	       0, MPI_COMM_WORLD);
    tot_size = 0;
    for (i = 0; i < n_nodes; i++)
      tot_size += sizes[i];

    fdata->n_particles = tot_size;
    fdata->coords =
      cdata = malloc(tot_size*3*sizeof(float));

    /* then fetch particle locations */
    for (pnode = 0; pnode < n_nodes; pnode++) {
      COMM_TRACE(fprintf(stderr, "node %d reports %d particles\n",
			 pnode, sizes[pnode]));
      if (sizes[pnode] > 0) {
	if (pnode == this_node) {
	  g = 0;
	  for (c = 0; c < n_cells; c++)
	    for (i = 0; i < cells[c].pList.n; i++) {
	      cdata[g++] = cells[c].pList.part[i].r.p[0];
	      cdata[g++] = cells[c].pList.part[i].r.p[1];
	      cdata[g++] = cells[c].pList.part[i].r.p[2];
	    }
	}
	else {
	  MPI_Recv(cdata, 3*sizes[pnode], MPI_FLOAT, pnode, REQ_GATHER,
		   MPI_COMM_WORLD, &status);
	}
	cdata += 3*sizes[pnode];
      }
    }

    free(sizes);
    COMM_TRACE(cdata = fdata->coords;
	       for (i = 0; i < fdata->n_particles; i++) {
		 printf("%d %f %f %f\n", i, cdata[3*i],
			cdata[3*i+1], cdata[3*i+2]);
	       });
    break;
  default: ;
  }
}

void mpi_gather_stats_slave(int job)
{
  int n_part;
  int i, c, g;
  float *cdata;

  switch (job) {
  case 0:
    MPI_Gather(&minimum_part_dist, 1, MPI_DOUBLE, NULL,
	       1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    break;
  case 1:
    /* get packed particle coords - similar to REQ_WHO_HAS */
    /* first collect number of particles on each node */
    MPI_Gather(&n_part, 1, MPI_INT, NULL, 1, MPI_INT,
	       0, MPI_COMM_WORLD);
    cdata = malloc(3*n_part*sizeof(float));
    g = 0;
    for (c = 0; c < n_cells; c++)
      for (i = 0; i < cells[c].pList.n; i++) {
	cdata[g++] = cells[c].pList.part[i].r.p[0];
	cdata[g++] = cells[c].pList.part[i].r.p[1];
	cdata[g++] = cells[c].pList.part[i].r.p[2];
      }
    MPI_Send(cdata, 3*n_part, MPI_FLOAT, 0, REQ_GATHER,
	     MPI_COMM_WORLD);

    free(cdata);
    break;
  default:;
  }
}

/*********************** MAIN LOOP for slaves ****************/

void mpi_loop()
{
  for (;;) {
    MPI_Bcast(request, 2, MPI_INT, 0, MPI_COMM_WORLD);
    COMM_TRACE(fprintf(stderr, "%d: processing %s %d...\n", this_node,
		       names[request[0]], request[1]));
 
    if ((request[0] < 0) || (request[0] >= REQ_MAXIMUM))
      continue;
    callbacks[request[0]](request[1]);
    COMM_TRACE(fprintf(stderr, "%d: finished %s %d\n", this_node,
		       names[request[0]], request[1]));
  }
}
