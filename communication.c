#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "communication.h"
#include "interaction_data.h"
#include "integrate.h"
#include "forces.h"
#include "global.h"
#include "debug.h"
/** \file communication.c
    Implementation of \ref communication.h "communication.h".
*/

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
#define REQ_ATTACH    3
/** Action number for \ref mpi_send_pos. */
#define REQ_SET_POS   4
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
void mpi_attach_particle_slave(int parm);
void mpi_send_pos_slave(int parm);
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
  mpi_attach_particle_slave,     /*  3: REQ_ATTACH */
  mpi_send_pos_slave,            /*  4: REQ_SET_POS */
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
  "ATTACH"    , /*  3 */
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

/**************** REQ_TERM ************/

static int terminated = 0;

void mpi_stop()
{
  request[0] = REQ_TERM;
  request[1] = 0;

  if (terminated)
    return;

  COMM_TRACE(fprintf(stderr, "0: issuing TERM\n"));
  MPI_Bcast(request, 2, MPI_INT, 0, MPI_COMM_WORLD);

  COMM_TRACE(fprintf(stderr, "%d: sent TERM\n", this_node));

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
  request[0] = REQ_BCAST_PAR;
  request[1] = i;

  COMM_TRACE(fprintf(stderr, "0: issuing BCAST_PAR %d\n", i));
  MPI_Bcast(request, 2, MPI_INT, 0, MPI_COMM_WORLD);

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
  int *sizes = malloc(sizeof(int)*n_nodes);
  int *pdata = NULL;
  int pdata_s = 0, i;
  int pnode;
  MPI_Status status;

  request[0] = REQ_WHO_HAS;
  request[1] = 0;

  COMM_TRACE(fprintf(stderr, "0: issuing WHO_HAS\n"));

  MPI_Bcast(request, 2, MPI_INT, 0, MPI_COMM_WORLD);

  /* first collect number of particles on each node */
  MPI_Gather(&n_particles, 1, MPI_INT, sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);

  /* then fetch particle locations */
  for (pnode = 0; pnode < n_nodes; pnode++) {
    COMM_TRACE(fprintf(stderr, "node %d reports %d particles\n",
		       pnode, sizes[pnode]));
    if (pnode == this_node) {
      for (i = 0; i < n_particles; i++)
	particle_node[particles[i].identity] = pnode;
    }
    else if (sizes[pnode] > 0) {
      if (pdata_s < sizes[pnode]) {
	pdata_s = sizes[pnode];
	pdata = realloc(pdata, sizeof(int)*pdata_s);
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
  int npart, i;
  int *sendbuf;
  MPI_Gather(&n_particles, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
  if (n_particles == 0)
    return;

  sendbuf = malloc(sizeof(int)*n_particles);
  npart = 0;
  for (i = 0; i < n_particles; i++)
    sendbuf[npart++] = particles[i].identity;
  MPI_Send(sendbuf, npart, MPI_INT, 0, REQ_WHO_HAS, MPI_COMM_WORLD);
  free(sendbuf);
}

/**************** REQ_ATTACH ***********/
void mpi_attach_particle(int part, int pnode)
{
  request[0] = REQ_ATTACH;
  request[1] = pnode;

  COMM_TRACE(fprintf(stderr, "0: issuing ATTACH %d %d\n", pnode, part));

  MPI_Bcast(request, 2, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&part, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (pnode == this_node)
    add_particle(part);
}

void mpi_attach_particle_slave(int pnode)
{
  int part;

  MPI_Bcast(&part, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  if (pnode == this_node)
    add_particle(part);
}

/****************** REQ_SET_POS ************/
void mpi_send_pos(int pnode, int part, double p[3])
{
  COMM_TRACE(fprintf(stderr, "0: issuing SET_POS %d %d\n", pnode, part));
  if (pnode == this_node) {
    int index = got_particle(part);
    particles[index].p[0] = p[0];
    particles[index].p[1] = p[1];
    particles[index].p[2] = p[2];
  }
  else {
    request[0] = REQ_SET_POS;
    request[1] = part;
    MPI_Bcast(request, 2, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Send(p, 3, MPI_DOUBLE, pnode, REQ_SET_POS, MPI_COMM_WORLD);
  }
}

void mpi_send_pos_slave(int part)
{
  int index = got_particle(part);
  MPI_Status status;
  if (index != -1) {
    MPI_Recv(particles[index].p, 3, MPI_DOUBLE, 0, REQ_SET_POS,
	     MPI_COMM_WORLD, &status);
  }
}

/****************** REQ_SET_V ************/
void mpi_send_v(int pnode, int part, double v[3])
{
  COMM_TRACE(fprintf(stderr, "0: issuing SET_V %d %d\n", pnode, part));
  if (pnode == this_node) {
    int index = got_particle(part);
    particles[index].v[0] = v[0];
    particles[index].v[1] = v[1];
    particles[index].v[2] = v[2];
  }
  else {
    request[0] = REQ_SET_V;
    request[1] = part;
    MPI_Bcast(request, 2, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Send(v, 3, MPI_DOUBLE, pnode, REQ_SET_V, MPI_COMM_WORLD);
  }
}

void mpi_send_v_slave(int part)
{
  int index = got_particle(part);
  MPI_Status status;
  if (index != -1) {
    MPI_Recv(particles[index].v, 3, MPI_DOUBLE, 0, REQ_SET_V,
	     MPI_COMM_WORLD, &status);
  }
}

/****************** REQ_SET_F ************/
void mpi_send_f(int pnode, int part, double F[3])
{
  COMM_TRACE(fprintf(stderr, "0: issuing SET_F %d %d\n", pnode, part));
  if (pnode == this_node) {
    int index = got_particle(part);
    particles[index].f[0] = F[0];
    particles[index].f[1] = F[1];
    particles[index].f[2] = F[2];
  }
  else {
    request[0] = REQ_SET_F;
    request[1] = part;
    MPI_Bcast(request, 2, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Send(F, 3, MPI_DOUBLE, pnode, REQ_SET_F, MPI_COMM_WORLD);
  }
}

void mpi_send_f_slave(int part)
{
  int index = got_particle(part);
  MPI_Status status;
  if (index != -1) {
    MPI_Recv(particles[index].f, 3, MPI_DOUBLE, 0, REQ_SET_F,
	     MPI_COMM_WORLD, &status);
  }
}

/********************* REQ_SET_Q ********/
void mpi_send_q(int pnode, int part, double q)
{
  COMM_TRACE(fprintf(stderr, "0: issuing SET_Q %d %d\n", pnode, part));
  if (pnode == this_node) {
    int index = got_particle(part);
    particles[index].q = q;
  }
  else {
    request[0] = REQ_SET_Q;
    request[1] = part;
    MPI_Bcast(request, 2, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Send(&q, 1, MPI_DOUBLE, pnode, REQ_SET_Q, MPI_COMM_WORLD);
  }
}

void mpi_send_q_slave(int part)
{
  int index = got_particle(part);
  MPI_Status status;
  if (index != -1) {
    MPI_Recv(&particles[index].q, 1, MPI_DOUBLE, 0, REQ_SET_Q,
	     MPI_COMM_WORLD, &status);
  }
}

/********************* REQ_SET_TYPE ********/
void mpi_send_type(int pnode, int part, int type)
{
  COMM_TRACE(fprintf(stderr, "0: issuing SET_TYPE %d %d\n", pnode, part));
  if (pnode == this_node) {
    int index = got_particle(part);
    particles[index].type = type;
  }
  else {
    request[0] = REQ_SET_TYPE;
    request[1] = part;
    MPI_Bcast(request, 2, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Send(&type, 1, MPI_INT, pnode, REQ_SET_TYPE, MPI_COMM_WORLD);
  }
}

void mpi_send_type_slave(int part)
{
  int index = got_particle(part);
  MPI_Status status;
  if (index != -1) {
    MPI_Recv(&particles[index].type, 1, MPI_INT, 0, REQ_SET_TYPE,
	     MPI_COMM_WORLD, &status);
  }
}

/********************* REQ_SET_BOND ********/
void mpi_send_bond(int pnode, int part, int *bond)
{
  int i,bond_size;
  COMM_TRACE(fprintf(stderr, "0: issuing SET_BOND to node %d (part %d bond type_num %d)\n"
		     , pnode, part, bond[0]));
  bond_size = bonded_ia_params[bond[0]].num + 1;
  if (pnode == this_node) {
    int ind = got_particle(part);
    realloc_part_bonds(ind, particles[ind].n_bonds+bond_size);
    for(i=0;i<bond_size;i++)
      particles[ind].bonds[particles[ind].n_bonds+i] = bond[i];
    particles[ind].n_bonds += bond_size;
  }
  else {
    request[0] = REQ_SET_BOND;
    request[1] = part;
    MPI_Bcast(request, 2, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Send(&bond_size, 1, MPI_INT, pnode, REQ_SET_BOND, MPI_COMM_WORLD);
    MPI_Send(bond, bond_size, MPI_INT, pnode, REQ_SET_BOND, MPI_COMM_WORLD);
  }
}

void mpi_send_bond_slave(int part)
{
  int ind = got_particle(part);
  int bond_size;
  MPI_Status status;
  if (ind != -1) {
    MPI_Recv(&bond_size, 1, MPI_INT, 0, REQ_SET_BOND,
	     MPI_COMM_WORLD, &status);
    realloc_part_bonds(ind, particles[ind].n_bonds+bond_size);
    MPI_Recv(&(particles[ind].bonds[particles[ind].n_bonds]), bond_size, 
	     MPI_INT, 0, REQ_SET_BOND, MPI_COMM_WORLD, &status);
    particles[ind].n_bonds += bond_size;
  }
}

/****************** REQ_GET_PART ************/
void mpi_recv_part(int pnode, int part, Particle *pdata)
{
  COMM_TRACE(fprintf(stderr, "0: issuing GET_PART %d %d\n", pnode, part));
  
  if (pnode == this_node) {
    int index = got_particle(part);
    
    memcpy(pdata, &particles[index], sizeof(Particle));
    pdata->max_bonds = pdata->n_bonds;
    if (pdata->n_bonds > 0) {
      pdata->bonds = malloc(sizeof(int)*pdata->n_bonds);
      memcpy(pdata->bonds, particles[index].bonds,
	     sizeof(int)*particles[index].n_bonds);
    }
    else
      pdata->bonds = NULL;
  }
  else {
    MPI_Status status;
    request[0] = REQ_GET_PART;
    request[1] = part;

    MPI_Bcast(request, 2, MPI_INT, 0, MPI_COMM_WORLD);
    pdata->identity = part;
    MPI_Recv(&pdata->type, 1, MPI_INT, pnode,
	     REQ_GET_PART, MPI_COMM_WORLD, &status);
    MPI_Recv(pdata->p, 3, MPI_DOUBLE, pnode,
	     REQ_GET_PART, MPI_COMM_WORLD, &status);
    MPI_Recv(pdata->i, 3, MPI_INT, pnode,
	     REQ_GET_PART, MPI_COMM_WORLD, &status);
    MPI_Recv(&pdata->q, 1, MPI_DOUBLE, pnode,
	     REQ_GET_PART, MPI_COMM_WORLD, &status);
    MPI_Recv(pdata->v, 3, MPI_DOUBLE, pnode,
	     REQ_GET_PART, MPI_COMM_WORLD, &status);
    MPI_Recv(pdata->f, 3, MPI_DOUBLE, pnode,
	     REQ_GET_PART, MPI_COMM_WORLD, &status);
    MPI_Recv(&pdata->n_bonds, 1, MPI_INT, pnode,
	     REQ_GET_PART, MPI_COMM_WORLD, &status);
    pdata->max_bonds = pdata->n_bonds;
    if (pdata->n_bonds > 0) {
      pdata->bonds = malloc(sizeof(int)*pdata->n_bonds);
      MPI_Recv(pdata->bonds, pdata->n_bonds, MPI_INT, pnode,
	       REQ_GET_PART, MPI_COMM_WORLD, &status);
    }
    else
      pdata->bonds = NULL;
  }
}

void mpi_recv_part_slave(int part)
{
  int index;
  Particle *pdata;

  if ((index = got_particle(part)) == -1)
    return;
  pdata = &particles[index];
  
  MPI_Send(&pdata->type, 1, MPI_INT, 0, REQ_GET_PART,
	   MPI_COMM_WORLD);
  MPI_Send(pdata->p, 3, MPI_DOUBLE, 0, REQ_GET_PART,
	   MPI_COMM_WORLD);
  MPI_Send(pdata->i, 3, MPI_INT, 0, REQ_GET_PART,
	   MPI_COMM_WORLD);
  MPI_Send(&pdata->q, 1, MPI_DOUBLE, 0, REQ_GET_PART,
	   MPI_COMM_WORLD);
  MPI_Send(pdata->v, 3, MPI_DOUBLE, 0, REQ_GET_PART,
	   MPI_COMM_WORLD);
  MPI_Send(pdata->f, 3, MPI_DOUBLE, 0, REQ_GET_PART,
	   MPI_COMM_WORLD);
  MPI_Send(&pdata->n_bonds, 1, MPI_INT, 0, REQ_GET_PART,
	   MPI_COMM_WORLD);
  if (pdata->n_bonds > 0)
    MPI_Send(pdata->bonds, pdata->n_bonds, MPI_INT, 0, REQ_GET_PART,
	     MPI_COMM_WORLD);
}

/*********************+ REQ_INTEGRATE ********/
void mpi_integrate(int n_steps)
{
  int task;
  request[0] = REQ_INTEGRATE;
  request[1] = n_steps;

  COMM_TRACE(fprintf(stderr, "0: issuing INTEGRATE %d\n", n_steps));
  MPI_Bcast(request, 2, MPI_INT, this_node, MPI_COMM_WORLD);

  task = request[1];
  if(task==0) {
    /* initialize integrator */
    integrate_vv_init();
  }
  else if(task > 0)
    /* => task = number of steps */
    integrate_vv(task);
  else if(task < 0)
    integrate_vv_exit();
}

void mpi_integrate_slave(int task)
{
  if(task==0) {
    integrate_vv_init();
  }
  else if(task > 0)
    integrate_vv(task);
  else if(task < 0)
    integrate_vv_exit();
  
  COMM_TRACE(fprintf(stderr, "%d: integration task %d done.\n",
		     this_node, task));
}

/*************** REQ_BCAST_IA ************/
void mpi_bcast_ia_params(int i, int j)
{
  request[0] = REQ_BCAST_IA;
  request[1] = i;

  COMM_TRACE(fprintf(stderr, "0: issuing BCAST_IA %d %d\n", i, j));
  MPI_Bcast(request, 2, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&j,  1, MPI_INT, 0, MPI_COMM_WORLD);

  if(j>=0) { /* non-bonded interaction parameters */
    /* INCOMPATIBLE WHEN NODES USE DIFFERENT ARCHITECTURES */
    MPI_Bcast(get_ia_param(i, j), sizeof(IA_parameters), MPI_BYTE,
	      0, MPI_COMM_WORLD);
  }
  else { /* bonded interaction parameters */
    /* INCOMPATIBLE WHEN NODES USE DIFFERENT ARCHITECTURES */
    MPI_Bcast(&(bonded_ia_params[i]),sizeof(Bonded_ia_parameters), MPI_BYTE,
	      0, MPI_COMM_WORLD);
  }
}

void mpi_bcast_ia_params_slave(int i)
{
  int j;
  MPI_Bcast(&j,  1, MPI_INT, 0, MPI_COMM_WORLD);

  if(j>=0) { /* non-bonded interaction parameters */
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
  request[0] = REQ_BCAST_IA_SIZE;
  request[1] = ns;

  COMM_TRACE(fprintf(stderr, "0: issuing BCAST_IA_SIZE %d\n", ns));
  MPI_Bcast(request, 2, MPI_INT, 0, MPI_COMM_WORLD);
}

void mpi_bcast_n_particle_types_slave(int ns)
{
  realloc_ia_params(ns);
}

/*************** REQ_GATHER ************/
void mpi_gather_stats(int job, void *result)
{
  int tot_size, i, pnode;
  int *sizes;
  float *cdata;
  float_packed_particle_data *fdata;
  MPI_Status status;

  request[0] = REQ_GATHER;
  request[1] = job;

  COMM_TRACE(fprintf(stderr, "0: issuing GATHER %d\n", job));
  MPI_Bcast(request, 2, MPI_INT, 0, MPI_COMM_WORLD);

  switch (job) {
  case 0:
    MPI_Gather(&minimum_part_dist, 1, MPI_DOUBLE, (double *)result,
	       1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    break;
  case 1:
    /* similar to REQ_WHO_HAS */
    sizes = malloc(sizeof(int)*n_nodes);
    fdata = (float_packed_particle_data *)result;

    /* first collect number of particles on each node */
    MPI_Gather(&n_particles, 1, MPI_INT, sizes, 1, MPI_INT,
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
	  for (i = 0; i < n_particles; i++) {
	    cdata[3*i    ] = particles[i].p[0];
	    cdata[3*i + 1] = particles[i].p[1];
	    cdata[3*i + 2] = particles[i].p[2];
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
  int i;
  float *cdata;

  switch (job) {
  case 0:
    MPI_Gather(&minimum_part_dist, 1, MPI_DOUBLE, NULL,
	       1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    break;
  case 1:
    /* similar to REQ_WHO_HAS */
    /* first collect number of particles on each node */
    MPI_Gather(&n_particles, 1, MPI_INT, NULL, 1, MPI_INT,
	       0, MPI_COMM_WORLD);
    cdata = malloc(3*n_particles*sizeof(float));
    for (i = 0; i < n_particles; i++) {
      cdata[3*i    ] = particles[i].p[0];
      cdata[3*i + 1] = particles[i].p[1];
      cdata[3*i + 2] = particles[i].p[2];
    }
    MPI_Send(cdata, 3*n_particles, MPI_FLOAT, 0, REQ_GATHER,
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
