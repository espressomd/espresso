#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "communication.h"
#include "integrate.h"
#include "global.h"
#include "debug.h"

int this_node = -1;
int nprocs = -1;

/**********************************************
 * slave callbacks.
 **********************************************/

SlaveCallback *callbacks[] = 
{
  mpi_stop_slave,                /*  0: REQ_TERM */
  mpi_bcast_parameter_slave,     /*  1: REQ_BCAST_PAR */
  mpi_who_has_slave,             /*  2: REQ_WHO_HAS */
  mpi_attach_particle_slave,     /*  3: REQ_ATTACH */
  mpi_send_pos_slave,            /*  4: REQ_SET_POS */
  NULL,                          /*  5: REQ_SET_V */
  NULL,                          /*  6: REQ_SET_F */
  mpi_send_q_slave,              /*  7: REQ_SET_Q */
  mpi_send_type_slave,           /*  8: REQ_SET_TYPE */
  mpi_recv_part_slave,           /*  9: REQ_GET_PART */
  mpi_integrate_slave,           /* 10: REQ_INTEGRATE */
  mpi_bcast_ia_params_slave      /* 11: REQ_BCAST_IA */ 
};


/**********************************************
 * procedures
 **********************************************/

void mpi_init(int *argc, char ***argv)
{
  MPI_Init(argc, argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &this_node);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
}

/**************** REQ_TERM ************/

int terminated = 0;

void mpi_stop()
{
  int req[2] = { REQ_TERM, 0 };

  if (terminated)
    return;

  MPI_Bcast(req, 2, MPI_INT, 0, MPI_COMM_WORLD);

  COMM_TRACE(fprintf(stderr, "%d: sent TERM\n", this_node));

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  terminated = 1;
}

void mpi_stop_slave(int param)
{
  COMM_TRACE(fprintf(stderr, "%d: exiting\n", this_node));

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  exit(-1);
}

/*************** REQ_BCAST_PAR ************/
void mpi_bcast_parameter(int i)
{
  int req[2] = { REQ_BCAST_PAR, i };

  /* start broadcasting. */
  MPI_Bcast(req, 2, MPI_INT, 0, MPI_COMM_WORLD);

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
int mpi_who_has(int part)
{
  int req[2]   = { REQ_WHO_HAS, part };
  int  sendbuf = got_particle(part);
  int *recvbuf = malloc(sizeof(int)*nprocs);
  int pnode;

  COMM_TRACE(fprintf(stderr, "ask for part %d %d\n", part, sendbuf));

  MPI_Bcast(req, 2, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&sendbuf, 1, MPI_INT,
	     recvbuf, 1, MPI_INT,
	     0, MPI_COMM_WORLD);

  for (pnode = 0; pnode < nprocs; pnode++) {
    COMM_TRACE(fprintf(stderr, "pnode %d answered %d\n",
		       pnode, recvbuf[pnode]));
    if (recvbuf[pnode] != -1)
      break;
  }
  free(recvbuf);
  
  return (pnode == nprocs) ? -1 : pnode;
}

void mpi_who_has_slave(int ident)
{
  int sendbuf = got_particle(ident);
  MPI_Gather(&sendbuf,  1, MPI_INT,
	     NULL, nprocs, MPI_INT,
	     0, MPI_COMM_WORLD);
  COMM_TRACE(fprintf(stderr, "%d: i got %d %d\n", this_node, ident, sendbuf));
}

/**************** REQ_ATTACH ***********/
void mpi_attach_particle(int part, int pnode)
{
  int req[2] = { REQ_ATTACH, pnode };

  COMM_TRACE(fprintf(stderr, "attach req %d %d\n", part, pnode));

  MPI_Bcast(req, 2, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&part, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (pnode == this_node)
    add_particle(part);

  if (part >= n_total_particles)
    n_total_particles = part + 1;
}

void mpi_attach_particle_slave(int pnode)
{
  int part;

  MPI_Bcast(&part, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  if (pnode == this_node)
    add_particle(part);

  if (part >= n_total_particles)
    n_total_particles = part + 1;
}

/****************** REQ_SET_POS ************/
void mpi_send_pos(int pnode, int part, double p[3])
{
  if (pnode == this_node) {
    int index = got_particle(part);
    particles[index].p[0] = p[0];
    particles[index].p[1] = p[1];
    particles[index].p[2] = p[2];
  }
  else {
    int req[2] = { REQ_SET_POS, part };
    MPI_Bcast(req, 2, MPI_INT, 0, MPI_COMM_WORLD);
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

/****************** REQ_SET_POS ************/
void mpi_recv_part(int pnode, int part, Particle *pdata)
{
  COMM_TRACE(fprintf(stderr, "part req %d %d\n", pnode, part));
  
  if (pnode == this_node) {
    int index = got_particle(part);
    
    memcpy(pdata, &particles[index], sizeof(Particle));
    pdata->max_bonds = pdata->n_bonds;
    pdata->bonds = malloc(sizeof(int)*pdata->n_bonds);
    memcpy(pdata->bonds, particles[index].bonds,
	   sizeof(int)*particles[index].n_bonds);
  }
  else {
    int req[2] = { REQ_GET_PART, part };
    MPI_Status status;

    MPI_Bcast(req, 2, MPI_INT, 0, MPI_COMM_WORLD);
    pdata->identity = part;
    MPI_Recv(&pdata->type, 1, MPI_INT, pnode,
	     REQ_GET_PART, MPI_COMM_WORLD, &status);
    MPI_Recv(pdata->p, 3, MPI_DOUBLE, pnode,
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
    pdata->bonds = malloc(sizeof(int)*pdata->n_bonds);
    MPI_Recv(pdata->bonds, pdata->n_bonds, MPI_INT, pnode,
	     REQ_GET_PART, MPI_COMM_WORLD, &status);
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
  MPI_Send(&pdata->q, 1, MPI_DOUBLE, 0, REQ_GET_PART,
	   MPI_COMM_WORLD);
  MPI_Send(pdata->v, 3, MPI_DOUBLE, 0, REQ_GET_PART,
	   MPI_COMM_WORLD);
  MPI_Send(pdata->f, 3, MPI_DOUBLE, 0, REQ_GET_PART,
	   MPI_COMM_WORLD);
  MPI_Send(&pdata->n_bonds, 1, MPI_INT, 0, REQ_GET_PART,
	   MPI_COMM_WORLD);
  MPI_Send(pdata->bonds, pdata->n_bonds, MPI_INT, 0, REQ_GET_PART,
	   MPI_COMM_WORLD);
}

/*********************+ REQ_INTEGRATE ********/
void mpi_integrate(int n_steps)
{
  int req[2] = { REQ_INTEGRATE, n_steps };
  int task;

  MPI_Bcast(req, 2, MPI_INT, this_node, MPI_COMM_WORLD);

  task = req[1];
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

/********************* REQ_SET_Q ********/
void mpi_send_q(int pnode, int part, double q)
{
  if (pnode == this_node) {
    int index = got_particle(part);
    particles[index].q = q;
  }
  else {
    int req[2] = { REQ_SET_Q, part };
    MPI_Bcast(req, 2, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Send(&q, 1, MPI_DOUBLE, pnode, REQ_SET_POS, MPI_COMM_WORLD);
  }
}

void mpi_send_q_slave(int part)
{
  int index = got_particle(part);
  MPI_Status status;
  if (index != -1) {
    MPI_Recv(&particles[index].q, 1, MPI_DOUBLE, 0, REQ_SET_POS,
	     MPI_COMM_WORLD, &status);
  }
}

/********************* REQ_SET_TYPE ********/
void mpi_send_type(int pnode, int part, int type)
{
  if (pnode == this_node) {
    int index = got_particle(part);
    particles[index].type = type;
  }
  else {
    int req[2] = { REQ_SET_Q, part };
    MPI_Bcast(req, 2, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Send(&type, 1, MPI_INT, pnode, REQ_SET_POS, MPI_COMM_WORLD);
  }
}

void mpi_send_type_slave(int part)
{
  int index = got_particle(part);
  MPI_Status status;
  if (index != -1) {
    MPI_Recv(&particles[index].type, 1, MPI_INT, 0, REQ_SET_POS,
	     MPI_COMM_WORLD, &status);
  }
}

/*************** REQ_BCAST_IA ************/
void mpi_bcast_ia_params(int i, int j)
{
  int req[2] = { REQ_BCAST_IA, i };

  /* start broadcasting. */
  MPI_Bcast(req, 2, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&j,  1, MPI_INT, 0, MPI_COMM_WORLD);
  /* INCOMPATIBLE WHEN NODES USE DIFFERENT ARCHITECTURES */
  MPI_Bcast(get_ia_param(i, j), sizeof(IA_parameters), MPI_BYTE,
	    0, MPI_COMM_WORLD);
}

void mpi_bcast_ia_params_slave(int i)
{
  int j;
  MPI_Bcast(&j,  1, MPI_INT, 0, MPI_COMM_WORLD);

  /* INCOMPATIBLE WHEN NODES USE DIFFERENT ARCHITECTURES */
  MPI_Bcast(safe_get_ia_param(i, j), sizeof(IA_parameters), MPI_BYTE,
	    0, MPI_COMM_WORLD);
}

/*********************** MAIN LOOP for slaves ****************/

void mpi_loop()
{
  int req[2];

  for (;;) {
    MPI_Bcast(req, 2, MPI_INT, 0, MPI_COMM_WORLD);
    COMM_TRACE(fprintf(stderr, "%d: processing request %d %d...", this_node,
		       req[0], req[1]));

    if ((req[0] < 0) || (req[0] >= REQ_MAXIMUM))
      continue;
    callbacks[req[0]](req[1]);
    COMM_TRACE(fprintf(stderr, "%d: finished\n", this_node));
  }
}
