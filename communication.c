#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "communication.h"
#include "interaction_data.h"
#include "integrate.h"
#include "global.h"
#include "debug.h"

int this_node = -1;
int nprocs = -1;

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

#define REQ_TERM      0
#define REQ_BCAST_PAR 1
#define REQ_WHO_HAS   2
#define REQ_ATTACH    3
#define REQ_SET_POS   4
#define REQ_SET_V     5
#define REQ_SET_F     6
#define REQ_SET_Q     7
#define REQ_SET_TYPE  8
#define REQ_GET_PART  9
#define REQ_INTEGRATE 10
#define REQ_BCAST_IA  11
#define REQ_MAXIMUM   12

void mpi_stop_slave(int parm);
void mpi_bcast_parameter_slave(int parm);
void mpi_who_has_slave(int parm);
void mpi_attach_particle_slave(int parm);
void mpi_send_pos_slave(int parm);
void mpi_send_v_slave(int parm);
void mpi_send_f_slave(int parm);
void mpi_send_q_slave(int parm);
void mpi_send_type_slave(int parm);
void mpi_recv_part_slave(int parm);
void mpi_integrate_slave(int parm);
void mpi_bcast_ia_params_slave(int i);

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
  mpi_recv_part_slave,           /*  9: REQ_GET_PART */
  mpi_integrate_slave,           /* 10: REQ_INTEGRATE */
  mpi_bcast_ia_params_slave      /* 11: REQ_BCAST_IA */ 
};

char *names[] = {
  "TERM",      /*  0 */
  "BCAST_PAR", /*  1 */
  "WHO_HAS",   /*  2 */
  "ATTACH",    /*  3 */
  "SET_POS",   /*  4 */
  "SET_V",     /*  5 */
  "SET_F",     /*  6 */
  "SET_Q",     /*  7 */
  "SET_TYPE",  /*  8 */
  "GET_PART",  /*  9 */
  "INTEGRATE", /* 10 */
  "BCAST_IA"   /* 11 */
};

/**********************************************
 * procedures
 **********************************************/

#ifdef MPI_CORE
void mpi_core(MPI_Comm *comm, int *errcode,...) {
  fprintf(stderr, "Aborting due to MPI error %d, forcing core dump\n", *errcode);
  *(int *)0 = 0;  
}
#endif

void mpi_init(int *argc, char ***argv)
{
#ifdef MPI_CORE
  MPI_Errhandler mpi_errh;
#endif

  MPI_Init(argc, argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &this_node);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

#ifdef MPI_CORE
  MPI_Errhandler_create((MPI_Handler_function *)mpi_core, &mpi_errh);
  MPI_Errhandler_set(MPI_COMM_WORLD, mpi_errh);
#endif

}

/**************** REQ_TERM ************/

int terminated = 0;

void mpi_stop()
{
  int req[2] = { REQ_TERM, 0 };

  if (terminated)
    return;

  COMM_TRACE(fprintf(stderr, "0: issuing TERM\n"));
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
  exit(0);
}

/*************** REQ_BCAST_PAR ************/
void mpi_bcast_parameter(int i)
{
  int req[2];
  req[0] = REQ_BCAST_PAR;
  req[1] = i;

  COMM_TRACE(fprintf(stderr, "0: issuing BCAST_PAR %d\n", i));
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
void mpi_who_has()
{
  int req[2] = { REQ_WHO_HAS, 0};
  int *sizes = malloc(sizeof(int)*nprocs);
  int *pdata = NULL;
  int pdata_s = 0, i;
  int pnode;
  MPI_Status status;

  COMM_TRACE(fprintf(stderr, "0: issuing WHO_HAS\n"));

  MPI_Bcast(req, 2, MPI_INT, 0, MPI_COMM_WORLD);

  /* first collect number of particles on each node */
  MPI_Gather(&n_particles, 1, MPI_INT, sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);

  /* then fetch particle locations */
  for (pnode = 0; pnode < nprocs; pnode++) {
    COMM_TRACE(fprintf(stderr, "node %d reports %d particles\n",
		       pnode, sizes[pnode]));
    if (pnode == this_node) {
      for (i = 0; i < n_particles; i++) {
	if (particles[i].identity != -1)
	  particle_node[particles[i].identity] = pnode;
      }
    }
    else {
      if (pdata_s < sizes[pnode]) {
	pdata_s = sizes[pnode];
	pdata = realloc(pdata, sizeof(int)*pdata_s);
      }
      MPI_Recv(pdata, sizes[pnode], MPI_INT, pnode, REQ_WHO_HAS,
	       MPI_COMM_WORLD, &status);
      for (i = 0; i < sizes[pnode]; i++)
	particle_node[particles[i].identity] = pnode;
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
  sendbuf = malloc(sizeof(int)*n_particles);
  npart = 0;
  for (i = 0; i < n_particles; i++)
    if (particles[i].identity != -1)
      sendbuf[npart++] = particles[i].identity;
  MPI_Send(sendbuf, npart, MPI_INT, 0, REQ_WHO_HAS, MPI_COMM_WORLD);
  free(sendbuf);
}

/**************** REQ_ATTACH ***********/
void mpi_attach_particle(int part, int pnode)
{
  int req[2];
  req[0] = REQ_ATTACH;
  req[1] = pnode;

  COMM_TRACE(fprintf(stderr, "0: issuing ATTACH %d %d\n", pnode, part));

  MPI_Bcast(req, 2, MPI_INT, 0, MPI_COMM_WORLD);
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
    int req[2];
    req[0] = REQ_SET_POS;
    req[1] = part;
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
    int req[2];
    req[0] = REQ_SET_V;
    req[1] = part;
    MPI_Bcast(req, 2, MPI_INT, 0, MPI_COMM_WORLD);
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
    int req[2];
    req[0] = REQ_SET_F;
    req[1] = part;
    MPI_Bcast(req, 2, MPI_INT, 0, MPI_COMM_WORLD);
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
    int req[2];
    req[0] = REQ_SET_Q;
    req[1] = part;
    MPI_Bcast(req, 2, MPI_INT, 0, MPI_COMM_WORLD);
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
    int req[2];
    req[0] = REQ_SET_TYPE;
    req[1] = part;
    MPI_Bcast(req, 2, MPI_INT, 0, MPI_COMM_WORLD);
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

/****************** REQ_GET_PART ************/
void mpi_recv_part(int pnode, int part, Particle *pdata)
{
  COMM_TRACE(fprintf(stderr, "0: issuing GET_PART %d %d\n", pnode, part));
  
  if (pnode == this_node) {
    int index = got_particle(part);
    
    memcpy(pdata, &particles[index], sizeof(Particle));
    pdata->max_bonds = pdata->n_bonds;
    pdata->bonds = malloc(sizeof(int)*pdata->n_bonds);
    memcpy(pdata->bonds, particles[index].bonds,
	   sizeof(int)*particles[index].n_bonds);
  }
  else {
    int req[2];
    MPI_Status status;
    req[0] = REQ_GET_PART;
    req[1] = part;

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
  int req[2];
  int task;
  req[0] = REQ_INTEGRATE;
  req[1] = n_steps;

  COMM_TRACE(fprintf(stderr, "0: issuing INTEGRATE %d\n", n_steps));
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

/*************** REQ_BCAST_IA ************/
void mpi_bcast_ia_params(int i, int j)
{
  int req[2];
  req[0] = REQ_BCAST_IA;
  req[1] = i;

  COMM_TRACE(fprintf(stderr, "0: issuing BCAST_IA %d %d\n", i, j));
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
    COMM_TRACE(fprintf(stderr, "%d: processing %s %d...\n", this_node,
		       names[req[0]], req[1]));
 
    if ((req[0] < 0) || (req[0] >= REQ_MAXIMUM))
      continue;
    callbacks[req[0]](req[1]);
    COMM_TRACE(fprintf(stderr, "%d: finished %s %d\n", this_node,
		       names[req[0]], req[1]));
  }
}
