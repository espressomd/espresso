#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "communication.h"
#include "integrate.h"
#include "global.h"

//#define DEBUG

int this_node = -1;
int nprocs = -1;
int processor_grid[3] = { -1, -1, -1};
int neighbors[6] = {0, 0, 0, 0, 0, 0};

void init_mpi(int *argc, char ***argv)
{
  MPI_Init(argc, argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &this_node);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
}

void stop_mpi()
{
  int req[2] = { REQ_TERM, 0 };
  MPI_Bcast(req, 2, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Finalize();
}

void mpi_set_system_parameters()
{
  int req[2] = { REQ_SYS_PAR, 0 };
  int i;

  /* calculate dependend parameters. */
  for(i=0;i<3;i++) local_box_l[i] = box_l[i]/(double)processor_grid[i]; 
  /* start broadcasting. */
  MPI_Bcast(req, 2, MPI_INT, 0, MPI_COMM_WORLD); 
  /* broadcast system parameters. */
  MPI_Bcast(&n_total_particles, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(box_l, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
  MPI_Bcast(local_box_l, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
  MPI_Bcast(processor_grid, 3, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast(&n_particle_types, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&n_interaction_types, 1, MPI_INT, 0, MPI_COMM_WORLD);
  return;
}

int mpi_who_has(int part)
{
  int req[2]   = { REQ_WHO_HAS, part };
  int  sendbuf = got_particle(part);
  int *recvbuf = malloc(sizeof(int)*nprocs);
  int pnode;

#ifdef DEBUG
  fprintf(stderr, "ask for part %d %d\n", part, sendbuf);
#endif
  MPI_Bcast(req, 2, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&sendbuf, 1, MPI_INT,
	     recvbuf, 1, MPI_INT,
	     0, MPI_COMM_WORLD);

  for (pnode = 0; pnode < nprocs; pnode++) {
#ifdef DEBUG
    fprintf(stderr, "pnode %d answered %d\n", pnode, recvbuf[pnode]);
#endif
    if (recvbuf[pnode] != -1)
      break;
  }
  free(recvbuf);
  
  return (pnode == nprocs) ? -1 : pnode;
}

void mpi_attach_particle(int part, int pnode)
{
#ifdef DEBUG
  fprintf(stderr, "attach req %d %d\n", part, pnode);
#endif
  if (pnode == this_node)
    add_particle(part);
  else {
    int req[2] = { REQ_ATTACH, pnode };
    MPI_Bcast(req, 2, MPI_INT, this_node, MPI_COMM_WORLD);
    MPI_Send(&part, 1, MPI_INT, pnode, REQ_ATTACH, MPI_COMM_WORLD);
  }
}

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

void mpi_recv_part(int pnode, int part, Particle *pdata)
{
  if (pnode == this_node) {
    int index = got_particle(part);
    memcpy(pdata, &particles[index], sizeof(Particle));
    // for now, until we have a deep copy...
    pdata->n_bond   = 0;
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
    // for now, until we have a deep copy...
    pdata->n_bond   = 0;
  }
}


void mpi_send_pdata(Particle *pdata)
{
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
}

void mpi_integrate(int n_steps)
{
  int req[2] = { REQ_INTEGRATE, n_steps };
  int task;
  int *recvbuf = malloc(sizeof(int)*nprocs);

  MPI_Bcast(req, 2, MPI_INT, 0, MPI_COMM_WORLD);

  task = req[1];
  if(task==0) {
    /* initialize integrator */
    /* send some required variables */
    MPI_Bcast(box_l, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
    MPI_Bcast(processor_grid, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
    integrate_vv_init();
  }
  else if(task > 0)
    /* => task = number of steps */
    integrate_vv(task);
  else if(task < 0)
    integrate_vv_exit();
}

void mpi_loop()
{
  int    req[2];
  int    sendbuf;
  int    part_num, index;
  int    task;
  MPI_Status status;

  for (;;) {
    MPI_Bcast(req, 2, MPI_INT, 0, MPI_COMM_WORLD);
    switch (req[0]) {
    case REQ_TERM:
      MPI_Finalize();
      return;
    case REQ_SYS_PAR:
      MPI_Bcast(&n_total_particles, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(box_l, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
      MPI_Bcast(local_box_l, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
      MPI_Bcast(processor_grid, 3, MPI_INT, 0, MPI_COMM_WORLD); 
      MPI_Bcast(&n_particle_types, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&n_interaction_types, 1, MPI_INT, 0, MPI_COMM_WORLD);
      break;
    case REQ_WHO_HAS:
      sendbuf = got_particle(req[1]);
      MPI_Gather(&sendbuf,  1, MPI_INT,
		 NULL, nprocs, MPI_INT,
		 0, MPI_COMM_WORLD);
#ifdef DEBUG
      fprintf(stderr, "%d: answered part %d %d\n", this_node, req[1], sendbuf);
#endif
      break;
    case REQ_ATTACH:
      if (req[1] == this_node) {
	MPI_Recv(&part_num, 1, MPI_INT, 0, REQ_ATTACH,
		 MPI_COMM_WORLD, &status);
	add_particle(part_num);
      }
      break;
    case REQ_SET_POS:
      index = got_particle(req[1]);
      if (index != -1) {
	MPI_Recv(particles[index].p, 3, MPI_DOUBLE, 0, REQ_SET_POS,
		 MPI_COMM_WORLD, &status);
      }
      break;
    case REQ_GET_PART:
      index = got_particle(req[1]);
      if (index != -1)
	mpi_send_pdata(&particles[index]);
      break;
    case REQ_INTEGRATE:
      task = req[1];
      if(task==0) {

	/* the local box_l would be also fine */
	MPI_Bcast(box_l, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
	MPI_Bcast(processor_grid, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD); 

	integrate_vv_init();
      }
      else if(task > 0)
	integrate_vv(task);
      else if(task < 0)
	integrate_vv_exit();

#ifdef DEBUG
      fprinf(stderr, "%d: integration task %d done.\n",
	     this_node, req[1]);
#endif
      break;
    default:
#ifdef DEBUG
      fprintf(stderr, "%d: got unknown request %d\n", this_node, req[0]);
#endif
      break;
    }
  }
}
