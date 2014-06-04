/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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
#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "utils.hpp"
#include "communication.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "integrate.hpp"
#include "cells.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "initialize.hpp"
#include "forces.hpp"
#include "rotation.hpp"
#include "p3m.hpp"
#include "statistics.hpp"
#include "energy.hpp"
#include "pressure.hpp"
#include "random.hpp"
#include "lj.hpp"
#include "lb.hpp"
#include "lb-boundaries.hpp"
#include "morse.hpp"
#include "buckingham.hpp"
#include "tab.hpp"
#include "overlap.hpp"
#include "ljcos.hpp"
#include "ljangle.hpp"
#include "gb.hpp"
#include "mmm1d.hpp"
#include "mmm2d.hpp"
#include "maggs.hpp"
#include "elc.hpp"
#include "iccp3m.hpp"
#include "statistics_chain.hpp"
#include "statistics_fluid.hpp"
#include "virtual_sites.hpp"
#include "topology.hpp"
#include "errorhandling.hpp"
#include "molforces.hpp"
#include "mdlc_correction.hpp"
#include "reaction.hpp"
#include "galilei.hpp"
#include "external_potential.hpp"
#include "statistics_correlation.hpp"
#include "cuda_interface.hpp"
#include "EspressoSystemInterface.hpp"
#include "statistics_observable.hpp"

int this_node = -1;
int n_nodes = -1;
MPI_Comm comm_cart;
/**********************************************
 * slave callbacks.
 **********************************************/
typedef void (SlaveCallback)(int node, int param);

// if you want to add a callback, add it here, and here only
#define CALLBACK_LIST \
  CB(mpi_stop_slave) \
  CB(mpi_bcast_parameter_slave) \
  CB(mpi_who_has_slave) \
  CB(mpi_bcast_event_slave) \
  CB(mpi_place_particle_slave) \
  CB(mpi_send_v_slave) \
  CB(mpi_send_f_slave) \
  CB(mpi_send_q_slave) \
  CB(mpi_send_type_slave) \
  CB(mpi_send_bond_slave) \
  CB(mpi_recv_part_slave) \
  CB(mpi_integrate_slave) \
  CB(mpi_bcast_ia_params_slave) \
  CB(mpi_bcast_n_particle_types_slave) \
  CB(mpi_gather_stats_slave) \
  CB(mpi_set_time_step_slave) \
  CB(mpi_get_particles_slave) \
  CB(mpi_bcast_coulomb_params_slave) \
  CB(mpi_bcast_collision_params_slave) \
  CB(mpi_send_ext_force_slave) \
  CB(mpi_send_ext_torque_slave) \
  CB(mpi_place_new_particle_slave) \
  CB(mpi_remove_particle_slave) \
  CB(mpi_bcast_constraint_slave) \
  CB(mpi_random_seed_slave) \
  CB(mpi_random_stat_slave) \
  CB(mpi_cap_forces_slave) \
  CB(mpi_bit_random_seed_slave) \
  CB(mpi_bit_random_stat_slave) \
  CB(mpi_get_constraint_force_slave) \
  CB(mpi_rescale_particles_slave) \
  CB(mpi_bcast_cell_structure_slave) \
  CB(mpi_send_quat_slave) \
  CB(mpi_send_omega_slave) \
  CB(mpi_send_torque_slave) \
  CB(mpi_send_mol_id_slave) \
  CB(mpi_bcast_nptiso_geom_slave) \
  CB(mpi_update_mol_ids_slave) \
  CB(mpi_sync_topo_part_info_slave) \
  CB(mpi_send_mass_slave) \
  CB(mpi_send_solvation_slave) \
  CB(mpi_gather_runtime_errors_slave) \
  CB(mpi_send_exclusion_slave) \
  CB(mpi_bcast_lb_params_slave) \
  CB(mpi_bcast_cuda_global_part_vars_slave) \
  CB(mpi_send_dip_slave) \
  CB(mpi_send_dipm_slave) \
  CB(mpi_send_fluid_slave) \
  CB(mpi_recv_fluid_slave) \
  CB(mpi_local_stress_tensor_slave) \
  CB(mpi_send_virtual_slave) \
  CB(mpi_iccp3m_iteration_slave) \
  CB(mpi_iccp3m_init_slave) \
  CB(mpi_send_rotational_inertia_slave) \
  CB(mpi_bcast_lbboundary_slave) \
  CB(mpi_send_mu_E_slave) \
  CB(mpi_bcast_max_mu_slave) \
  CB(mpi_send_vs_relative_slave) \
  CB(mpi_recv_fluid_populations_slave) \
  CB(mpi_send_fluid_populations_slave) \
  CB(mpi_recv_fluid_boundary_flag_slave) \
  CB(mpi_set_particle_temperature_slave) \
  CB(mpi_set_particle_gamma_slave) \
  CB(mpi_kill_particle_motion_slave) \
  CB(mpi_kill_particle_forces_slave) \
  CB(mpi_system_CMS_slave) \
  CB(mpi_system_CMS_velocity_slave) \
  CB(mpi_galilei_transform_slave) \
  CB(mpi_setup_reaction_slave) \
  CB(mpi_send_rotation_slave) \
  CB(mpi_external_potential_broadcast_slave) \
  CB(mpi_external_potential_tabulated_read_potential_file_slave) \
  CB(mpi_external_potential_sum_energies_slave) \
  CB(mpi_observable_lb_radial_velocity_profile_slave) \

// create the forward declarations
#define CB(name) void name(int node, int param);
CALLBACK_LIST

// create the list of callbacks
#undef CB
#define CB(name) name,
static SlaveCallback *slave_callbacks[] = {
  CALLBACK_LIST
};

const int N_CALLBACKS = sizeof(slave_callbacks)/sizeof(SlaveCallback*);

// create the list of names
#undef CB
#define CB(name) #name,

const char *names[] = {
  CALLBACK_LIST
};

// tag which is used by MPI send/recv inside the slave functions
#define SOME_TAG 42


/** The requests are compiled statically here, so that after a crash
    you can get the last issued request from the debugger. */ 
static int request[3];

/**********************************************
 * procedures
 **********************************************/

#ifdef MPI_CORE
void mpi_core(MPI_Comm *comm, int *errcode,...) {
  int len;
  char string[1024];
  MPI_Error_string(*errcode, string, &len);
  fprintf(stderr, "%d: Aborting due to MPI error %d: %s. Forcing core dump.\n", this_node, *errcode, string);
  core();
}
#endif

void mpi_init(int *argc, char ***argv)
{
#ifdef MPI_CORE
  MPI_Errhandler mpi_errh;
#endif

  MPI_Init(argc, argv);

  MPI_Comm_size(MPI_COMM_WORLD, &n_nodes);

  int periodic[3]={1,1,1}, reorder = 1;
  MPI_Dims_create(n_nodes, 3, node_grid);

  MPI_Cart_create(MPI_COMM_WORLD, 3, node_grid, periodic, reorder, &comm_cart);

  MPI_Comm_rank(comm_cart, &this_node);

  MPI_Cart_coords(comm_cart, this_node, 3, node_pos);

#ifdef MPI_CORE
  MPI_Comm_create_errhandler((MPI_Handler_function *)mpi_core, &mpi_errh);
  MPI_Comm_set_errhandler(comm_cart, mpi_errh);
#endif
}

static void mpi_call(SlaveCallback cb, int node, int param) {
  // find req number in callback array
  int reqcode;
  for (reqcode = 0; reqcode < N_CALLBACKS; reqcode++) {
    if (cb == slave_callbacks[reqcode]) break;
  }

  if (reqcode >= N_CALLBACKS) {
    fprintf(stderr, "%d: INTERNAL ERROR: unknown callback %d called\n", this_node, reqcode);
    errexit();
  }

  request[0] = reqcode;
  request[1] = node;
  request[2] = param;

  COMM_TRACE(fprintf(stderr, "%d: issuing %s %d %d\n",
		     this_node, names[reqcode], node, param));
#ifdef ASYNC_BARRIER
  MPI_Barrier(comm_cart);
#endif
  MPI_Bcast(request, 3, MPI_INT, 0, comm_cart);
  COMM_TRACE(fprintf(stderr, "%d: finished sending.\n", this_node));
}

/**************** REQ_TERM ************/

static int terminated = 0;

void mpi_stop()
{
  if (terminated)
    return;

  mpi_call(mpi_stop_slave, -1, 0);

  MPI_Barrier(comm_cart);
  MPI_Finalize();
  regular_exit = 1;
  terminated = 1;
}

void mpi_stop_slave(int node, int param)
{
  COMM_TRACE(fprintf(stderr, "%d: exiting\n", this_node));

  MPI_Barrier(comm_cart);
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
	      MPI_INT, 0, comm_cart);
    break;
  case TYPE_BOOL:
    MPI_Bcast((int *)fields[i].data, 1,
	      MPI_INT, 0, comm_cart);
    break;
  case TYPE_DOUBLE:
    MPI_Bcast((double *)fields[i].data, fields[i].dimension,
	      MPI_DOUBLE, 0, comm_cart);
    break;
  default: break;
  }

  on_parameter_change(i);
}

int mpi_bcast_parameter(int i)
{
  mpi_call(mpi_bcast_parameter_slave, -1, i);

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
  int *sizes = (int*)malloc(sizeof(int)*n_nodes);
  int *pdata = NULL;
  int pdata_s = 0, i, c;
  int pnode;
  int n_part;

  mpi_call(mpi_who_has_slave, -1, 0);

  n_part = cells_get_n_particles();
  /* first collect number of particles on each node */
  MPI_Gather(&n_part, 1, MPI_INT, sizes, 1, MPI_INT, 0, comm_cart);

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
      MPI_Recv(pdata, sizes[pnode], MPI_INT, pnode, SOME_TAG,
	       comm_cart, MPI_STATUS_IGNORE);
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
  MPI_Gather(&n_part, 1, MPI_INT, NULL, 0, MPI_INT, 0, comm_cart);
  if (n_part == 0)
    return;

  sendbuf = (int*)malloc(sizeof(int)*n_part);
  npart = 0;
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    for (i = 0; i < cell->n; i++)
      sendbuf[npart++] = cell->part[i].p.identity;
  }
  MPI_Send(sendbuf, npart, MPI_INT, 0, SOME_TAG, comm_cart);
  free(sendbuf);
}

/**************** REQ_CHTOPL ***********/
void mpi_bcast_event(int event)
{
  mpi_call(mpi_bcast_event_slave, -1, event);
  mpi_bcast_event_slave(-1, event);
}

void mpi_bcast_event_slave(int node, int event)
{
  switch (event) {
#ifdef ELECTROSTATICS
#ifdef P3M
  case P3M_COUNT_CHARGES:
    p3m_count_charged_particles();
    break;
#endif
  case MAGGS_COUNT_CHARGES:
    maggs_count_charged_particles();
    break; 
#endif
  case SORT_PARTICLES:
    local_sort_particles();
    break;
  case CHECK_PARTICLES:
    check_particles();
    break;

#ifdef DP3M
  case P3M_COUNT_DIPOLES:
    dp3m_count_magnetic_particles();
    break;
#endif

  default:;
  }
}

/****************** REQ_PLACE/REQ_PLACE_NEW ************/

void mpi_place_particle(int pnode, int part, double p[3])
{
  mpi_call(mpi_place_particle_slave, pnode, part);

  if (pnode == this_node)
    local_place_particle(part, p, 0);
  else
    MPI_Send(p, 3, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);

  on_particle_change();
}

void mpi_place_particle_slave(int pnode, int part)
{
  double p[3];
  
  if (pnode == this_node) {
    MPI_Recv(p, 3, MPI_DOUBLE, 0, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
    local_place_particle(part, p, 0);
  }

  on_particle_change();
}

void mpi_place_new_particle(int pnode, int part, double p[3])
{
  mpi_call(mpi_place_new_particle_slave, pnode, part);
  added_particle(part);

  if (pnode == this_node)
    local_place_particle(part, p, 1);
  else
    MPI_Send(p, 3, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);

  on_particle_change();
}


void mpi_place_new_particle_slave(int pnode, int part)
{
  double p[3];
  
  added_particle(part);

  if (pnode == this_node) {
    MPI_Recv(p, 3, MPI_DOUBLE, 0, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
    local_place_particle(part, p, 1);
  }

  on_particle_change();
}

/****************** REQ_SET_V ************/
void mpi_send_v(int pnode, int part, double v[3])
{
  mpi_call(mpi_send_v_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    memcpy(p->m.v, v, 3*sizeof(double));
  }
  else
    MPI_Send(v, 3, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);

  on_particle_change();
}

void mpi_send_v_slave(int pnode, int part)
{
  if (pnode == this_node) {
    Particle *p = local_particles[part];
        MPI_Recv(p->m.v, 3, MPI_DOUBLE, 0, SOME_TAG,
	     comm_cart, MPI_STATUS_IGNORE);
  }

  on_particle_change();
}

/****************** REQ_SET_F ************/
void mpi_send_f(int pnode, int part, double F[3])
{
  mpi_call(mpi_send_f_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    memcpy(p->f.f, F, 3*sizeof(double));
  }
  else
    MPI_Send(F, 3, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);

  on_particle_change();
}

void mpi_send_f_slave(int pnode, int part)
{
  if (pnode == this_node) {
    Particle *p = local_particles[part];
        MPI_Recv(p->f.f, 3, MPI_DOUBLE, 0, SOME_TAG,
	     comm_cart, MPI_STATUS_IGNORE);
  }

  on_particle_change();
}

/********************* REQ_SET_Q ********/
void mpi_send_q(int pnode, int part, double q)
{
#ifdef ELECTROSTATICS
  mpi_call(mpi_send_q_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->p.q = q;
  }
  else {
    MPI_Send(&q, 1, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
#endif
}

void mpi_send_q_slave(int pnode, int part)
{
#ifdef ELECTROSTATICS
  if (pnode == this_node) {
    Particle *p = local_particles[part];
        MPI_Recv(&p->p.q, 1, MPI_DOUBLE, 0, SOME_TAG,
	     comm_cart, MPI_STATUS_IGNORE);
  }

  on_particle_change();
#endif
}

/********************* REQ_SET_MU_E ********/
void mpi_send_mu_E(int pnode, int part, double mu_E[3])
{
#ifdef LB_ELECTROHYDRODYNAMICS
  mpi_call(mpi_send_mu_E_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->p.mu_E[0] = mu_E[0];
    p->p.mu_E[1] = mu_E[1];
    p->p.mu_E[2] = mu_E[2];
  }
  else {
    MPI_Send(&mu_E, 3, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
#endif
}

void mpi_send_mu_E_slave(int pnode, int part)
{
#ifdef LB_ELECTROHYDRODYNAMICS
  if (pnode == this_node) {
    Particle *p = local_particles[part];
        MPI_Recv(&p->p.mu_E, 3, MPI_DOUBLE, 0, SOME_TAG,
	     comm_cart, MPI_STATUS_IGNORE);
  }

  on_particle_change();
#endif
}

/********************* REQ_SET_SOLV ********/
void mpi_send_solvation(int pnode, int part, double* solvation)
{
#ifdef SHANCHEN
  int ii;
  mpi_call(mpi_send_solvation_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    for(ii=0;ii<2*LB_COMPONENTS;ii++)
       p->p.solvation[ii]= solvation[ii];
  }
  else {
    MPI_Send(&solvation, LB_COMPONENTS, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
#endif
}

void mpi_send_solvation_slave(int pnode, int part)
{
#ifdef SHANCHEN
  if (pnode == this_node) {
    Particle *p = local_particles[part];
        MPI_Recv(&p->p.solvation, 2*LB_COMPONENTS, MPI_DOUBLE, 0, SOME_TAG,
	     comm_cart, MPI_STATUS_IGNORE);
  }

  on_particle_change();
#endif
}


/********************* REQ_SET_M ********/
void mpi_send_mass(int pnode, int part, double mass)
{
#ifdef MASS
  mpi_call(mpi_send_mass_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->p.mass = mass;
  }
  else {
    MPI_Send(&mass, 1, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
#endif
}

void mpi_send_mass_slave(int pnode, int part)
{
#ifdef MASS
  if (pnode == this_node) {
    Particle *p = local_particles[part];
        MPI_Recv(&p->p.mass, 1, MPI_DOUBLE, 0, SOME_TAG,
	     comm_cart, MPI_STATUS_IGNORE);
  }

  on_particle_change();
#endif
}

/********************* REQ_SET_RINERTIA ********/

void mpi_send_rotational_inertia(int pnode, int part, double rinertia[3])
{
#ifdef ROTATIONAL_INERTIA
  mpi_call(mpi_send_rotational_inertia_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->p.rinertia[0] = rinertia[0];
    p->p.rinertia[1] = rinertia[1];
    p->p.rinertia[2] = rinertia[2];
  }
  else {
    MPI_Send(rinertia, 3, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
#endif
}

void mpi_send_rotational_inertia_slave(int pnode, int part)
{
#ifdef ROTATIONAL_INERTIA
  if (pnode == this_node) {
    Particle *p = local_particles[part];
        MPI_Recv(p->p.rinertia, 3, MPI_DOUBLE, 0, SOME_TAG,
	     comm_cart, MPI_STATUS_IGNORE);
  }

  on_particle_change();
#endif
}

/********************* REQ_SET_TYPE ********/
void mpi_send_type(int pnode, int part, int type)
{
  mpi_call(mpi_send_type_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->p.type = type;
  }
  else
    MPI_Send(&type, 1, MPI_INT, pnode, SOME_TAG, comm_cart);

  on_particle_change();
}

void mpi_send_type_slave(int pnode, int part)
{
  if (pnode == this_node) {
    Particle *p = local_particles[part];
        MPI_Recv(&p->p.type, 1, MPI_INT, 0, SOME_TAG,
	     comm_cart, MPI_STATUS_IGNORE);
  }

  on_particle_change();
}

/********************* REQ_SET_MOLID ********/
void mpi_send_mol_id(int pnode, int part, int mid)
{
  mpi_call(mpi_send_mol_id_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->p.mol_id = mid;
  }
  else
    MPI_Send(&mid, 1, MPI_INT, pnode, SOME_TAG, comm_cart);

  on_particle_change();
}

void mpi_send_mol_id_slave(int pnode, int part)
{
  if (pnode == this_node) {
    Particle *p = local_particles[part];
        MPI_Recv(&p->p.mol_id, 1, MPI_INT, 0, SOME_TAG,
	     comm_cart, MPI_STATUS_IGNORE);
  }

  on_particle_change();
}

/********************* REQ_SET_QUAT ********/

void mpi_send_quat(int pnode, int part, double quat[4])
{
#ifdef ROTATION
  mpi_call(mpi_send_quat_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->r.quat[0] = quat[0];
    p->r.quat[1] = quat[1];
    p->r.quat[2] = quat[2];
    p->r.quat[3] = quat[3];
    convert_quat_to_quatu(p->r.quat, p->r.quatu);
#ifdef DIPOLES
    convert_quatu_to_dip(p->r.quatu, p->p.dipm, p->r.dip);
#endif
  }
  else {
    MPI_Send(quat, 4, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
#endif
}

void mpi_send_quat_slave(int pnode, int part)
{
#ifdef ROTATION
  if (pnode == this_node) {
    Particle *p = local_particles[part];
        MPI_Recv(p->r.quat, 4, MPI_DOUBLE, 0, SOME_TAG,
	     comm_cart, MPI_STATUS_IGNORE);
    convert_quat_to_quatu(p->r.quat, p->r.quatu);
#ifdef DIPOLES
    convert_quatu_to_dip(p->r.quatu, p->p.dipm, p->r.dip);
#endif
  }

  on_particle_change();
#endif
}

/********************* REQ_SET_OMEGA ********/

void mpi_send_omega(int pnode, int part, double omega[3])
{
#ifdef ROTATION
  mpi_call(mpi_send_omega_slave, pnode, part);

  if (pnode == this_node) {
   Particle *p = local_particles[part];
/*  memcpy(p->omega, omega, 3*sizeof(double));*/
    p->m.omega[0] = omega[0];
    p->m.omega[1] = omega[1];
    p->m.omega[2] = omega[2];
  }
  else {
    MPI_Send(omega, 3, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
#endif
}

void mpi_send_omega_slave(int pnode, int part)
{
#ifdef ROTATION
  if (pnode == this_node) {
    Particle *p = local_particles[part];
        MPI_Recv(p->m.omega, 3, MPI_DOUBLE, 0, SOME_TAG,
	     comm_cart, MPI_STATUS_IGNORE);
  }

  on_particle_change();
#endif
}

/********************* REQ_SET_TORQUE ********/

void mpi_send_torque(int pnode, int part, double torque[3])
{
#ifdef ROTATION
  mpi_call(mpi_send_torque_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->f.torque[0] = torque[0];
    p->f.torque[1] = torque[1];
    p->f.torque[2] = torque[2];
  }
  else {
    MPI_Send(torque, 3, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
#endif
}

void mpi_send_torque_slave(int pnode, int part)
{
#ifdef ROTATION
  if (pnode == this_node) {
    Particle *p = local_particles[part];
        MPI_Recv(p->f.torque, 3, MPI_DOUBLE, 0, SOME_TAG,
	     comm_cart, MPI_STATUS_IGNORE);
  }

  on_particle_change();
#endif
}

/********************* REQ_SET_DIP ********/

void mpi_send_dip(int pnode, int part, double dip[3])
{
#ifdef DIPOLES
  mpi_call(mpi_send_dip_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->r.dip[0] = dip[0];
    p->r.dip[1] = dip[1];
    p->r.dip[2] = dip[2];
#ifdef ROTATION
    convert_dip_to_quat(p->r.dip, p->r.quat, &p->p.dipm);
    convert_quat_to_quatu(p->r.quat, p->r.quatu);
#else
    p->p.dipm = sqrt(p->r.dip[0]*p->r.dip[0] + p->r.dip[1]*p->r.dip[1] + p->r.dip[2]*p->r.dip[2]);
#endif
  }
  else {
    MPI_Send(dip, 3, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
#endif
}

void mpi_send_dip_slave(int pnode, int part)
{
#ifdef DIPOLES
  if (pnode == this_node) {
    Particle *p = local_particles[part];
        MPI_Recv(p->r.dip, 3, MPI_DOUBLE, 0, SOME_TAG,
	     comm_cart, MPI_STATUS_IGNORE);
#ifdef ROTATION
    convert_dip_to_quat(p->r.dip, p->r.quat, &p->p.dipm);
    convert_quat_to_quatu(p->r.quat, p->r.quatu);
#else
    p->p.dipm = sqrt(p->r.dip[0]*p->r.dip[0] + p->r.dip[1]*p->r.dip[1] + p->r.dip[2]*p->r.dip[2]);
#endif
  }

  on_particle_change();
#endif
}

/********************* REQ_SET_DIPM ********/

void mpi_send_dipm(int pnode, int part, double dipm)
{
#ifdef DIPOLES
  mpi_call(mpi_send_dipm_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->p.dipm = dipm;
#ifdef ROTATION
    convert_quatu_to_dip(p->r.quatu, p->p.dipm, p->r.dip);
#endif
  }
  else {
    MPI_Send(&dipm, 1, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
#endif
}

void mpi_send_dipm_slave(int pnode, int part)
{
#ifdef DIPOLES
  if (pnode == this_node) {
    Particle *p = local_particles[part];
        MPI_Recv(&p->p.dipm, 1, MPI_DOUBLE, 0, SOME_TAG,
	     comm_cart, MPI_STATUS_IGNORE);
#ifdef ROTATION
    convert_quatu_to_dip(p->r.quatu, p->p.dipm, p->r.dip);
#endif
  }

  on_particle_change();
#endif
}

/********************* REQ_SET_ISVI ********/

void mpi_send_virtual(int pnode, int part, int isVirtual)
{
#ifdef VIRTUAL_SITES
  mpi_call(mpi_send_virtual_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->p.isVirtual = isVirtual;
  }
  else {
    MPI_Send(&isVirtual, 1, MPI_INT, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
#endif
}

void mpi_send_virtual_slave(int pnode, int part)
{
#ifdef VIRTUAL_SITES
  if (pnode == this_node) {
    Particle *p = local_particles[part];
        MPI_Recv(&p->p.isVirtual, 1, MPI_INT, 0, SOME_TAG,
	     comm_cart, MPI_STATUS_IGNORE);
  }

  on_particle_change();
#endif
}

/********************* REQ_SET_BOND ********/

void mpi_send_vs_relative(int pnode, int part, int vs_relative_to, double vs_distance)
{
#ifdef VIRTUAL_SITES_RELATIVE
  mpi_call(mpi_send_vs_relative_slave, pnode, part);

  // If the particle is on the node on which this function was called
  // set the values locally
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->p.vs_relative_to_particle_id = vs_relative_to;
    p->p.vs_relative_distance = vs_distance;
  }
  else {
    MPI_Send(&vs_relative_to, 1, MPI_INT, pnode, SOME_TAG, comm_cart);
    MPI_Send(&vs_distance, 1, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
#endif
}

void mpi_send_vs_relative_slave(int pnode, int part)
{
#ifdef VIRTUAL_SITES_RELATIVE
  if (pnode == this_node) {
    Particle *p = local_particles[part];
        MPI_Recv(&p->p.vs_relative_to_particle_id, 1, MPI_INT, 0, SOME_TAG,
	     comm_cart, MPI_STATUS_IGNORE);
    MPI_Recv(&p->p.vs_relative_distance, 1, MPI_DOUBLE, 0, SOME_TAG,
	     comm_cart, MPI_STATUS_IGNORE);
  }

  on_particle_change();
#endif
}

// ********************************

void mpi_send_rotation(int pnode, int part, int rot)
{
#ifdef ROTATION_PER_PARTICLE
  mpi_call(mpi_send_rotation_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->p.rotation = rot;
  }
  else {
    MPI_Send(&rot, 1, MPI_INT, pnode, SOME_TAG, MPI_COMM_WORLD);
  }

  on_particle_change();
#endif
}

void mpi_send_rotation_slave(int pnode, int part)
{
#ifdef ROTATION_PER_PARTICLE
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Status status;
    MPI_Recv(&p->p.rotation, 1, MPI_INT, 0, SOME_TAG,
	     MPI_COMM_WORLD, &status);
  }

  on_particle_change();
#endif
}

void mpi_observable_lb_radial_velocity_profile() 
{
#ifdef LB
  mpi_call(mpi_observable_lb_radial_velocity_profile_slave, 0, 0);
#endif
}
void mpi_observable_lb_radial_velocity_profile_slave(int pnode, int part)
{
#ifdef LB
  mpi_observable_lb_radial_velocity_profile_slave_implementation(); 
#endif
}


/********************* REQ_SET_BOND ********/


int mpi_send_bond(int pnode, int part, int *bond, int _delete)
{
  int bond_size, stat=0;
  
  mpi_call(mpi_send_bond_slave, pnode, part);

  bond_size = (bond) ? bonded_ia_params[bond[0]].num + 1 : 0;

  if (pnode == this_node) {
    stat = local_change_bond(part, bond, _delete);
    on_particle_change();
    return stat;
  }
  /* else */
  MPI_Send(&bond_size, 1, MPI_INT, pnode, SOME_TAG, comm_cart);
  if (bond_size)
    MPI_Send(bond, bond_size, MPI_INT, pnode, SOME_TAG, comm_cart);
  MPI_Send(&_delete, 1, MPI_INT, pnode, SOME_TAG, comm_cart);
  MPI_Recv(&stat, 1, MPI_INT, pnode, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
  on_particle_change();
  return stat;
}

void mpi_send_bond_slave(int pnode, int part)
{
  int bond_size=0, *bond, _delete=0, stat;
  
  if (pnode == this_node) {
    MPI_Recv(&bond_size, 1, MPI_INT, 0, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
    if (bond_size) {
      bond = (int *)malloc(bond_size*sizeof(int));
      MPI_Recv(bond, bond_size, MPI_INT, 0, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
    }
    else
      bond = NULL;
    MPI_Recv(&_delete, 1, MPI_INT, 0, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
    stat = local_change_bond(part, bond, _delete);
    if (bond)
      free(bond);
    MPI_Send(&stat, 1, MPI_INT, 0, SOME_TAG, comm_cart);
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
  
  /* fetch fixed data */
  if (pnode == this_node)
    memcpy(pdata, local_particles[part], sizeof(Particle));
  else {
    mpi_call(mpi_recv_part_slave, pnode, part);
    MPI_Recv(pdata, sizeof(Particle), MPI_BYTE, pnode,
	     SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
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
               SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
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
               SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
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

  MPI_Send(p, sizeof(Particle), MPI_BYTE, 0, SOME_TAG,
	   comm_cart);
  if (p->bl.n > 0)
    MPI_Send(p->bl.e, p->bl.n, MPI_INT, 0, SOME_TAG,
	     comm_cart);
#ifdef EXCLUSIONS
  if (p->el.n > 0)
    MPI_Send(p->el.e, p->el.n, MPI_INT, 0, SOME_TAG,
	     comm_cart);
#endif
}

/****************** REQ_REM_PART ************/
void mpi_remove_particle(int pnode, int part)
{
  mpi_call(mpi_remove_particle_slave, pnode, part);
  mpi_remove_particle_slave(pnode, part);
}

void mpi_remove_particle_slave(int pnode, int part)
{
  if (part != -1) {
    n_part--;

    if (pnode == this_node)
      local_remove_particle(part);
    
    remove_all_bonds_to(part);
  }
  else
    local_remove_all_particles();

  on_particle_change();
}

/********************* REQ_INTEGRATE ********/
int mpi_integrate(int n_steps, int reuse_forces)
{
  if (!correlations_autoupdate) {
    mpi_call(mpi_integrate_slave, n_steps, reuse_forces);
    integrate_vv(n_steps, reuse_forces);
    COMM_TRACE(fprintf(stderr, "%d: integration task %d done.\n", \
                       this_node, n_steps));
    return check_runtime_errors();
  } else {
    for (int i=0; i<n_steps; i++) {
      mpi_call(mpi_integrate_slave, 1, reuse_forces);
      integrate_vv(1, reuse_forces);
      reuse_forces = 0; // makes even less sense after the first time step
      COMM_TRACE(fprintf(stderr, "%d: integration task %d done.\n",     \
                         this_node, i));
      if (check_runtime_errors())
        return check_runtime_errors();
      autoupdate_correlations();
    }
  }

  return 0;
}

void mpi_integrate_slave(int n_steps, int reuse_forces)
{
  integrate_vv(n_steps, reuse_forces);
  COMM_TRACE(fprintf(stderr, "%d: integration for %d n_steps with %d reuse_forces done.\n", this_node, n_steps, reuse_forces));

  check_runtime_errors();
}

/*************** REQ_BCAST_IA ************/
void mpi_bcast_ia_params(int i, int j)
{
  mpi_call(mpi_bcast_ia_params_slave, i, j);
#ifdef TABULATED
  int tablesize = tabulated_forces.max;
#endif
  if (j>=0) {
    /* non-bonded interaction parameters */
    /* INCOMPATIBLE WHEN NODES USE DIFFERENT ARCHITECTURES */
    MPI_Bcast(get_ia_param(i, j), sizeof(IA_parameters), MPI_BYTE,
	      0, comm_cart);

    copy_ia_params(get_ia_param(j, i), get_ia_param(i, j));

#ifdef TABULATED
    /* If there are tabulated forces broadcast those as well */
    if ( get_ia_param(i,j)->TAB_maxval > 0) {
      /* First let all nodes know the new size for force and energy tables */
      MPI_Bcast(&tablesize, 1, MPI_INT, 0, comm_cart);

      /* Communicate the data */
      MPI_Bcast(tabulated_forces.e,tablesize, MPI_DOUBLE, 0 , comm_cart);
      MPI_Bcast(tabulated_energies.e,tablesize, MPI_DOUBLE, 0 , comm_cart);
    }    
#endif
  }
  else {
    /* bonded interaction parameters */
    /* INCOMPATIBLE WHEN NODES USE DIFFERENT ARCHITECTURES */
    MPI_Bcast(&(bonded_ia_params[i]),sizeof(Bonded_ia_parameters), MPI_BYTE,
	      0, comm_cart);
#ifdef TABULATED
    /* For tabulated potentials we have to send the tables extra */
    if(bonded_ia_params[i].type == BONDED_IA_TABULATED) {
      int size = bonded_ia_params[i].p.tab.npoints;
      MPI_Bcast(bonded_ia_params[i].p.tab.f, size, MPI_DOUBLE, 0 , comm_cart);
      MPI_Bcast(bonded_ia_params[i].p.tab.e, size, MPI_DOUBLE, 0 , comm_cart);
    }
#endif
#ifdef OVERLAPPED
    /* For overlapped potentials we have to send the ovelapped-functions extra */
    if(bonded_ia_params[i].type == BONDED_IA_OVERLAPPED) {
      int size = bonded_ia_params[i].p.overlap.noverlaps;
      MPI_Bcast(bonded_ia_params[i].p.overlap.para_a, size, MPI_DOUBLE, 0 , comm_cart);
      MPI_Bcast(bonded_ia_params[i].p.overlap.para_b, size, MPI_DOUBLE, 0 , comm_cart);
      MPI_Bcast(bonded_ia_params[i].p.overlap.para_c, size, MPI_DOUBLE, 0 , comm_cart);
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
	      0, comm_cart);

    copy_ia_params(get_ia_param(j, i), get_ia_param(i, j));

#ifdef TABULATED
    {
      int tablesize=0;
      /* If there are tabulated forces broadcast those as well */
      if ( get_ia_param(i,j)->TAB_maxval > 0) {
	/* Determine the new size for force and energy tables */
	MPI_Bcast(&tablesize,1,MPI_INT,0,comm_cart);
	/* Allocate sizes accordingly */
	realloc_doublelist(&tabulated_forces, tablesize);
	realloc_doublelist(&tabulated_energies, tablesize);
	/* Now communicate the data */
	MPI_Bcast(tabulated_forces.e,tablesize, MPI_DOUBLE, 0 , comm_cart);
	MPI_Bcast(tabulated_energies.e,tablesize, MPI_DOUBLE, 0 , comm_cart);
      }
    }
#endif
  } else { /* bonded interaction parameters */
    make_bond_type_exist(i); /* realloc bonded_ia_params on slave nodes! */
    MPI_Bcast(&(bonded_ia_params[i]),sizeof(Bonded_ia_parameters), MPI_BYTE,
	      0, comm_cart);
#ifdef TABULATED
    /* For tabulated potentials we have to send the tables extra */
    if(bonded_ia_params[i].type == BONDED_IA_TABULATED) {
      int size = bonded_ia_params[i].p.tab.npoints;
      /* alloc force and energy tables on slave nodes! */
      bonded_ia_params[i].p.tab.f = (double*)malloc(size*sizeof(double));
      bonded_ia_params[i].p.tab.e = (double*)malloc(size*sizeof(double));
      MPI_Bcast(bonded_ia_params[i].p.tab.f, size, MPI_DOUBLE, 0 , comm_cart);
      MPI_Bcast(bonded_ia_params[i].p.tab.e, size, MPI_DOUBLE, 0 , comm_cart);
    }
#endif
#ifdef OVERLAPPED
    /* For overlapped potentials we have to send the ovelapped-functions extra */
    if(bonded_ia_params[i].type == BONDED_IA_OVERLAPPED) {
      int size = bonded_ia_params[i].p.overlap.noverlaps;
      /* alloc overlapped parameter arrays on slave nodes! */
      bonded_ia_params[i].p.overlap.para_a = (double*)malloc(size*sizeof(double));
      bonded_ia_params[i].p.overlap.para_b = (double*)malloc(size*sizeof(double));
      bonded_ia_params[i].p.overlap.para_c = (double*)malloc(size*sizeof(double));
      MPI_Bcast(bonded_ia_params[i].p.overlap.para_a, size, MPI_DOUBLE, 0 , comm_cart);
      MPI_Bcast(bonded_ia_params[i].p.overlap.para_b, size, MPI_DOUBLE, 0 , comm_cart);
      MPI_Bcast(bonded_ia_params[i].p.overlap.para_c, size, MPI_DOUBLE, 0 , comm_cart);
    }
#endif
  }

  on_short_range_ia_change();
}

/*************** REQ_BCAST_IA_SIZE ************/

void mpi_bcast_n_particle_types(int ns)
{
  mpi_call(mpi_bcast_n_particle_types_slave, -1, ns);
  mpi_bcast_n_particle_types_slave(-1, ns);

}

void mpi_bcast_n_particle_types_slave(int pnode, int ns)
{
  realloc_ia_params(ns);
}

/*************** REQ_GATHER ************/
void mpi_gather_stats(int job, void *result, void *result_t, void *result_nb, void *result_t_nb)
{
  switch (job) {
  case 1:
    mpi_call(mpi_gather_stats_slave, -1, 1);
    energy_calc((double*)result);
    break;
  case 2:
    /* calculate and reduce (sum up) virials for 'analyze pressure' or
       'analyze stress_tensor' */
    mpi_call(mpi_gather_stats_slave, -1, 2);
    pressure_calc((double*)result, (double*)result_t,
                  (double*)result_nb, (double*)result_t_nb,0);
    break;
  case 3:
    mpi_call(mpi_gather_stats_slave, -1, 3);
    pressure_calc((double*)result, (double*)result_t, 
                  (double*)result_nb, (double*)result_t_nb,1);
    break;
  case 4:
    mpi_call(mpi_gather_stats_slave, -1, 4);
    predict_momentum_particles((double*)result);
    break;
#ifdef LB
  case 5:
    mpi_call(mpi_gather_stats_slave, -1, 5);
    lb_calc_fluid_mass((double*)result);
    break;
  case 6: 
    mpi_call(mpi_gather_stats_slave, -1, 6);
    lb_calc_fluid_momentum((double*)result);
    break;
  case 7:
    mpi_call(mpi_gather_stats_slave, -1, 7);
    lb_calc_fluid_temp((double*)result);
    break;
#ifdef LB_BOUNDARIES
  case 8:
    mpi_call(mpi_gather_stats_slave, -1, 8);
    lb_collect_boundary_forces((double*)result);
    break;
#endif
#endif
  default:
    fprintf(stderr, "%d: INTERNAL ERROR: illegal request %d for mpi_gather_stats_slave\n", this_node, job);
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
    /* calculate and reduce (sum up) virials for 'analyze pressure' or 'analyze stress_tensor'*/
    pressure_calc(NULL,NULL,NULL,NULL,0);
    break;
  case 3:
    /* calculate and reduce (sum up) virials, revert velocities half a timestep for 'analyze p_inst' */
    pressure_calc(NULL,NULL,NULL,NULL,1);
    break;
  case 4:
    predict_momentum_particles(NULL);
    break;
#ifdef LB
  case 5:
    lb_calc_fluid_mass(NULL);
    break;
  case 6:
    lb_calc_fluid_momentum(NULL);
    break;
  case 7:
    lb_calc_fluid_temp(NULL);
    break;
#ifdef LB_BOUNDARIES
  case 8:
    lb_collect_boundary_forces(NULL);
    break;
#endif
#endif
  default:
    fprintf(stderr, "%d: INTERNAL ERROR: illegal request %d for mpi_gather_stats_slave\n", this_node, job);
    errexit();
  }
}

/*************** REQ_GET_LOCAL_STRESS_TENSOR ************/
void mpi_local_stress_tensor(DoubleList *TensorInBin, int bins[3], int periodic[3], double range_start[3], double range[3]) {
  
  int i,j;
  DoubleList *TensorInBin_;
  PTENSOR_TRACE(fprintf(stderr,"%d: mpi_local_stress_tensor: Broadcasting local_stress_tensor parameters\n",this_node));


  mpi_call(mpi_local_stress_tensor_slave,-1,0);

  TensorInBin_ = (DoubleList*)malloc(bins[0]*bins[1]*bins[2]*sizeof(DoubleList));
  for ( i = 0 ; i < bins[0]*bins[1]*bins[2]; i++ ) {
    init_doublelist(&TensorInBin_[i]);
    alloc_doublelist(&TensorInBin_[i],9);
    for ( j = 0 ; j < 9 ; j++ ) {
      TensorInBin_[i].e[j] = TensorInBin[i].e[j];
      TensorInBin[i].e[j] = 0;
    }
  }
  
  MPI_Bcast(bins, 3, MPI_INT, 0, comm_cart);
  MPI_Bcast(periodic, 3, MPI_INT, 0, comm_cart);
  MPI_Bcast(range_start, 3, MPI_DOUBLE, 0, comm_cart);
  MPI_Bcast(range, 3, MPI_DOUBLE, 0, comm_cart);
  
  PTENSOR_TRACE(fprintf(stderr,"%d: mpi_local_stress_tensor: Call local_stress_tensor_calc\n",this_node));
  local_stress_tensor_calc(TensorInBin_, bins, periodic, range_start, range);
  
  PTENSOR_TRACE(fprintf(stderr,"%d: mpi_local_stress_tensor: Reduce local stress tensors with MPI_Reduce\n",this_node));
  for (i=0;i<bins[0]*bins[1]*bins[2];i++) {
    MPI_Reduce(TensorInBin_[i].e, TensorInBin[i].e, 9, MPI_DOUBLE, MPI_SUM, 0, comm_cart); 
  }
  
}

void mpi_local_stress_tensor_slave(int ana_num, int job) {
  int bins[3] = {0,0,0};
  int periodic[3]= {0,0,0};
  double range_start[3]= {0,0,0};
  double range[3]= {0,0,0};
  DoubleList *TensorInBin;
  int i, j;

  MPI_Bcast(bins, 3, MPI_INT, 0, comm_cart);
  MPI_Bcast(periodic, 3, MPI_INT, 0, comm_cart);
  MPI_Bcast(range_start, 3, MPI_DOUBLE, 0, comm_cart);
  MPI_Bcast(range, 3, MPI_DOUBLE, 0, comm_cart);
  
  TensorInBin = (DoubleList*)malloc(bins[0]*bins[1]*bins[2]*sizeof(DoubleList));
  for ( i = 0 ; i < bins[0]*bins[1]*bins[2]; i++ ) {
    init_doublelist(&TensorInBin[i]);
    alloc_doublelist(&TensorInBin[i],9);
    for ( j = 0 ; j < 9 ; j++ ) {
      TensorInBin[i].e[j] = 0.0;
    }
  }
  
  local_stress_tensor_calc(TensorInBin, bins, periodic, range_start, range);

  for (i=0;i<bins[0]*bins[1]*bins[2];i++) {
    MPI_Reduce(TensorInBin[i].e, NULL, 9, MPI_DOUBLE, MPI_SUM, 0, comm_cart);
    PTENSOR_TRACE(fprintf(stderr,"%d: mpi_local_stress_tensor: Tensor sent in bin %d is {",this_node,i));
    for (j=0; j<9;j++) {
      PTENSOR_TRACE(fprintf(stderr,"%f ",TensorInBin[i].e[j]));
    }
    PTENSOR_TRACE(fprintf(stderr,"}\n"));
  }

  for ( i = 0 ; i < bins[0]*bins[1]*bins[2] ; i++ ) {
    realloc_doublelist(&TensorInBin[i],0);
  }
  free(TensorInBin);
}

/*************** REQ_GETPARTS ************/
void mpi_get_particles(Particle *result, IntList *bi)
{
  IntList local_bi;
  int local_part;
  int tot_size, i, g, pnode;
  int *sizes;
  Cell *cell;
  int c;
  
  mpi_call(mpi_get_particles_slave, -1, bi != NULL);

  sizes = (int*)malloc(sizeof(int)*n_nodes);
  local_part = cells_get_n_particles();

  /* first collect number of particles on each node */
  MPI_Gather(&local_part, 1, MPI_INT, sizes, 1, MPI_INT, 0, comm_cart);
  tot_size = 0;
  for (i = 0; i < n_nodes; i++)
    tot_size += sizes[i];

  if (tot_size!=n_part) {
    fprintf(stderr,"%d: ERROR: mpi_get_particles: n_part %d, but I counted %d. Exiting...\n",
	    this_node, n_part, tot_size);
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
	MPI_Recv(&result[g], sizes[pnode]*sizeof(Particle), MPI_BYTE, pnode, SOME_TAG,
		 comm_cart, MPI_STATUS_IGNORE);
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

#ifdef DIPOLES
  COMM_TRACE(for(i = 0; i < tot_size; i++) {
    printf("%d: %d -> %d %d  (%f, %f, %f) (%f, %f, %f)\n", this_node, i, result[i].p.identity, result[i].p.type,
	   result[i].r.p[0], result[i].r.p[1], result[i].r.p[2], result[i].r.dip[0],
	   result[i].r.dip[1], result[i].r.dip[2]);
  });
#endif

  /* gather bonding information */
  if (bi) {
    int *bonds;

    init_intlist(bi);
    MPI_Gather(&local_bi.n, 1, MPI_INT, sizes, 1, MPI_INT, 0, comm_cart);
    for (pnode = 0; pnode < n_nodes; pnode++) {
      if (sizes[pnode] > 0) {
	realloc_intlist(bi, bi->n + sizes[pnode]);

	if (pnode == this_node)
	  memcpy(&bi->e[bi->n], local_bi.e, sizes[pnode]*sizeof(int));
	else
	  MPI_Recv(&bi->e[bi->n], sizes[pnode], MPI_INT, pnode, SOME_TAG,
		   comm_cart, MPI_STATUS_IGNORE);

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
	     0, comm_cart);

  if (n_part > 0) {
    IntList local_bi;

    /* get (unsorted) particle informations as an array of type 'particle' */
    /* then get the particle information */
    result = (Particle*)malloc(n_part*sizeof(Particle));

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
    MPI_Send(result, n_part*sizeof(Particle), MPI_BYTE, 0, SOME_TAG, comm_cart);
    free(result);

    if (bi) {
      COMM_TRACE(fprintf(stderr, "%d: sending bonds\n", this_node));

      MPI_Gather(&local_bi.n, 1, MPI_INT, NULL, 1, MPI_INT, 0, comm_cart);      
      if (local_bi.n > 0)
	MPI_Send(local_bi.e, local_bi.n, MPI_INT, 0, SOME_TAG, comm_cart);
      realloc_intlist(&local_bi, 0);
    }
  }
  else {
    if (bi) {
      /* inform master node that we do not have bonds (as we don't have particles) */
      g = 0;
      MPI_Gather(&g, 1, MPI_INT, NULL, 1, MPI_INT, 0, comm_cart);      
    }
  }
}

/*************** REQ_SET_TIME_STEP ************/
void mpi_set_time_step(double time_s)
{
  double old_ts = time_step;

  mpi_call(mpi_set_time_step_slave, -1, 0);

  time_step = time_s;

  time_step_squared=time_step * time_step;
  time_step_squared_half = time_step_squared /2.;
  time_step_half= time_step / 2.;

  MPI_Bcast(&time_step, 1, MPI_DOUBLE, 0, comm_cart);

  rescale_velocities(time_step / old_ts);
  on_parameter_change(FIELD_TIMESTEP);
}

void mpi_set_time_step_slave(int node, int i)
{
  double old_ts = time_step;

  MPI_Bcast(&time_step, 1, MPI_DOUBLE, 0, comm_cart);
  rescale_velocities(time_step / old_ts);
  on_parameter_change(FIELD_TIMESTEP);
  time_step_squared=time_step * time_step;
  time_step_squared_half = time_step_squared /2.;
  time_step_half= time_step / 2.;
}

/*************** REQ_BCAST_COULOMB ************/
void mpi_bcast_coulomb_params()
{
#if  defined(ELECTROSTATICS) || defined(DIPOLES)
  mpi_call(mpi_bcast_coulomb_params_slave, 1, 0);
  mpi_bcast_coulomb_params_slave(-1, 0);
#endif
}

void mpi_bcast_coulomb_params_slave(int node, int parm)
{   
#if defined(ELECTROSTATICS) || defined(DIPOLES)
  MPI_Bcast(&coulomb, sizeof(Coulomb_parameters), MPI_BYTE, 0, comm_cart);

#ifdef ELECTROSTATICS
  switch (coulomb.method) {
  case COULOMB_NONE:
    break;
#ifdef P3M
  case COULOMB_ELC_P3M:
    MPI_Bcast(&elc_params, sizeof(ELC_struct), MPI_BYTE, 0, comm_cart);
    // fall through
  case COULOMB_P3M_GPU:
  case COULOMB_P3M:
    MPI_Bcast(&p3m.params, sizeof(p3m_parameter_struct), MPI_BYTE, 0, comm_cart);
    break;
#endif
  case COULOMB_DH:
  case COULOMB_DH_PW:
    MPI_Bcast(&dh_params, sizeof(Debye_hueckel_params), MPI_BYTE, 0, comm_cart);
    break;
  case COULOMB_MMM1D:
    MPI_Bcast(&mmm1d_params, sizeof(MMM1D_struct), MPI_BYTE, 0, comm_cart);
    break;
  case COULOMB_MMM2D:
    MPI_Bcast(&mmm2d_params, sizeof(MMM2D_struct), MPI_BYTE, 0, comm_cart);
    break;  
  case COULOMB_MAGGS:
    MPI_Bcast(&maggs, sizeof(MAGGS_struct), MPI_BYTE, 0, comm_cart); 
    break;
  case COULOMB_RF:
  case COULOMB_INTER_RF:
    MPI_Bcast(&rf_params, sizeof(Reaction_field_params), MPI_BYTE, 0, comm_cart);
    break;
  default:
    fprintf(stderr, "%d: INTERNAL ERROR: cannot bcast coulomb params for unknown method %d\n", this_node, coulomb.method);
    errexit();
  }
#endif

#ifdef DIPOLES
  switch (coulomb.Dmethod) {
  case DIPOLAR_NONE:
    break;
#ifdef DP3M
  case DIPOLAR_MDLC_P3M:
    MPI_Bcast(&dlc_params, sizeof(DLC_struct), MPI_BYTE, 0, comm_cart);
    // fall through
  case DIPOLAR_P3M:
    MPI_Bcast(&dp3m.params, sizeof(p3m_parameter_struct), MPI_BYTE, 0, comm_cart);
    break;
#endif
  case DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA :
   break;
 case  DIPOLAR_MDLC_DS:
     //fall trough
 case  DIPOLAR_DS:
    break;   
  default:
    fprintf(stderr, "%d: INTERNAL ERROR: cannot bcast dipolar params for unknown method %d\n", this_node, coulomb.Dmethod);
    errexit();
  }

#endif  
  
  on_coulomb_change();
#endif
}

/*************** REQ_BCAST_COULOMB ************/
void mpi_bcast_collision_params()
{
#ifdef COLLISION_DETECTION
  mpi_call(mpi_bcast_collision_params_slave, 1, 0);
  mpi_bcast_collision_params_slave(-1, 0);
#endif
}

void mpi_bcast_collision_params_slave(int node, int parm)
{   
#ifdef COLLISION_DETECTION
  MPI_Bcast(&collision_params, sizeof(Collision_parameters), MPI_BYTE, 0, comm_cart);

  recalc_forces = 1;
#endif
}

/****************** REQ_SET_EXT ************/

void mpi_send_ext_torque(int pnode, int part, int flag, int mask, double torque[3])
{
#ifdef EXTERNAL_FORCES
  #ifdef ROTATION
    int s_buf[2];
    mpi_call(mpi_send_ext_torque_slave, pnode, part);

    if (pnode == this_node) {
      Particle *p = local_particles[part];
      /* mask out old flags */
      p->l.ext_flag &= ~mask;
      /* set new values */
      p->l.ext_flag |= flag;

      if (mask & PARTICLE_EXT_TORQUE) 
        memcpy(p->l.ext_torque, torque, 3*sizeof(double));
    }
    else {
      s_buf[0] = flag; s_buf[1] = mask;
      MPI_Send(s_buf, 2, MPI_INT, pnode, SOME_TAG, comm_cart);
      if (mask & PARTICLE_EXT_TORQUE)
        MPI_Send(torque, 3, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
    }

    on_particle_change();
  #endif
#endif
}

void mpi_send_ext_torque_slave(int pnode, int part)
{
#ifdef EXTERNAL_FORCES
  #ifdef ROTATION
    int s_buf[2]={0,0};
    if (pnode == this_node) {
      Particle *p = local_particles[part];
          MPI_Recv(s_buf, 2, MPI_INT, 0, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
      /* mask out old flags */
      p->l.ext_flag &= ~s_buf[1];
      /* set new values */
      p->l.ext_flag |= s_buf[0];
      
      if (s_buf[1] & PARTICLE_EXT_TORQUE)
        MPI_Recv(p->l.ext_torque, 3, MPI_DOUBLE, 0, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
    }

    on_particle_change();
  #endif
#endif
}

void mpi_send_ext_force(int pnode, int part, int flag, int mask, double force[3])
{
#ifdef EXTERNAL_FORCES
  int s_buf[2];
  mpi_call(mpi_send_ext_force_slave, pnode, part);

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
    MPI_Send(s_buf, 2, MPI_INT, pnode, SOME_TAG, comm_cart);
    if (mask & PARTICLE_EXT_FORCE)
      MPI_Send(force, 3, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
#endif
}

void mpi_send_ext_force_slave(int pnode, int part)
{
#ifdef EXTERNAL_FORCES
  int s_buf[2]={0,0};
  if (pnode == this_node) {
    Particle *p = local_particles[part];
        MPI_Recv(s_buf, 2, MPI_INT, 0, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
    /* mask out old flags */
    p->l.ext_flag &= ~s_buf[1];
    /* set new values */
    p->l.ext_flag |= s_buf[0];
    
    if (s_buf[1] & PARTICLE_EXT_FORCE)
      MPI_Recv(p->l.ext_force, 3, MPI_DOUBLE, 0, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
  }

  on_particle_change();
#endif
}

/*************** REQ_BCAST_CONSTR ************/
void mpi_bcast_constraint(int del_num)
{
#ifdef CONSTRAINTS
  mpi_call(mpi_bcast_constraint_slave, 0, del_num);

  if (del_num == -1) {
    /* bcast new constraint */
    MPI_Bcast(&constraints[n_constraints-1], sizeof(Constraint), MPI_BYTE, 0, comm_cart);
  }
  else if (del_num == -2) {
    /* delete all constraints */
    n_constraints = 0;
    constraints = (Constraint*)realloc(constraints,n_constraints*sizeof(Constraint));
  }
  else {
    memcpy(&constraints[del_num],&constraints[n_constraints-1],sizeof(Constraint));
    n_constraints--;
    constraints = (Constraint*)realloc(constraints,n_constraints*sizeof(Constraint));
  }

  on_constraint_change();
#endif
}

void mpi_bcast_constraint_slave(int node, int parm)
{   
#ifdef CONSTRAINTS
  if(parm == -1) {
    n_constraints++;
    constraints = (Constraint*)realloc(constraints,n_constraints*sizeof(Constraint));
    MPI_Bcast(&constraints[n_constraints-1], sizeof(Constraint), MPI_BYTE, 0, comm_cart);
  }
  else if (parm == -2) {
    /* delete all constraints */
    n_constraints = 0;
    constraints = (Constraint*)realloc(constraints,n_constraints*sizeof(Constraint));
  }
  else {
    memcpy(&constraints[parm],&constraints[n_constraints-1],sizeof(Constraint));
    n_constraints--;
    constraints = (Constraint*)realloc(constraints,n_constraints*sizeof(Constraint));    
  }

  on_constraint_change();
#endif
}

/*************** REQ_BCAST_LBBOUNDARY ************/
void mpi_bcast_lbboundary(int del_num)
{
#if defined(LB_BOUNDARIES) || defined(LB_BOUNDARIES_GPU)
  mpi_call(mpi_bcast_lbboundary_slave, 0, del_num);
  
  if (del_num == -1) {
    /* bcast new boundaries */
    MPI_Bcast(&lb_boundaries[n_lb_boundaries-1], sizeof(LB_Boundary), MPI_BYTE, 0, comm_cart);
  }
  else if (del_num == -2) {
    /* delete all boundaries */
    n_lb_boundaries = 0;
    lb_boundaries = (LB_Boundary*) realloc(lb_boundaries,n_lb_boundaries*sizeof(LB_Boundary));
  }
#if defined(LB_BOUNDARIES_GPU)
  else if (del_num == -3) {
  //nothing, GPU code just requires to call on_lbboundary_change()
  }
#endif
  else {
    memcpy(&lb_boundaries[del_num],&lb_boundaries[n_lb_boundaries-1],sizeof(LB_Boundary));
    n_lb_boundaries--;
    lb_boundaries = (LB_Boundary*) realloc(lb_boundaries,n_lb_boundaries*sizeof(LB_Boundary));
  }

  on_lbboundary_change();
#endif
}

void mpi_bcast_lbboundary_slave(int node, int parm)
{   
#if defined(LB_BOUNDARIES) || defined(LB_BOUNDARIES_GPU)

#if defined(LB_BOUNDARIES)
  if(parm == -1) {
    n_lb_boundaries++;
    lb_boundaries = (LB_Boundary*) realloc(lb_boundaries,n_lb_boundaries*sizeof(LB_Boundary));
    MPI_Bcast(&lb_boundaries[n_lb_boundaries-1], sizeof(LB_Boundary), MPI_BYTE, 0, comm_cart);
  }
  else if (parm == -2) {
    /* delete all boundaries */
    n_lb_boundaries = 0;
    lb_boundaries = (LB_Boundary*) realloc(lb_boundaries,n_lb_boundaries*sizeof(LB_Boundary));
  }
#if defined(LB_BOUNDARIES_GPU)
  else if (parm == -3) {
  //nothing, GPU code just requires to call on_lbboundary_change()
  }
#endif
  else {
    memcpy(&lb_boundaries[parm],&lb_boundaries[n_lb_boundaries-1],sizeof(LB_Boundary));
    n_lb_boundaries--;
    lb_boundaries = (LB_Boundary*) realloc(lb_boundaries,n_lb_boundaries*sizeof(LB_Boundary));    
  }
#endif

  on_lbboundary_change();
#endif
}

/*************** REQ_RANDOM_SEED ************/
void mpi_random_seed(int cnt, long *seed) {
  long this_idum = print_random_seed();

  mpi_call(mpi_random_seed_slave, -1, cnt);

  if (cnt==0) {
    RANDOM_TRACE(printf("%d: Have seed %ld\n",this_node,this_idum));
    MPI_Gather(&this_idum,1,MPI_LONG,seed,1,MPI_LONG,0,comm_cart); }
  else {
    MPI_Scatter(seed,1,MPI_LONG,&this_idum,1,MPI_LONG,0,comm_cart);
    RANDOM_TRACE(printf("%d: Received seed %ld\n",this_node,this_idum));
    init_random_seed(this_idum);
  }
}

void mpi_random_seed_slave(int pnode, int cnt) {
  long this_idum = print_random_seed();

  if (cnt==0) {
    RANDOM_TRACE(printf("%d: Have seed %ld\n",this_node,this_idum));
    MPI_Gather(&this_idum,1,MPI_LONG,NULL,0,MPI_LONG,0,comm_cart); }
  else {
    MPI_Scatter(NULL,1,MPI_LONG,&this_idum,1,MPI_LONG,0,comm_cart);
    RANDOM_TRACE(printf("%d: Received seed %ld\n",this_node,this_idum));
    init_random_seed(this_idum);
  }
}

/*************** REQ_RANDOM_STAT ************/
void mpi_random_stat(int cnt, RandomStatus *stat) {
  RandomStatus this_stat = print_random_stat();

  mpi_call(mpi_random_stat_slave, -1, cnt);

  if (cnt==0) {
    RANDOM_TRACE(printf("%d: Have status %ld/%ld/...\n",this_node,this_stat.idum,this_stat.iy));
    MPI_Gather(&this_stat,1*sizeof(RandomStatus),MPI_BYTE,stat,1*sizeof(RandomStatus),MPI_BYTE,0,comm_cart); }
  else {
    MPI_Scatter(stat,1*sizeof(RandomStatus),MPI_BYTE,&this_stat,1*sizeof(RandomStatus),MPI_BYTE,0,comm_cart);
    RANDOM_TRACE(printf("%d: Received status %ld/%ld/...\n",this_node,this_stat.idum,this_stat.iy));
    init_random_stat(this_stat);
  }
}

void mpi_random_stat_slave(int pnode, int cnt) {
  RandomStatus this_stat = print_random_stat();

  if (cnt==0) {
    RANDOM_TRACE(printf("%d: Have status %ld/%ld/...\n",this_node,this_stat.idum,this_stat.iy));
    MPI_Gather(&this_stat,1*sizeof(RandomStatus),MPI_BYTE,NULL,0,MPI_BYTE,0,comm_cart); }
  else {
    MPI_Scatter(NULL,0,MPI_BYTE,&this_stat,1*sizeof(RandomStatus),MPI_BYTE,0,comm_cart);
    RANDOM_TRACE(printf("%d: Received status %ld/%ld/...\n",this_node,this_stat.idum,this_stat.iy));
    init_random_stat(this_stat);
  }
}


/*************** REQ_BCAST_LJFORCECAP ************/
/*************** REQ_BCAST_LJANGLEFORCECAP ************/
/*************** REQ_BCAST_MORSEFORCECAP ************/
/*************** REQ_BCAST_BUCKFORCECAP ************/
/*************** REQ_BCAST_TABFORCECAP ************/
void mpi_cap_forces(double fc)
{
  force_cap = fc;
  mpi_call(mpi_cap_forces_slave, 1, 0);
  mpi_cap_forces_slave(1, 0);
}

void mpi_cap_forces_slave(int node, int parm)
{
#ifdef LENNARD_JONES
  MPI_Bcast(&force_cap, 1, MPI_DOUBLE, 0, comm_cart);
  calc_lj_cap_radii();
#ifdef LENNARD_JONES_GENERIC
  calc_ljgen_cap_radii();
#endif
#ifdef LJCOS2
  calc_ljcos2_cap_radii();
#endif
  on_short_range_ia_change();
#endif
#ifdef LJ_ANGLE
  MPI_Bcast(&force_cap, 1, MPI_DOUBLE, 0, comm_cart);
  calc_ljangle_cap_radii();
  on_short_range_ia_change();
#endif
#ifdef MORSE
  MPI_Bcast(&force_cap, 1, MPI_DOUBLE, 0, comm_cart);
  calc_morse_cap_radii();
  on_short_range_ia_change();
#endif
#ifdef BUCKINGHAM
  MPI_Bcast(&force_cap, 1, MPI_DOUBLE, 0, comm_cart);
  calc_buck_cap_radii();
  on_short_range_ia_change();
#endif
#ifdef TABULATED
  MPI_Bcast(&force_cap, 1, MPI_DOUBLE, 0, comm_cart);
/* to do: check if "check_tab_forcecap" is still useful since force capping
   is defined globally now -- the cap for other forces would be removed too! 
  check_tab_forcecap(force_cap);
*/
  on_short_range_ia_change();
#endif
}

/*************** REQ_GET_CONSFOR ************/
void mpi_get_constraint_force(int cons, double force[3])
{
#ifdef CONSTRAINTS
  mpi_call(mpi_get_constraint_force_slave, -1, cons);
  MPI_Reduce(constraints[cons].part_rep.f.f, force, 3, MPI_DOUBLE, MPI_SUM, 0, comm_cart);
#endif
}

void mpi_get_constraint_force_slave(int node, int parm)
{
#ifdef CONSTRAINTS
  MPI_Reduce(constraints[parm].part_rep.f.f, NULL, 3, MPI_DOUBLE, MPI_SUM, 0, comm_cart);
#endif
}

/*************** REQ_BIT_RANDOM_SEED ************/
void mpi_bit_random_seed(int cnt, int *seed) {
  int this_idum = print_bit_random_seed();

  mpi_call(mpi_bit_random_seed_slave, -1, cnt);

  if (cnt==0) {
    RANDOM_TRACE(printf("%d: Have seed %d\n",this_node,this_idum));
    MPI_Gather(&this_idum,1,MPI_INT,seed,1,MPI_INT,0,comm_cart); }
  else {
    MPI_Scatter(seed,1,MPI_INT,&this_idum,1,MPI_INT,0,comm_cart);
    RANDOM_TRACE(printf("%d: Received seed %d\n",this_node,this_idum));
    init_bit_random_generator(this_idum);
  }
}

void mpi_bit_random_seed_slave(int pnode, int cnt) {
  int this_idum = print_bit_random_seed();

  if (cnt==0) {
    RANDOM_TRACE(printf("%d: Have seed %d\n",this_node,this_idum));
    MPI_Gather(&this_idum,1,MPI_INT,NULL,0,MPI_INT,0,comm_cart); }
  else {
    MPI_Scatter(NULL,1,MPI_INT,&this_idum,1,MPI_INT,0,comm_cart);
    RANDOM_TRACE(printf("%d: Received seed %d\n",this_node,this_idum));
    init_bit_random_generator(this_idum);
  }
}

/*************** REQ_BIT_RANDOM_STAT ************/
void mpi_bit_random_stat(int cnt, BitRandomStatus *stat) {
  BitRandomStatus this_stat = print_bit_random_stat();

  mpi_call(mpi_bit_random_stat_slave, -1, cnt);

  if (cnt==0) {
    RANDOM_TRACE(printf("%d: Have status %d/%d/...\n",this_node,this_stat.random_pointer_1,this_stat.random_pointer_2));
    MPI_Gather(&this_stat,1*sizeof(BitRandomStatus),MPI_BYTE,stat,1*sizeof(BitRandomStatus),MPI_BYTE,0,comm_cart); }
  else {
    MPI_Scatter(stat,1*sizeof(BitRandomStatus),MPI_BYTE,&this_stat,1*sizeof(BitRandomStatus),MPI_BYTE,0,comm_cart);
    RANDOM_TRACE(printf("%d: Received status %d/%d/...\n",this_node,this_stat.random_pointer_1,this_stat.random_pointer_2));
    init_bit_random_stat(this_stat);
  }
}

void mpi_bit_random_stat_slave(int pnode, int cnt) {
  BitRandomStatus this_stat = print_bit_random_stat();

  if (cnt==0) {
    RANDOM_TRACE(printf("%d: Have status %d/%d/...\n",this_node,this_stat.random_pointer_1,this_stat.random_pointer_2));
    MPI_Gather(&this_stat,1*sizeof(BitRandomStatus),MPI_BYTE,NULL,0,MPI_BYTE,0,comm_cart); }
  else {
    MPI_Scatter(NULL,0,MPI_BYTE,&this_stat,1*sizeof(BitRandomStatus),MPI_BYTE,0,comm_cart);
    RANDOM_TRACE(printf("%d: Received status %d/%d/...\n",this_node,this_stat.random_pointer_1,this_stat.random_pointer_2));
    init_bit_random_stat(this_stat);
  }
}

/****************** REQ_RESCALE_PART ************/

void mpi_rescale_particles(int dir, double scale) {
  int pnode;

  mpi_call(mpi_rescale_particles_slave, -1, dir);
  for (pnode = 0; pnode < n_nodes; pnode++) {
    if (pnode == this_node) {
      local_rescale_particles(dir, scale); }
    else {
      MPI_Send(&scale, 1, MPI_DOUBLE, pnode, SOME_TAG, comm_cart); }
  }
  on_particle_change();
}

void mpi_rescale_particles_slave(int pnode, int dir) {
  double scale=0.0; 
  MPI_Recv(&scale, 1, MPI_DOUBLE, 0, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
  local_rescale_particles(dir, scale);
  on_particle_change();
}

/*************** REQ_BCAST_CS *****************/

void mpi_bcast_cell_structure(int cs)
{
  mpi_call(mpi_bcast_cell_structure_slave, -1, cs);
  cells_re_init(cs);
}

void mpi_bcast_cell_structure_slave(int pnode, int cs)
{
  cells_re_init(cs);
}

/*************** REQ_BCAST_NPTISO_GEOM *****************/

void mpi_bcast_nptiso_geom()
{
  mpi_call(mpi_bcast_nptiso_geom_slave, -1 , 0);
  mpi_bcast_nptiso_geom_slave(-1,0);

}

void mpi_bcast_nptiso_geom_slave(int node,int parm)
{
  MPI_Bcast(&nptiso.geometry, 1,MPI_INT, 0, comm_cart);
  MPI_Bcast(&nptiso.dimension, 1,MPI_INT, 0, comm_cart);
  MPI_Bcast(&nptiso.cubic_box, 1,MPI_INT, 0, comm_cart);
  MPI_Bcast(&nptiso.non_const_dim, 1,MPI_INT, 0, comm_cart);

}

/***************REQ_UPDATE_MOL_IDS *********************/

void mpi_update_mol_ids()
{
  mpi_call(mpi_update_mol_ids_slave, -1, 0);
  mpi_update_mol_ids_slave(-1, 0);
}

void mpi_update_mol_ids_slave(int node,int parm)
{
  update_mol_ids_setchains();
}

/******************* REQ_SYNC_TOPO ********************/
int mpi_sync_topo_part_info() {
  int i;
  int molsize=0;
  int moltype=0;
  int n_mols=0;
  
  mpi_call(mpi_sync_topo_part_info_slave,-1,0);
  n_mols = n_molecules;
  MPI_Bcast(&n_mols,1,MPI_INT,0,comm_cart);

  for ( i = 0 ; i < n_molecules ; i++) {
    molsize = topology[i].part.n;
    moltype = topology[i].type;

#ifdef MOLFORCES
    MPI_Bcast(&(topology[i].trap_flag),1,MPI_INT,0,comm_cart);
    MPI_Bcast(topology[i].trap_center,3,MPI_DOUBLE,0,comm_cart);
    MPI_Bcast(&(topology[i].trap_spring_constant),1,MPI_DOUBLE,0,comm_cart);
    MPI_Bcast(&(topology[i].drag_constant),1,MPI_DOUBLE,0,comm_cart);
    MPI_Bcast(&(topology[i].noforce_flag),1,MPI_INT,0,comm_cart);
    MPI_Bcast(&(topology[i].isrelative),1,MPI_INT,0,comm_cart);
    MPI_Bcast(&(topology[i].favcounter),1,MPI_INT,0,comm_cart);
    if (topology[i].favcounter == -1)
      MPI_Bcast(topology[i].fav,3,MPI_DOUBLE,0,comm_cart);
    /* check if any molecules are trapped */
    if  ((topology[i].trap_flag != 32) && (topology[i].noforce_flag != 32)) {
      IsTrapped = 1;
    }
#endif

    MPI_Bcast(&molsize,1,MPI_INT,0,comm_cart);
    MPI_Bcast(&moltype,1,MPI_INT,0,comm_cart);
    MPI_Bcast(topology[i].part.e,topology[i].part.n,MPI_INT,0,comm_cart);
    MPI_Bcast(&topology[i].type,1,MPI_INT,0,comm_cart);
    
  }
  
  sync_topo_part_info();

  return 1;
}

void mpi_sync_topo_part_info_slave(int node,int parm ) {
  int i;
  int molsize=0;
  int moltype=0;
  int n_mols=0;

  MPI_Bcast(&n_mols,1,MPI_INT,0,comm_cart);
  realloc_topology(n_mols);
  for ( i = 0 ; i < n_molecules ; i++) {

#ifdef MOLFORCES
    MPI_Bcast(&(topology[i].trap_flag),1,MPI_INT,0,comm_cart);
    MPI_Bcast(topology[i].trap_center,3,MPI_DOUBLE,0,comm_cart);
    MPI_Bcast(&(topology[i].trap_spring_constant),1,MPI_DOUBLE,0,comm_cart);
    MPI_Bcast(&(topology[i].drag_constant),1,MPI_DOUBLE,0,comm_cart);
    MPI_Bcast(&(topology[i].noforce_flag),1,MPI_INT,0,comm_cart);
    MPI_Bcast(&(topology[i].isrelative),1,MPI_INT,0,comm_cart);
    MPI_Bcast(&(topology[i].favcounter),1,MPI_INT,0,comm_cart);
    if (topology[i].favcounter == -1)
      MPI_Bcast(topology[i].fav,3,MPI_DOUBLE,0,comm_cart);
    /* check if any molecules are trapped */
    if  ((topology[i].trap_flag != 32) && (topology[i].noforce_flag != 32)) {
      IsTrapped = 1;
    }
#endif

    MPI_Bcast(&molsize,1,MPI_INT,0,comm_cart);
    MPI_Bcast(&moltype,1,MPI_INT,0,comm_cart);
    topology[i].type = moltype;
    realloc_intlist(&topology[i].part,topology[i].part.n = molsize);

    MPI_Bcast(topology[i].part.e,topology[i].part.n,MPI_INT,0,comm_cart);
    MPI_Bcast(&topology[i].type,1,MPI_INT,0,comm_cart);

  }
  

  sync_topo_part_info();

}

/******************* REQ_BCAST_LBPAR ********************/

void mpi_bcast_lb_params(int field)
{
#ifdef LB
  mpi_call(mpi_bcast_lb_params_slave, -1, field);
  mpi_bcast_lb_params_slave(-1, field);
#endif
}

void mpi_bcast_lb_params_slave(int node, int field)
{
#ifdef LB
  MPI_Bcast(&lbpar, sizeof(LB_Parameters), MPI_BYTE, 0, comm_cart);
  on_lb_params_change(field);
#endif
}


/******************* REQ_BCAST_CUDA_GLOBAL_PART_VARS ********************/

void mpi_bcast_cuda_global_part_vars() {
#ifdef CUDA
  mpi_call(mpi_bcast_cuda_global_part_vars_slave, 1, 0); // third parameter is meaningless
  mpi_bcast_cuda_global_part_vars_slave(-1,0);
#endif
}

void mpi_bcast_cuda_global_part_vars_slave(int node, int dummy)
{
#ifdef CUDA
  MPI_Bcast(gpu_get_global_particle_vars_pointer_host(), sizeof(CUDA_global_part_vars), MPI_BYTE, 0, comm_cart);
  espressoSystemInterface.requestParticleStructGpu();
#endif
}



/******************* REQ_GET_ERRS ********************/

int mpi_gather_runtime_errors(char **errors)
{
  // Tell other processors to send their erros
  mpi_call(mpi_gather_runtime_errors_slave, -1, 0);

  // If no processor encountered an error, return
  if (!check_runtime_errors())
    return ES_OK;

  // gather the maximum length of the error messages
  int *errcnt = (int*)malloc(n_nodes*sizeof(int));
  MPI_Gather(&n_error_msg, 1, MPI_INT, errcnt, 1, MPI_INT, 0, comm_cart);

  for (int node = 0; node < n_nodes; node++) {
    if (errcnt[node] > 0) {
      errors[node] = (char *)malloc(errcnt[node]);

      if (node == 0)
	strcpy(errors[node], error_msg);
      else 
	MPI_Recv(errors[node], errcnt[node], MPI_CHAR, node, 0, comm_cart, MPI_STATUS_IGNORE);
    }
    else
      errors[node] = NULL;
  }

  /* reset error message on master node */
  error_msg = (char*)realloc(error_msg, n_error_msg = 0);

  free(errcnt);

  return ES_ERROR;
}

void mpi_gather_runtime_errors_slave(int node, int parm)
{
  if (!check_runtime_errors())
    return;

  MPI_Gather(&n_error_msg, 1, MPI_INT, NULL, 0, MPI_INT, 0, comm_cart);
  if (n_error_msg > 0) {
    MPI_Send(error_msg, n_error_msg, MPI_CHAR, 0, 0, comm_cart);
    /* reset error message on slave node */
    error_msg = (char*)realloc(error_msg, n_error_msg = 0);
  }
}

/********************* REQ_SET_EXCL ********/
void mpi_send_exclusion(int part1, int part2, int _delete)
{
#ifdef EXCLUSIONS
  mpi_call(mpi_send_exclusion_slave, part1, part2);

  MPI_Bcast(&_delete, 1, MPI_INT, 0, comm_cart);
  local_change_exclusion(part1, part2, _delete);
  on_particle_change();
#endif
}

void mpi_send_exclusion_slave(int part1, int part2)
{
#ifdef EXCLUSIONS
  int _delete=0;
  MPI_Bcast(&_delete, 1, MPI_INT, 0, comm_cart);  
  local_change_exclusion(part1, part2, _delete);
  on_particle_change();
#endif
}

/************** REQ_SET_FLUID **************/
void mpi_send_fluid(int node, int index, double rho, double *j, double *pi) {
#ifdef LB
  if (node==this_node) {
    lb_calc_n_from_rho_j_pi(index, rho, j, pi);
  } else {
    double data[10] = { rho, j[0], j[1], j[2], pi[0], pi[1], pi[2], pi[3], pi[4], pi[5] };
    mpi_call(mpi_send_fluid_slave, node, index);
    MPI_Send(data, 10, MPI_DOUBLE, node, SOME_TAG, comm_cart);
  }
#endif
}

void mpi_send_fluid_slave(int node, int index) {
#ifdef LB
  if (node==this_node) {
    double data[10];
    MPI_Recv(data, 10, MPI_DOUBLE, 0, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);

    lb_calc_n_from_rho_j_pi(index, data[0], &data[1], &data[4]);
  }
#endif
}

/************** REQ_GET_FLUID **************/
void mpi_recv_fluid(int node, int index, double *rho, double *j, double *pi) {
#ifdef LB
  if (node==this_node) {
    lb_calc_local_fields(index, rho, j, pi);
  } else {
    double data[10];
    mpi_call(mpi_recv_fluid_slave, node, index);
        MPI_Recv(data, 10, MPI_DOUBLE, node, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
    *rho = data[0];
    j[0] = data[1];
    j[1] = data[2];
    j[2] = data[3];
    pi[0] = data[4];
    pi[1] = data[5];
    pi[2] = data[6];
    pi[3] = data[7];
    pi[4] = data[8];
    pi[5] = data[9];

  }
#endif
}

void mpi_recv_fluid_slave(int node, int index) {
#ifdef LB
  if (node==this_node) {
    double data[10];
    lb_calc_local_fields(index, &data[0], &data[1], &data[4]);
    MPI_Send(data, 10, MPI_DOUBLE, 0, SOME_TAG, comm_cart);
  }
#endif
}

/************** REQ_LB_GET_BOUNDARY_FLAG **************/
void mpi_recv_fluid_boundary_flag(int node, int index, int *boundary) {
#ifdef LB_BOUNDARIES
  if (node==this_node) {
    lb_local_fields_get_boundary_flag(index, boundary);
  } else {
    int data = 0;
    mpi_call(mpi_recv_fluid_boundary_flag_slave, node, index);
    MPI_Recv(&data, 1, MPI_INT, node, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
    *boundary = data;
  }
#endif
}

void mpi_recv_fluid_boundary_flag_slave(int node, int index) {
#ifdef LB_BOUNDARIES
  if (node==this_node) {
    int data;
    lb_local_fields_get_boundary_flag(index, &data);
    MPI_Send(&data, 1, MPI_INT, 0, SOME_TAG, comm_cart);
  }
#endif
}

/********************* REQ_ICCP3M_ITERATION ********/
int mpi_iccp3m_iteration(int dummy)
{
#ifdef ELECTROSTATICS
  mpi_call(mpi_iccp3m_iteration_slave, -1, 0);

  iccp3m_iteration();

  COMM_TRACE(fprintf(stderr, "%d: iccp3m iteration task %d done.\n", this_node, dummy));

  return check_runtime_errors();
#else
  return 0;
#endif

}

void mpi_iccp3m_iteration_slave(int dummy, int dummy2)
{
#ifdef ELECTROSTATICS
  iccp3m_iteration();
  COMM_TRACE(fprintf(stderr, "%d: iccp3m iteration task %d done.\n", dummy, dummy2));

  check_runtime_errors();
#endif
}

/********************* REQ_ICCP3M_INIT********/
int mpi_iccp3m_init(int n_induced_charges)
{
#ifdef ELECTROSTATICS
  /* nothing has to be done on the master node, this 
   * passes only the number of induced charges, in order for
   * slaves to allocate memory */
  
  mpi_call(mpi_iccp3m_init_slave, -1, n_induced_charges);
   
  bcast_iccp3m_cfg();

  COMM_TRACE(fprintf(stderr, "%d: iccp3m init task %d done.\n", this_node, n_induced_charges));

  return check_runtime_errors();
#else
  return 0;
#endif

}

void mpi_iccp3m_init_slave(int node, int dummy)
{
#ifdef ELECTROSTATICS
  COMM_TRACE(fprintf(stderr, "%d: iccp3m iteration task %d done.\n", this_node, dummy));

  if(iccp3m_initialized==0){
    iccp3m_init();
    iccp3m_initialized=1;
 }

  bcast_iccp3m_cfg();
  
  check_runtime_errors();
#endif
}

void mpi_recv_fluid_populations(int node, int index, double *pop) {
#ifdef LB
  if (node==this_node) {
    lb_get_populations(index, pop);
  } else {
    mpi_call(mpi_recv_fluid_populations_slave, node, index);
    MPI_Recv(pop, 19, MPI_DOUBLE, node, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
  }
  lbpar.resend_halo=1;
#endif
}

void mpi_recv_fluid_populations_slave(int node, int index) {
#ifdef LB
  if (node==this_node) {
    double data[19];
    lb_get_populations(index, data);
    MPI_Send(data, 19, MPI_DOUBLE, 0, SOME_TAG, comm_cart);
  }
  lbpar.resend_halo=1;
#endif
}

void mpi_send_fluid_populations(int node, int index, double *pop) {
#ifdef LB
  if (node==this_node) {
    lb_set_populations(index, pop);
  } else {
    mpi_call(mpi_send_fluid_populations_slave, node, index);
    MPI_Send(pop, 19, MPI_DOUBLE, node, SOME_TAG, comm_cart);
  }
#endif
}

void mpi_send_fluid_populations_slave(int node, int index) {
#ifdef LB
  if (node==this_node) {
    double data[19];
    MPI_Recv(data, 19, MPI_DOUBLE, 0, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
    lb_set_populations(index, data);
  }
#endif
}

/****************************************************/

void mpi_bcast_max_mu() {
#ifdef DIPOLES
  mpi_call(mpi_bcast_max_mu_slave, -1, 0);
  
  calc_mu_max();
  
#endif
}

void mpi_bcast_max_mu_slave(int node, int dummy) {
#ifdef DIPOLES
  
  calc_mu_max();
 
#endif
}

#ifdef LANGEVIN_PER_PARTICLE

/******************** REQ_SEND_PARTICLE_T ********************/
void mpi_set_particle_temperature(int pnode, int part, double _T)
{
  mpi_call(mpi_set_particle_temperature_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    /* here the setting actually happens, if the particle belongs to the local node */
    p->p.T = _T;
  }
  else {
    MPI_Send(&_T, 1, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
}
#endif

void mpi_set_particle_temperature_slave(int pnode, int part)
{
#ifdef LANGEVIN_PER_PARTICLE
  double s_buf = 0.;
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Status status;
    MPI_Recv(&s_buf, 1, MPI_DOUBLE, 0, SOME_TAG, comm_cart, &status);
    /* here the setting happens for nonlocal nodes */
    p->p.T = s_buf;
  }

  on_particle_change();
#endif
}

#ifdef LANGEVIN_PER_PARTICLE
void mpi_set_particle_gamma(int pnode, int part, double gamma)
{
  mpi_call(mpi_set_particle_gamma_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    /* here the setting actually happens, if the particle belongs to the local node */
    p->p.gamma = gamma;
  }
  else {
    MPI_Send(&gamma, 1, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
}
#endif

void mpi_set_particle_gamma_slave(int pnode, int part)
{
#ifdef LANGEVIN_PER_PARTICLE
  double s_buf = 0.;
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Status status;
    MPI_Recv(&s_buf, 1, MPI_DOUBLE, 0, SOME_TAG, comm_cart, &status);
    /* here the setting happens for nonlocal nodes */
    p->p.gamma = s_buf;
  }

  on_particle_change();
#endif
}

/***** GALILEI TRANSFORM AND ASSOCIATED FUNCTIONS ****/

void mpi_kill_particle_motion( int rotation )
{
  mpi_call(mpi_kill_particle_motion_slave, -1, rotation );
  local_kill_particle_motion( rotation );
  on_particle_change();
}

void mpi_kill_particle_motion_slave(int pnode, int rotation )
{
  local_kill_particle_motion( rotation );
  on_particle_change();
}

void mpi_kill_particle_forces( int torque )
{
  mpi_call(mpi_kill_particle_forces_slave, -1, torque );
  local_kill_particle_forces( torque );
  on_particle_change();
}

void mpi_kill_particle_forces_slave(int pnode, int torque)
{
  local_kill_particle_forces( torque );
  on_particle_change();
}

void mpi_system_CMS() {
  int pnode;
  double data[4];
  double rdata[4];
  double *pdata = rdata;

  data[0] = 0.0; 
  data[1] = 0.0;
  data[2] = 0.0;
  data[3] = 0.0;

  mpi_call(mpi_system_CMS_slave, -1, 0);

  for (pnode = 0; pnode < n_nodes; pnode++) {
    if (pnode==this_node) {
      local_system_CMS( pdata );
      data[0] += rdata[0];
      data[1] += rdata[1];
      data[2] += rdata[2];
      data[3] += rdata[3];
    } else {
      MPI_Recv(rdata, 4, MPI_DOUBLE, MPI_ANY_SOURCE, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
      data[0] += rdata[0];
      data[1] += rdata[1];
      data[2] += rdata[2];
      data[3] += rdata[3];
    }
  }

  gal.cms[0] = data[0]/data[3];
  gal.cms[1] = data[1]/data[3];
  gal.cms[2] = data[2]/data[3];
}

void mpi_system_CMS_slave(int node, int index) {
  double rdata[4];
  double *pdata = rdata;
  local_system_CMS( pdata );
  MPI_Send(rdata, 4, MPI_DOUBLE, 0, SOME_TAG, comm_cart);
}

void mpi_system_CMS_velocity() {
  int pnode;
  double data[4];
  double rdata[4];
  double *pdata = rdata;

  data[0] = 0.0; 
  data[1] = 0.0;
  data[2] = 0.0;
  data[3] = 0.0;

  mpi_call(mpi_system_CMS_velocity_slave, -1, 0);

  for (pnode = 0; pnode < n_nodes; pnode++) {
    if (pnode==this_node) {
      local_system_CMS_velocity( pdata );
      data[0] += rdata[0];
      data[1] += rdata[1];
      data[2] += rdata[2];
      data[3] += rdata[3];
    } else {
      MPI_Recv(rdata, 4, MPI_DOUBLE, MPI_ANY_SOURCE, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
      data[0] += rdata[0];
      data[1] += rdata[1];
      data[2] += rdata[2];
      data[3] += rdata[3];
    }
  }

  gal.cms_vel[0] = data[0]/data[3];
  gal.cms_vel[1] = data[1]/data[3];
  gal.cms_vel[2] = data[2]/data[3];
}

void mpi_system_CMS_velocity_slave(int node, int index) {
  double rdata[4];
  double *pdata = rdata;
  local_system_CMS_velocity( pdata );
  MPI_Send(rdata, 4, MPI_DOUBLE, 0, SOME_TAG, comm_cart);
}

void mpi_galilei_transform()
{
  double cmsvel[3];

  mpi_system_CMS_velocity();
  memcpy(cmsvel, gal.cms_vel, 3*sizeof(double));

  mpi_call(mpi_galilei_transform_slave, -1, 0);
  MPI_Bcast(cmsvel, 3, MPI_DOUBLE, 0, comm_cart);

  local_galilei_transform( cmsvel );

  on_particle_change();
}

void mpi_galilei_transform_slave(int pnode, int i)
{
  double cmsvel[3];
  MPI_Bcast(cmsvel, 3, MPI_DOUBLE, 0, comm_cart);

  local_galilei_transform( cmsvel );
  on_particle_change();
}

/******************** REQ_CATALYTIC_REACTIONS ********************/

void mpi_setup_reaction()
{
#ifdef CATALYTIC_REACTIONS
  mpi_call(mpi_setup_reaction_slave, -1, 0);
  local_setup_reaction();
#endif
}

void mpi_setup_reaction_slave(int pnode, int i)
{
#ifdef CATALYTIC_REACTIONS
  local_setup_reaction();
#endif
}

/*********************** MAIN LOOP for slaves ****************/

void mpi_loop()
{
  for (;;) {
#ifdef ASYNC_BARRIER
    MPI_Barrier(comm_cart);
#endif
    MPI_Bcast(request, 3, MPI_INT, 0, comm_cart);
    COMM_TRACE(fprintf(stderr, "%d: processing %s %d %d...\n", this_node,
		       names[request[0]], request[1], request[2]));
    if ((request[0] < 0) || (request[0] >= N_CALLBACKS)) {
      fprintf(stderr, "%d: INTERNAL ERROR: unknown request %d\n", this_node, request[0]);
      errexit();
    }
    slave_callbacks[request[0]](request[1], request[2]);
    COMM_TRACE(fprintf(stderr, "%d: finished %s %d %d\n", this_node,
		       names[request[0]], request[1], request[2]));
  
  }
}

void mpi_external_potential_broadcast(int number) {
  mpi_call(mpi_external_potential_broadcast_slave, 0, number);
  MPI_Bcast(&external_potentials[number], sizeof(ExternalPotential), MPI_BYTE, 0, comm_cart);
  MPI_Bcast(external_potentials[number].scale, external_potentials[number].n_particle_types, MPI_DOUBLE, 0, comm_cart);
}

void mpi_external_potential_broadcast_slave(int node, int number) {
  ExternalPotential E;
  MPI_Bcast(&E, sizeof(ExternalPotential), MPI_BYTE, 0, comm_cart);
  ExternalPotential* new_;
  generate_external_potential(&new_);
  external_potentials[number] = E;
  external_potentials[number].scale=(double*) malloc(external_potentials[number].n_particle_types);
  MPI_Bcast(external_potentials[number].scale, external_potentials[number].n_particle_types, MPI_DOUBLE, 0, comm_cart);
}

void mpi_external_potential_tabulated_read_potential_file(int number) {
  mpi_call(mpi_external_potential_tabulated_read_potential_file_slave, 0, number);
  external_potential_tabulated_read_potential_file(number);
}

void mpi_external_potential_tabulated_read_potential_file_slave(int node, int number) {
  external_potential_tabulated_read_potential_file(number);
}

void mpi_external_potential_sum_energies() {
  mpi_call(mpi_external_potential_sum_energies_slave, 0, 0);
  double* energies = (double*) malloc(n_external_potentials*sizeof(double));
  for (int i=0; i<n_external_potentials; i++) {
    energies[i]=external_potentials[i].energy;
  }
  double* energies_sum = (double*) malloc(n_external_potentials*sizeof(double)); 
  MPI_Reduce(energies, energies_sum, n_external_potentials, MPI_DOUBLE, MPI_SUM, 0, comm_cart); 
  for (int i=0; i<n_external_potentials; i++) {
    external_potentials[i].energy=energies_sum[i];
  }
  free(energies);
  free(energies_sum);
}


void mpi_external_potential_sum_energies_slave(int dummy1, int dummy2) {
  double* energies = (double*)malloc(n_external_potentials*sizeof(double));
  for (int i=0; i<n_external_potentials; i++) {
    energies[i]=external_potentials[i].energy;
  }
  MPI_Reduce(energies, 0, n_external_potentials, MPI_DOUBLE, MPI_SUM, 0, comm_cart); 
  free(energies);
}
