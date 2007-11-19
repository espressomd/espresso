// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "domain_decomposition.h"
#include "rattle.h"

int n_rigidbonds = 0;

#ifdef BOND_CONSTRAINT

/** \name Private functions */
/************************************************************/
/*@{*/

/** Calculates the corrections required for each of the particle coordinates
    according to the RATTLE algorithm. Invoked from \ref correct_pos_shake()*/
void compute_pos_corr_vec(int *repeat_);

/** Positional Corrections are added to the current particle positions. Invoked from \ref correct_pos_shake() */
void app_pos_correction();

/* Checks if verlet list rebuilt is required.*/
void check_rebuild_verletlist();

/** Transfers temporarily the current forces from f.f[3] of the \ref Particle
    structure to r.p_old[3] location and also intializes velocity correction
    vector. Invoked from \ref correct_vel_shake()*/
void transfer_force_init_vel();

/** Calculates corrections of the  current particle velocities according to RATTLE
    algorithm. Invoked from \ref correct_vel_shake()*/
void compute_vel_corr_vec(int *repeat_);

/** Velocity corrections are added to the current particle velocities. Invoked from
    \ref correct_vel_shake()*/
void apply_vel_corr();


/**Invoked from \ref correct_vel_shake(). Put back the forces from r.p_old to f.f*/
void revert_force();

/**For debugging purpose--prints the bond lengths between particles that have
rigid_bonds*/
void print_bond_len();

/*@}*/

/*Initialize old positions (particle positions at previous time step)
  of the particles*/
void save_old_pos()
{
  int c, i, j, np;
  Particle *p;
  Cell *cell;
  for (c = 0; c < local_cells.n; c++)
  {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      for(j=0;j<3;j++)
        p[i].r.p_old[j]=p[i].r.p[j];
    } //for i loop
  }// for c loop

  for (c = 0; c < ghost_cells.n; c++)
  {
    cell = ghost_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
       for(j=0;j<3;j++)
          p[i].r.p_old[j]=p[i].r.p[j];
    }
  }
}

/**Initialize the correction vector. The correction vector is stored in f.f of particle strcuture. */
void init_correction_vector()
{

  int c, i, j, np;
  Particle *p;
  Cell *cell;

  for (c = 0; c < local_cells.n; c++)
  {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      for(j=0;j<3;j++)
        p[i].f.f[j] = 0.0;
     } //for i loop
  }// for c loop

  for (c = 0; c < ghost_cells.n; c++)
  {
    cell = ghost_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
       for(j=0;j<3;j++)
          p[i].f.f[j] = 0.0;
    }
  }
}

void init_shake_p()
{

  int c, i, j, np;
  Particle *p;
  Cell *cell;

  for (c = 0; c < local_cells.n; c++)
  {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      for(j=0;j<3;j++)
        p[i].r.shake_p=0.0;
     } //for i loop
  }// for c loop
}

/**Compute positional corrections*/
void compute_pos_corr_vec(int *repeat_)
{
  Bonded_ia_parameters *ia_params;
  int i, j, k, c, np, cnt=-1;
  Cell *cell;
  Particle *p, *p1, *p2;
  double r_ij_t[3], r_ij[3], r_ij_dot, G, pos_corr, r_ij2;

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      p1 = &p[i];
      k=0;
      while(k<p1->bl.n) {
	ia_params = &bonded_ia_params[p1->bl.e[k++]];
	if( ia_params->type == BONDED_IA_RIGID_BOND ) {
	  cnt++;
	  p2 = local_particles[p1->bl.e[k++]];
	  if (!p2) {
	    char *errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
	    ERROR_SPRINTF(errtxt,"{051 rigid bond broken between particles %d and %d (particles not stored on the same node)} ",
		    p1->p.identity, p1->bl.e[k-1]);
	    return;
	  }

	  get_mi_vector(r_ij  , p1->r.p    , p2->r.p    );
	  r_ij2 = sqrlen(r_ij);
	  if(fabs(1.0 - r_ij2/ia_params->p.rigid_bond.d2) > ia_params->p.rigid_bond.p_tol) {
            get_mi_vector(r_ij_t, p1->r.p_old, p2->r.p_old);
	    r_ij_dot = scalar(r_ij_t, r_ij);
	    G = 0.50*(ia_params->p.rigid_bond.d2 - r_ij2 )/r_ij_dot;
#ifdef MASS
	    G /= (PMASS(*p1)+PMASS(*p2));
#else
	    G /= 2;
#endif
	    for (j=0;j<3;j++) {
	      pos_corr = G*r_ij_t[j];
	      p1->f.f[j] += pos_corr*PMASS(*p2);
	      p2->f.f[j] -= pos_corr*PMASS(*p1);
	      local_particles[p1->p.identity]->r.shake_p+=pos_corr*r_ij_t[j];
	    }
	    /*Increase the 'repeat' flag by one */
	      *repeat_ = *repeat_ + 1;
	  }
	}
	else
	  /* skip bond partners of nonrigid bond */
          k+=ia_params->num;
      } //while loop
    } //for i loop
  } //for c loop
}

/**Apply corrections to each particle**/
void app_pos_correction()
{
  int c, i, j, np;
  Particle *p, *p1;
  Cell *cell;


     /*Apply corrections*/
     for (c = 0; c < local_cells.n; c++)
     {
      cell = local_cells.cell[c];
      p  = cell->part;
      np = cell->n;
      for(i = 0; i < np; i++) {
        p1 = &p[i];
	for (j=0;j<3;j++) {
	   p1->r.p[j] += p1->f.f[j];
	   p1->m.v[j] += p1->f.f[j];
	}

	/**Completed for one particle*/
       } //for i loop
     } //for c loop
}

/* Check if verlet list is required to be re-built*/
void check_rebuild_verletlist()
{
     int i, c, np;
     Cell *cell;
     Particle *p;
     double skin2 = SQR(skin/2.0);
     rebuild_verletlist=0;
     for (c = 0; c < local_cells.n; c++)
     {
      cell = local_cells.cell[c];
      p  = cell->part;
      np = cell->n;
      for(i = 0; i < np; i++) {
        /* Verlet criterion check */
	if(distance2(p[i].r.p, p[i].l.p_old) > skin2 ) rebuild_verletlist = 1;
       } //for i loop
     } //for c loop
     announce_rebuild_vlist();
}


void correct_pos_shake()
{
   int    repeat_,  cnt=0;
   int repeat=1;

   while (repeat!=0 && cnt<SHAKE_MAX_ITERATIONS)
   {
     init_correction_vector();
     repeat_ = 0;
     compute_pos_corr_vec(&repeat_);
     ghost_communicator(&cell_structure.collect_ghost_force_comm);
     app_pos_correction();
     /**Ghost Positions Update*/
     ghost_communicator(&cell_structure.update_ghost_pos_comm);
     if(this_node==0)
         MPI_Reduce(&repeat_, &repeat, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
     else
         MPI_Reduce(&repeat_, NULL, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
     MPI_Bcast(&repeat, 1, MPI_INT, 0, MPI_COMM_WORLD);

       cnt++;
   }// while(repeat) loop
   if (cnt >= SHAKE_MAX_ITERATIONS) {
     char *errtxt = runtime_error(100 + TCL_INTEGER_SPACE);
     ERROR_SPRINTF(errtxt, "{053 RATTLE failed to converge after %d iterations} ", cnt);
   }

   check_rebuild_verletlist();
}

/**The forces are transfered temporarily from f.f member of particle structure to r.p_old,
    which is idle now and initialize the velocity correction vector to zero at f.f[3]
    of Particle structure*/
void transfer_force_init_vel()
{

  int c, i, j, np;
  Particle *p;
  Cell *cell;

  for (c = 0; c < local_cells.n; c++)
    {
      cell = local_cells.cell[c];
      p  = cell->part;
      np = cell->n;
      for(i = 0; i < np; i++) {
	for(j=0;j<3;j++)
	  {
	    p[i].r.p_old[j]=p[i].f.f[j];
	    p[i].f.f[j]=0.0;
	  }
      } //for i loop
    }// for c loop

  for (c = 0; c < ghost_cells.n; c++)
    {
      cell = ghost_cells.cell[c];
      p  = cell->part;
      np = cell->n;
      for(i = 0; i < np; i++) {
	for(j=0;j<3;j++)
	  {
	    p[i].r.p_old[j]=p[i].f.f[j];
	    p[i].f.f[j]=0.0;
	  }
      }
    }
}

/** Velocity correction vectors are computed*/
void compute_vel_corr_vec(int *repeat_)
{
  Bonded_ia_parameters *ia_params;
  int i, j, k, c, np;
  Cell *cell;
  Particle *p, *p1, *p2;
  double v_ij[3], r_ij[3], K, vel_corr;

  for (c = 0; c < local_cells.n; c++)
    {
      cell = local_cells.cell[c];
      p  = cell->part;
      np = cell->n;
      for(i = 0; i < np; i++) {
        p1 = &p[i];
        k=0;
	while(k<p1->bl.n)
	  {
	    ia_params = &bonded_ia_params[p1->bl.e[k++]];
	    if( ia_params->type == BONDED_IA_RIGID_BOND )
	      {
		p2 = local_particles[p1->bl.e[k++]];
		if (!p2) {
		  char *errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
		  ERROR_SPRINTF(errtxt,"{054 rigid bond broken between particles %d and %d (particles not stored on the same node)} ",
			  p1->p.identity, p1->bl.e[k-1]);
		  return;
		}

		vecsub(p1->m.v, p2->m.v, v_ij);
                get_mi_vector(r_ij, p1->r.p, p2->r.p);
		if(fabs(scalar(v_ij, r_ij)) > ia_params->p.rigid_bond.v_tol)
	        {
		  K = scalar(v_ij, r_ij)/ia_params->p.rigid_bond.d2;
#ifdef MASS
		  K /= (PMASS(*p1) + PMASS(*p2));
#else
		  K /= 2.0;
#endif
		for (j=0;j<3;j++)
		  {
		    vel_corr = K*r_ij[j];
		    p1->f.f[j] -= vel_corr*PMASS(*p2);
		    p2->f.f[j] += vel_corr*PMASS(*p1);
		    local_particles[p1->p.identity]->r.shake_p-=vel_corr*r_ij[j];

		  }
		  *repeat_ = *repeat_ + 1 ;
		}
	      }
	    else
	      k += ia_params->num;
	  } //while loop
      } //for i loop
    } //for c loop
}

/**Apply velocity corrections*/
void apply_vel_corr()
{

  int c, i, j, np;
  Particle *p, *p1;
  Cell *cell;

     /*Apply corrections*/
     for (c = 0; c < local_cells.n; c++)
     {
      cell = local_cells.cell[c];
      p  = cell->part;
      np = cell->n;
      for(i = 0; i < np; i++) {
        p1 = &p[i];
	for (j=0;j<3;j++) {
	   p1->m.v[j] += p1->f.f[j];

	}
	/**Completed for one particle*/
       } //for i loop
     } //for c loop

}

/**Put back the forces from r.p_old to f.f*/
void revert_force()
{
  int c, i, j, np;
  Particle *p;
  Cell *cell;
  for (c = 0; c < local_cells.n; c++)
  {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      for(j=0;j<3;j++)
      	p[i].f.f[j]=p[i].r.p_old[j];

     } //for i loop
  }// for c loop


  for (c = 0; c < ghost_cells.n; c++)
  {
    cell = ghost_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++)
    {
       for(j=0;j<3;j++)
  	  p[i].f.f[j]=p[i].r.p_old[j];
    }
  }
}

void correct_vel_shake()
{
   int    repeat_, repeat=1, cnt=0;
   /**transfer the current forces to r.p_old of the particle structure so that
   velocity corrections can be stored temporarily at the f.f[3] of the particle
   structure  */
   transfer_force_init_vel();
   while (repeat!=0 && cnt<SHAKE_MAX_ITERATIONS)
   {
     init_correction_vector();
     repeat_ = 0;
     compute_vel_corr_vec(&repeat_);
     ghost_communicator(&cell_structure.collect_ghost_force_comm);
     apply_vel_corr();
     ghost_communicator(&cell_structure.update_ghost_pos_comm);
     if(this_node == 0)
       MPI_Reduce(&repeat_, &repeat, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
     else
       MPI_Reduce(&repeat_, NULL, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

     MPI_Bcast(&repeat, 1, MPI_INT, 0, MPI_COMM_WORLD);
     cnt++;
   }

   if (cnt >= SHAKE_MAX_ITERATIONS) {
         fprintf(stderr, "%d: VEL CORRECTIONS IN RATTLE failed to converge after %d iterations !!\n", this_node, cnt);
	 errexit();
   }
   /**Puts back the forces from r.p_old to f.f[3]*/
   revert_force();
}


void print_bond_len()
{
  int i, k, c, np;
  Cell *cell;
  Particle *p;
  Bonded_ia_parameters *b_ia;
  double r_ij[3];
  printf("%d: ", this_node);
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
       k=0;
       while(k<p[i].bl.n)
       {
             b_ia = &bonded_ia_params[p[i].bl.e[k]];
	     if(b_ia->type == BONDED_IA_RIGID_BOND)
             {
	       Particle *p2 = local_particles[p[i].bl.e[k++]];
	       if (!p2) {
		 char *errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
		 ERROR_SPRINTF(errtxt,"{056 rigid bond broken between particles %d and %d (particles not stored on the same node)} ", p[i].p.identity, p[i].bl.e[k-1]);
		 return;
	       }

	       get_mi_vector(r_ij, p[i].r.p , p2->r.p);
	       printf(" bl (%d %d): %f\t", p[i].p.identity, p2->p.identity,sqrlen(r_ij));
             }
	     else
	       k += b_ia->num;
       } //while k loop
    } //for i loop
  }// for c loop
  printf("\n");
}

/*****************************************************************************
 *   setting parameters
 *****************************************************************************/
int rigid_bond_set_params(int bond_type, double d, double p_tol, double v_tol)
{
  if(bond_type < 0)
    return TCL_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].p.rigid_bond.d2 = d*d;
  bonded_ia_params[bond_type].p.rigid_bond.p_tol = 2.0*p_tol;
  bonded_ia_params[bond_type].p.rigid_bond.v_tol = v_tol*time_step;
  bonded_ia_params[bond_type].type = BONDED_IA_RIGID_BOND;
  bonded_ia_params[bond_type].num = 1;
  n_rigidbonds += 1;
  mpi_bcast_ia_params(bond_type, -1);
  mpi_bcast_parameter(FIELD_RIGIDBONDS);

  return TCL_OK;
}

#endif
