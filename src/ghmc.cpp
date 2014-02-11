/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
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
/** \file ghmc.cpp

    For more information see \ref ghmc.hpp
 */

#include <cmath>
#include "utils.hpp"

#include "ghmc.hpp"
#include "random.hpp"
#include "particle_data.hpp"
#include "communication.hpp"
#include "energy.hpp"
#include "cells.hpp"
#include "virtual_sites.hpp"
#include "statistics.hpp"

/************************************************************/

//momentum flip flag
int ghmc_mflip = GHMC_MFLIP_OFF;
//temperature scaling flag
int ghmc_tscale = GHMC_TSCALE_OFF;
//result of mc decision
int ghmc_mc_res;

#ifdef GHMC
//MC statistics variables
int ghmc_att, ghmc_acc;
//partial momentum update parameters
double cosp, sinp;
//inverse temperature
double beta;
//thermostat data struct
Ghmc ghmcdata = { 0, 0, 0.0, 0.0};
#endif

/************************************************************/
/* local prototypes                                         */
/************************************************************/

void hamiltonian_calc(int ekin_update_flag);

double calc_local_temp();

void calc_kinetic(double *ek_trans , double *ek_rot);

void simple_momentum_update();

void partial_momentum_update();

void tscale_momentum_update();

void momentum_flip();
/************************************************************/

/************************************************************/
/** \name Privat Functions */
/************************************************************/
/*@{*/

#ifdef GHMC

void save_last_state()
{
	
	int i, c, np;
	Particle *part;
	
	//part = local_cells.cell[0]->part;
	//fprintf(stderr,"%d: save part %d: px_ls before %f, px before %f\n",this_node,part[0].p.identity,part[0].l.r_ls.p[0],part[0].r.p[0]);
	//fprintf(stderr,"%d: save part %d: mx_ls before %f, mx before %f\n",this_node,part[0].p.identity,part[0].l.m_ls.v[0],part[0].m.v[0]);
	
  for (c = 0; c < local_cells.n; c++) {
    np   = local_cells.cell[c]->n;
    part = local_cells.cell[c]->part;
    for (i = 0; i < np; i++) {
			memcpy(&part[i].l.r_ls, &part[i].r, sizeof(ParticlePosition));
			memcpy(&part[i].l.m_ls, &part[i].m, sizeof(ParticleMomentum));
		}
	}
	
		//part = local_cells.cell[0]->part;
		//fprintf(stderr,"%d: save part %d: px_ls after %f, px after %f\n",this_node,part[0].p.identity,part[0].l.r_ls.p[0],part[0].r.p[0]);
		//fprintf(stderr,"%d: save part %d: mx_ls after %f, mx after %f\n",this_node,part[0].p.identity,part[0].l.m_ls.v[0],part[0].m.v[0]);
		
}

void load_last_state()
{
	
	int i, c, np;
	Particle *part;
	
	//part = local_cells.cell[0]->part;
	//fprintf(stderr,"%d: load part %d: px_ls before %f, px before %f\n",this_node,part[0].p.identity,part[0].l.r_ls.p[0],part[0].r.p[0]);
	//fprintf(stderr,"%d: load part %d: mx_ls before %f, mx before %f\n",this_node,part[0].p.identity,part[0].l.m_ls.v[0],part[0].m.v[0]);
	
  for (c = 0; c < local_cells.n; c++) {
    np   = local_cells.cell[c]->n;
    part = local_cells.cell[c]->part;
    for (i = 0; i < np; i++) {
			memcpy(&part[i].r, &part[i].l.r_ls, sizeof(ParticlePosition));
			memcpy(&part[i].m, &part[i].l.m_ls, sizeof(ParticleMomentum));
		}
	}
  //part = local_cells.cell[0]->part;
  //fprintf(stderr,"%d: load part %d: px_ls after %f, px after %f\n",this_node,part[0].p.identity,part[0].l.r_ls.p[0],part[0].r.p[0]);
  //fprintf(stderr,"%d: load part %d: mx_ls after %f, mx after %f\n",this_node,part[0].p.identity,part[0].l.m_ls.v[0],part[0].m.v[0]);
  
}


void hamiltonian_calc(int ekin_update_flag)
{
  
  /* if ekin_update_flag = 0, we calculate all energies with \ref master_energy_calc().
   if ekin_update_flag = 1, we only updated momenta, so there we only need to recalculate 
   kinetic energy with \ref calc_kinetic(). */
  
  int i;
  double result = 0.0;
  double ekt, ekr;
  
  INTEG_TRACE(fprintf(stderr,"%d: hamiltonian_calc:\n",this_node));
  
  //initialize energy struct 
  if (total_energy.init_status == 0) {
    init_energies(&total_energy);
    //if we are initializing energy we have to calculate all energies anyway
    ekin_update_flag = 0;
  }
  
  //calculate energies
  if (ekin_update_flag == 0)
    master_energy_calc();
  else
    calc_kinetic(&ekt, &ekr);

  //sum up energies on master node, and update ghmcdata struct
  if (this_node==0) {
    ghmcdata.hmlt_old = ghmcdata.hmlt_new;
    for (i = ekin_update_flag; i < total_energy.data.n; i++) {
      result += total_energy.data.e[i];
    }
    if (ekin_update_flag == 1)
      result += ekt+ekr;
    ghmcdata.hmlt_new=result;
  }

}

//get local temperature - here for debbuging purposes
double calc_local_temp()
{
	
	int i, j, c, np, tot_np = 0;
  Particle *part;
	
	double temp = 0.0;
	
	for (c = 0; c < local_cells.n; c++) {
    np   = local_cells.cell[c]->n;
		tot_np += np;
    part = local_cells.cell[c]->part;
    for (i = 0; i < np; i++) {
      for (j = 0; j < 3; j++) {
				temp += PMASS(part[i])*SQR(part[i].m.v[j]/time_step);
				#ifdef ROTATION
					temp += SQR(part[i].m.omega[j]);
				#endif	
			}
    }
	}
	#ifdef ROTATION
		temp /= 6*tot_np;
	#else
		temp /= 3*tot_np;
	#endif	
  return temp;
	
}

void calc_kinetic(double *ek_trans , double *ek_rot)
{
	
	int i, c, np;
  Particle *part;
	double et = 0.0, er = 0.0;
		
	for (c = 0; c < local_cells.n; c++) {
    np   = local_cells.cell[c]->n;
		
    part = local_cells.cell[c]->part;
    for (i = 0; i < np; i++) {
			#ifdef VIRTUAL_SITES
				if (ifParticleIsVirtual(&part[i])) continue;
			#endif
			
			/* kinetic energy */
			et += (SQR(part[i].m.v[0]) + SQR(part[i].m.v[1]) + SQR(part[i].m.v[2]))*PMASS(part[i]);

			/* rotational energy */
			#ifdef ROTATION
			#ifdef ROTATIONAL_INERTIA
				er += SQR(part[i].m.omega[0])*part[i].p.rinertia[0] +
								SQR(part[i].m.omega[1])*part[i].p.rinertia[1] +
								SQR(part[i].m.omega[2])*part[i].p.rinertia[2];
			#else
				er += SQR(part[i].m.omega[0]) + SQR(part[i].m.omega[1]) + SQR(part[i].m.omega[2]);
			#endif
			#endif	
    }
	}
	
  MPI_Allreduce(MPI_IN_PLACE, &et, 1, MPI_DOUBLE, MPI_SUM, comm_cart);
  MPI_Allreduce(MPI_IN_PLACE, &er, 1, MPI_DOUBLE, MPI_SUM, comm_cart);
	
	et /= (2.0*time_step*time_step);
	#ifdef ROTATION
		er /= 2.0;
	#endif
		
	*ek_trans=et;
	*ek_rot=er;

}


/* momentum update step of ghmc */
void ghmc_momentum_update()
{
	
	INTEG_TRACE(fprintf(stderr,"%d: ghmc_momentum_update:\n",this_node));
  
  /* currently, when temperature scaling is enabled the procedure tscale_momentum_update is used,
   although it may seem like a code duplication, this separation is maintained for now until this feature is tested thoroughly*/
  
  if (ghmc_tscale == GHMC_TSCALE_ON)
    tscale_momentum_update();
  else
    partial_momentum_update();
  
  int ekin_update_flag = 1;
  hamiltonian_calc(ekin_update_flag);
  
}

/* momentum update step of ghmc with temperature scaling */
void tscale_momentum_update()
{
  
  int i, j, c, np;
  Particle *part;
  
  //fprintf(stderr,"%d: temp before update: %f\n",this_node,calc_local_temp());

  save_last_state();
  
  simple_momentum_update();

  //fprintf(stderr,"%d: temp after simple update: %f\n",this_node,calc_local_temp());
  
  double tempt, tempr;
  calc_kinetic(&tempt, &tempr);
  tempt /= (1.5*n_part);
  tempr /= (1.5*n_part);
  
  double scalet = sqrt(temperature / tempt);
  #ifdef ROTATION    
    double scaler = sqrt(temperature / tempr);
  #endif  
            
  for (c = 0; c < local_cells.n; c++) {
    np   = local_cells.cell[c]->n;
    part = local_cells.cell[c]->part;
    for (i = 0; i < np; i++) {
      for (j = 0; j < 3; j++) {
        part[i].m.v[j] *= scalet;
        #ifdef ROTATION
          part[i].m.omega[j] *= scaler;
        #endif  
      }
    }
  }
  
  //fprintf(stderr,"%d: temp after scale: %f\n",this_node,calc_local_temp());
  
	for (c = 0; c < local_cells.n; c++) {
    np   = local_cells.cell[c]->n;
    part = local_cells.cell[c]->part;
    for (i = 0; i < np; i++) {
			#ifdef VIRTUAL_SITES
				if (ifParticleIsVirtual(&part[i])) continue;
			#endif
			
      for (j = 0; j < 3; j++) {
				part[i].m.v[j] = cosp*(part[i].l.m_ls.v[j])+sinp*(part[i].m.v[j]);
				#ifdef ROTATION
					part[i].m.omega[j] = cosp*(part[i].l.m_ls.omega[j])+sinp*(part[i].m.omega[j]);
				#endif	
			}
    }
	}
	
	//fprintf(stderr,"%d : temp after partial update: %f\n",this_node,calc_local_temp());
	
}

/* momentum update step of ghmc */
void simple_momentum_update()
{
	
	int i, j, c, np;
  Particle *part;
	double sigmat, sigmar;
	

	sigmat = sqrt(temperature); sigmar = sqrt(temperature);
	for (c = 0; c < local_cells.n; c++) {
    np   = local_cells.cell[c]->n;
    part = local_cells.cell[c]->part;
    for (i = 0; i < np; i++) {
			#ifdef VIRTUAL_SITES
				if (ifParticleIsVirtual(&part[i])) continue;
			#endif
			
			#ifdef MASS
				sigmat = sqrt(temperature / PMASS(part[i]));
			#endif
      for (j = 0; j < 3; j++) {
				part[i].m.v[j] = sigmat*gaussian_random()*time_step;
				#ifdef ROTATION
					#ifdef ROTATIONAL_INERTIA
						sigmar = sqrt(temperature / part[i].p.rinertia[j]);
					#endif
					part[i].m.omega[j] = sigmar*gaussian_random();
				#endif	
			}
    }
	}
	
	//fprintf(stderr,"%d: temp after simple update: %f\n",this_node,calc_local_temp());
	
}

/* partial momentum update step for ghmc */
void partial_momentum_update()
{

	int i, j, c, np;
  Particle *part;
	double sigmat, sigmar;

	//fprintf(stderr,"%d: temp before partial update: %f. expected: %f\n",this_node,calc_local_temp(),temperature);
	sigmat = sqrt(temperature); sigmar = sqrt(temperature);
	for (c = 0; c < local_cells.n; c++) {
    np   = local_cells.cell[c]->n;
    part = local_cells.cell[c]->part;
    for (i = 0; i < np; i++) {
			#ifdef MASS
				sigmat = sqrt(temperature / PMASS(part[i]));
			#endif
      for (j = 0; j < 3; j++) {
				part[i].m.v[j] = cosp*(part[i].m.v[j])+sinp*(sigmat*gaussian_random()*time_step);
				#ifdef ROTATION
					#ifdef ROTATIONAL_INERTIA
						sigmar = sqrt(temperature / part[i].p.rinertia[j]);
					#endif
					part[i].m.omega[j] = cosp*(part[i].m.omega[j])+sinp*(sigmar*gaussian_random());
				#endif	
			}
    }
	}
	//fprintf(stderr,"%d: temp after partial update: %f\n", this_node, calc_local_temp());
}

/* momentum flip for ghmc */
void momentum_flip()
{

	int i, j, c, np;
  Particle *part;

	INTEG_TRACE(fprintf(stderr,"%d: ghmc_momentum_flip:\n",this_node));
	
	for (c = 0; c < local_cells.n; c++) {
    np   = local_cells.cell[c]->n;
    part = local_cells.cell[c]->part;
    for (i = 0; i < np; i++) {
      for (j = 0; j < 3; j++) {
				part[i].m.v[j] = -(part[i].m.v[j]);
				#ifdef ROTATION
					part[i].m.omega[j] = -(part[i].m.omega[j]);
				#endif	
			}
    }
	}
	
}
/*@}*/

#endif


/************************************************************/
/** \name Exported Functions */
/************************************************************/
/*@{*/

#ifdef GHMC

void thermo_init_ghmc()
{
		ghmc_att=0;
		ghmc_acc=0;
    
		cosp = ghmc_phi;
    sinp = sin(acos(ghmc_phi));
		
		ghmcdata.hmlt_old=0;
		ghmcdata.hmlt_new=0;
		beta=1.0/temperature;
		
    THERMO_TRACE(fprintf(stderr,"%d: thermo_init_ghmc: ghmc_csp=%f, ghmc_snp=%f \n",this_node,cosp,sinp));
}

void ghmc_init()
{
		INTEG_TRACE(fprintf(stderr,"%d: ghmc_init:\n",this_node));
		
	  ghmcdata.att=0;
		ghmcdata.acc=0;

		save_last_state();
}

void ghmc_close()
{
	
		INTEG_TRACE(fprintf(stderr,"%d: ghmc_close:\n",this_node));
	  ghmc_att += ghmcdata.att;
		ghmc_acc += ghmcdata.acc;
		
}

/* monte carlo step of ghmc - evaluation stage */
void ghmc_mc()
{
			INTEG_TRACE(fprintf(stderr,"%d: ghmc_mc:\n",this_node));
			
		  double boltzmann;
			
      int ekin_update_flag = 0;
			hamiltonian_calc(ekin_update_flag);
			
      //make MC decision only on the master
      if (this_node==0) {
      
        ghmcdata.att++;
      
        //metropolis algorithm
        boltzmann = ghmcdata.hmlt_new - ghmcdata.hmlt_old;
        if (boltzmann < 0)
          boltzmann = 1.0;
        else if (boltzmann > 30)
          boltzmann = 0.0;
        else
          boltzmann = exp(-beta*boltzmann);
        
        //fprintf(stderr,"old hamiltonian : %f, new hamiltonian: % f, boltzmann factor: %f\n",ghmcdata.hmlt_old,ghmcdata.hmlt_new,boltzmann);

        if ( d_random() < boltzmann) {
          ghmcdata.acc++;
          ghmc_mc_res = GHMC_MOVE_ACCEPT;
        } else {
          ghmc_mc_res = GHMC_MOVE_REJECT;
        }
        
      }
      
      //let all nodes know about the MC decision result
      mpi_bcast_parameter(FIELD_GHMC_RES);
      
      if (ghmc_mc_res == GHMC_MOVE_ACCEPT) {
        save_last_state();
        //fprintf(stderr,"%d: mc move accepted\n",this_node);
      } else {
        load_last_state();
        //fprintf(stderr,"%d: mc move rejected\n",this_node);
        
        //if the move is rejected we might need to resort particles according to the loaded configurations
        cells_resort_particles(CELL_GLOBAL_EXCHANGE);
        invalidate_obs();
        
        if (ghmc_mflip == GHMC_MFLIP_ON) {
          momentum_flip();
        } else if (ghmc_mflip == GHMC_MFLIP_RAND) {
          if (d_random() < 0.5) momentum_flip();
        }
      }      
        
      //fprintf(stderr,"%d: temp after mc: %f\n",this_node,calc_local_temp());
}

/*@}*/

#endif