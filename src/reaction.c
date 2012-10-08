/*
  Copyright (C) 2010,2011,2012 The ESPResSo project
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
/** \file reaction.c
 *
 */

#include "utils.h"
#include "reaction.h"
#include "initialize.h"
#include "forces.h"


#ifdef REACTIONS
void setup_reaction() {
  MPI_Bcast(&reaction.reactant_type, 1, MPI_INT, 0, comm_cart);
  MPI_Bcast(&reaction.product_type, 1, MPI_INT, 0, comm_cart);
  MPI_Bcast(&reaction.catalyzer_type, 1, MPI_INT, 0, comm_cart);
  MPI_Bcast(&reaction.range, 1, MPI_DOUBLE, 0, comm_cart);
  MPI_Bcast(&reaction.rate, 1, MPI_DOUBLE, 0, comm_cart);
  MPI_Bcast(&reaction.back_rate, 1, MPI_DOUBLE, 0, comm_cart);
    
  make_particle_type_exist(reaction.catalyzer_type);
  make_particle_type_exist(reaction.product_type);
  make_particle_type_exist(reaction.reactant_type);

  IA_parameters *data = get_ia_param_safe(reaction.reactant_type, reaction.catalyzer_type);
  
  if(!data) {    
	  char *error_msg = runtime_error(128);
	  ERROR_SPRINTF(error_msg, "{106 interaction parameters for reaction could not be set} ");
  }
  
  data->REACTION_range = reaction.range;

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(reaction.reactant_type, reaction.catalyzer_type);
}

void integrate_reaction() {
  int c, np, n, i;
  Particle * p1, *p2, **pairs;
  Cell * cell;
  double dist2, vec21[3], rand;

  if(reaction.rate > 0) {
    int reactants = 0, products = 0;
    int tot_reactants = 0, tot_products = 0;
    double ratexp, back_ratexp;

    ratexp = exp(-time_step*reaction.rate);
    
    on_observable_calc();

    for (c = 0; c < local_cells.n; c++) {
      /* Loop cell neighbors */
      for (n = 0; n < dd.cell_inter[c].n_neighbors; n++) {
        pairs = dd.cell_inter[c].nList[n].vList.pair;
        np = dd.cell_inter[c].nList[n].vList.n;
        
        /* verlet list loop */
        for(i = 0; i < 2 * np; i += 2) {
          p1 = pairs[i];   //pointer to particle 1
          p2 = pairs[i+1]; //pointer to particle 2
          
          if( (p1->p.type == reaction.reactant_type &&  p2->p.type == reaction.catalyzer_type) || (p2->p.type == reaction.reactant_type &&  p1->p.type == reaction.catalyzer_type) ) {
            get_mi_vector(vec21, p1->r.p, p2->r.p);
            dist2 = sqrlen(vec21);
            
            if(dist2 < reaction.range * reaction.range) {
           	  rand = d_random();
           	  
           		if(rand > ratexp) {
printf("DEBUG: integrate_reaction; react\n"); //TODO delete
           		  if(p1->p.type == reaction.reactant_type) {
						      p1->p.type = reaction.product_type;
					      }
					      else {
						      p2->p.type = reaction.product_type;
					      }
              }
           	}
          }  
        }
      }
    }

    if (reaction.back_rate < 0) { // we have to determine it dynamically 
      /* we count now how many reactants and products are in the sim box */
      for (c = 0; c < local_cells.n; c++) {
        cell = local_cells.cell[c];
        p1   = cell->part;
        np  = cell->n;
        
        for(i = 0; i < np; i++) {
          if(p1[i].p.type == reaction.reactant_type)
            reactants++;
          else if(p1[i].p.type == reaction.product_type)
            products++;
        }	
      }
      
      MPI_Allreduce(&reactants, &tot_reactants, 1, MPI_INT, MPI_SUM, comm_cart);
      MPI_Allreduce(&products, &tot_products, 1, MPI_INT, MPI_SUM, comm_cart);

      back_ratexp = ratexp * tot_reactants / tot_products ; //with this the asymptotic ratio reactant/product becomes 1/1 and the catalyzer volume only determines the time that it takes to reach this
    }
    else { //use the back reaction rate supplied by the user
  	  back_ratexp = exp(-time_step*reaction.back_rate);
    }
    
    if(back_ratexp < 1 ) {
      for (c = 0; c < local_cells.n; c++) {
        cell = local_cells.cell[c];
        p1   = cell->part;
        np  = cell->n;
        
        for(i = 0; i < np; i++) {
          if(p1[i].p.type == reaction.product_type) {
            rand = d_random();
            
			      if(rand > back_ratexp) {
printf("DEBUG: integrate_reaction; back_react\n"); //TODO delete
			        p1[i].p.type=reaction.reactant_type;
			      }
          }
        }
      }
    }

    on_particle_change();
  }
}
#endif /* ifdef REACTIONS */
