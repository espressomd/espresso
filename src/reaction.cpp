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
/** \file reaction.cpp
 *
 */

#include "utils.hpp"
#include "reaction.hpp"
#include "initialize.hpp"
#include "forces.hpp"
#include "errorhandling.hpp"

reaction_struct reaction;

#ifdef CATALYTIC_REACTIONS

void reactions_sanity_checks()
{
  char *errtext;

  if(reaction.ct_rate != 0.0) {

    if( dd.use_vList == 0 || cell_structure.type != CELL_STRUCTURE_DOMDEC) {
      errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"{105 The CATALYTIC_REACTIONS feature requires verlet lists and domain decomposition} ");
    }

    if(max_cut < reaction.range) {
      errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"{106 Reaction range of %f exceeds maximum cutoff of %f} ", reaction.range, max_cut);
    }
  }
}


void local_setup_reaction() {
  
  /* Make available the various reaction parameters */
  MPI_Bcast(&reaction.reactant_type, 1, MPI_INT, 0, comm_cart);
  MPI_Bcast(&reaction.product_type, 1, MPI_INT, 0, comm_cart);
  MPI_Bcast(&reaction.catalyzer_type, 1, MPI_INT, 0, comm_cart);
  MPI_Bcast(&reaction.range, 1, MPI_DOUBLE, 0, comm_cart);
  MPI_Bcast(&reaction.ct_rate, 1, MPI_DOUBLE, 0, comm_cart);
  MPI_Bcast(&reaction.eq_rate, 1, MPI_DOUBLE, 0, comm_cart);
  MPI_Bcast(&reaction.sing_mult, 1, MPI_INT, 0, comm_cart);

  /* Create the various reaction related types (categories) */
  make_particle_type_exist(reaction.catalyzer_type);
  make_particle_type_exist(reaction.product_type);
  make_particle_type_exist(reaction.reactant_type);

  /* Make ESPResSo aware that reactants and catalyst are interacting species */
  IA_parameters *data = get_ia_param_safe(reaction.reactant_type, reaction.catalyzer_type);
  
  if(!data) {    
	  char *error_msg = runtime_error(128);
	  ERROR_SPRINTF(error_msg, "{106 interaction parameters for reaction could not be set} ");
  }

  /* Used for the range of the verlet lists */
  data->REACTION_range = reaction.range;

  /* Broadcast interaction parameters */
  mpi_bcast_ia_params(reaction.reactant_type, reaction.catalyzer_type);
}

void integrate_reaction() {
  int c, np, n, i,
      check_catalyzer;
  Particle *p1, *p2, **pairs;
  Cell *cell;
  double dist2, vec21[3],
         ct_ratexp, eq_ratexp,
         rand, bernoulli;

  if(reaction.ct_rate > 0.0) {

    /* Determine the reaction rate */
    ct_ratexp = exp(-time_step*reaction.ct_rate);
    
    on_observable_calc();

    for (c = 0; c < local_cells.n; c++) {

      /* Take into account only those cell neighbourhoods for which
         the central cell contains a catalyzer particle */

      check_catalyzer = 0;

      cell = local_cells.cell[c];
      p1   = cell->part;
      np  = cell->n;
      
      for(i = 0; i < np; i++) {
        if(p1[i].p.type == reaction.catalyzer_type) {
          check_catalyzer = 1;
          break;
        }
      }	

      /* If the central cell contains a catalyzer particle, ...*/
      if ( check_catalyzer != 0 ) {

        /* Loop cell neighbors */
        for (n = 0; n < dd.cell_inter[c].n_neighbors; n++) {
          pairs = dd.cell_inter[c].nList[n].vList.pair;
          np = dd.cell_inter[c].nList[n].vList.n;

          /* Verlet list loop */
          for(i = 0; i < 2 * np; i += 2) {
            p1 = pairs[i];   //pointer to particle 1
            p2 = pairs[i+1]; //pointer to particle 2

            if( (p1->p.type == reaction.reactant_type &&  p2->p.type == reaction.catalyzer_type) || (p2->p.type == reaction.reactant_type &&  p1->p.type == reaction.catalyzer_type) ) {
              get_mi_vector(vec21, p1->r.p, p2->r.p);
              dist2 = sqrlen(vec21);
              
              /* Count the number of times a reactant particle is
                 checked against a catalyst */
              if(dist2 < reaction.range * reaction.range) {

                if(p1->p.type == reaction.reactant_type) {
						       p1->p.catalyzer_count++;
					      }
					      else {
						       p2->p.catalyzer_count++;
					      }
             	}
            }  
          }
        }
      }
    }

    /* Now carry out the reaction on the particles which are tagged */
    for (c = 0; c < local_cells.n; c++) {
      cell = local_cells.cell[c];
      p1   = cell->part;
      np  = cell->n;
        
      for(i = 0; i < np; i++) {
        if(p1[i].p.type == reaction.reactant_type){

          if(p1[i].p.catalyzer_count > 0 ){

            if( reaction.sing_mult == 0 ) 
            {
              rand = d_random();

              bernoulli = pow(ct_ratexp, p1[i].p.catalyzer_count);

              if(rand > bernoulli) {
                p1[i].p.type = reaction.product_type; }
            
            }
            else /* We only consider each reactant once */
            {
              rand = d_random();

              if(rand > ct_ratexp) {
                p1[i].p.type = reaction.product_type; }
            }
            
            p1[i].p.catalyzer_count = 0;
          }
        }
      }
    }

    /* We only need to do something when the equilibrium 
       reaction rate constant is nonzero */
    if (reaction.eq_rate > 0.0) { 

  	  eq_ratexp = exp(-time_step*reaction.eq_rate);

      for (c = 0; c < local_cells.n; c++) {
        cell = local_cells.cell[c];
        p1   = cell->part;
        np  = cell->n;
    
        /* Convert products into reactants and vice versa
           according to the specified rate constant */
        for(i = 0; i < np; i++) {

          if(p1[i].p.type == reaction.product_type) {
            rand = d_random();
            
  	        if(rand > eq_ratexp) {
  	          p1[i].p.type=reaction.reactant_type;
  	        }
          }
          else if(p1[i].p.type == reaction.reactant_type) {
            rand = d_random();
            
  	        if(rand > eq_ratexp) {
  	          p1[i].p.type=reaction.product_type;
  	        }
          }

        }
      }
    }

    on_particle_change();
  }
}
#endif
