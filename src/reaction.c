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

#include "communication.h"
#include "utils.h"
#include "interaction_data.h"
#include "reaction.h"

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
#endif /* ifdef REACTIONS */
