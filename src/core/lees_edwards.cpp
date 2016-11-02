/*
  Copyright (C) 2010,2011,2012,2016 The ESPResSo project
  Copyright (C) 2009,2010 
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
#include "lees_edwards.hpp"
#include "domain_decomposition.hpp"
#include "lees_edwards_domain_decomposition.hpp"
#include "grid.hpp"

/* global state variables */
double lees_edwards_offset      = 0.0;
double lees_edwards_rate        = 0.0;

#ifdef LEES_EDWARDS
int lees_edwards_count          =   0;

/* local state variables */
double lees_edwards_prev_set_at = 0.0;
double lees_edwards_prev_offset = 0.0;


/** React to news of a broadcast change in the \ref lees_edwards_offset */
void lees_edwards_step_boundaries(){

 double delta;

 delta = lees_edwards_offset - lees_edwards_prev_offset;
 while( delta >  0.5 * box_l[0] ) delta -= box_l[0];
 while( delta < -0.5 * box_l[0] ) delta += box_l[0];
 
 while( lees_edwards_offset >  0.5*box_l[0]) { lees_edwards_count++; lees_edwards_offset -= box_l[0];}
 while( lees_edwards_offset < -0.5*box_l[0]) { lees_edwards_count--; lees_edwards_offset += box_l[0];}
   
 /* update the apparent shear rate */
 if( sim_time - lees_edwards_prev_set_at > 0.0 )
    lees_edwards_rate     = delta * time_step / (sim_time - lees_edwards_prev_set_at);
 else
    lees_edwards_rate     = 0.0;

 /* save some state */
 lees_edwards_prev_set_at = sim_time;
 lees_edwards_prev_offset = lees_edwards_offset;

 /* request a new verlet list */
 rebuild_verletlist    = 1;
 
 /* request a redo of particle-cell assignment */
 resort_particles      = 1;
 
 /* Only part of the comms system needs to be rebuilt,
    but it is still very slow. */
 cells_on_geometry_change( CELL_FLAG_LEES_EDWARDS );

 return;
}

#endif
