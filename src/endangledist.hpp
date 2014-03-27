/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
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

/** \file endangledist.hpp
 *  Routines which apply an angle potential between two particles and a wall constraint
 *  At distmax the angle potential is slowly switched on to a maximum at distmin
 *  phi0 is constant but could easily be implemented to depend on the distance 
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:johnston@mpip-mainz.mpg.de">Karen</a>
*/

/************************************************************/
#ifndef _ENDANGLEDIST_H
#define _ENDANGLEDIST_H

#include "particle_data.hpp"
#include "interaction_data.hpp"

#ifdef BOND_ENDANGLEDIST

/// set parameters for endangledist potential
int endangledist_set_params(int bond_type, double bend, double phi0 ,double distmin, double distmax);

///
int calc_endangledist_pair_force(Particle *p1, Particle *p2,
				 Bonded_ia_parameters *iaparams,
				 double dx[3], double force1[3],
				 double force2[3]);

///
int endangledist_pair_energy(Particle *p1, Particle *p2,
			     Bonded_ia_parameters *iaparams,
			     double dx[3], double *_energy);

#endif /* BOND_ENDANGLEDIST */
#endif /* ENDANGLEDIST_H */
