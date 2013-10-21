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
/** \file npt.hpp
    exports for the NPT code, which otherwise is really spread all over...
*/

#ifndef NPT_H
#define NPT_H

#include "utils.hpp"

/************************************************
 * data types
 ************************************************/

/** Structure to hold all variables related to the isotropic NpT-integration scheme. */
typedef struct { 
  /** mass of a virtual piston representing the shaken box */
  double piston;
  /** inverse of piston */
  double inv_piston;
  /** isotropic volume.  Note that we use the term volume throughout
      although for a 2d or 1d system we mean Area and Length
      respectively */
  double volume;

  /** desired pressure to which the algorithm strives to */
  double p_ext;
  /** instantaneous pressure the system currently has */
  double p_inst;
  /** instantaneous pressure averaged over current integration cycle */
  double p_inst_av;
  /** difference between \ref p_ext and \ref p_inst */
  double p_diff;
  /** virial (short-range) components of \ref p_inst */
  double p_vir[3];
  /** ideal gas components of \ref p_inst, derived from the velocities */
  double p_vel[3];
  /** flag which indicates if \ref p_vel may (0) or may not (1) be used
      in offline pressure calculations such as 'analyze p_inst' */ 
  int invalidate_p_vel;
  /** geometry information for the npt integrator.  Holds the vector
      <dir, dir ,dir> where a positive value for dir indicates that
      box movement is allowed in that direction. To check whether a
      given direction is turned on use bitwise comparison with \ref
      nptgeom_dir */
  int geometry;
  /** bitwise comparison values corresponding to different directions*/
  int nptgeom_dir[3];
  /** The number of dimensions in which npt boxlength motion is coupled to particles */
  int dimension;
  /** Set this flag if you want all box dimensions to be identical. Needed for electrostatics and magnetostatics.  
      If the value of dimension is less than 3 then box length motion in one or more
      directions will be decoupled from the particle motion */
  int cubic_box;
  /** An index to one of the non_constant dimensions. handy if you just want the variable box_l */
  int non_const_dim;
} nptiso_struct;
extern nptiso_struct nptiso;

/** Allowable values for nptiso.geometry*/
#define NPTGEOM_XDIR 1
#define NPTGEOM_YDIR 2
#define NPTGEOM_ZDIR 4

#endif
