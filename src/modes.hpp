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
#ifndef MODES_H
#define MODES_H

/** \file modes.hpp

    PLEASE INSERT DESCRIPTION
*/

#include "utils.hpp"
#include "statistics.hpp"

#ifdef FFTW

#ifdef _Complex_I
#warning the complex data type is predefined on your system, hoping it is compatible to a double[2]
#endif
#include <fftw3.h>
#define FFTW_REAL(x) (((double *)(x))[0])
#define FFTW_IMAG(x) (((double *)(x))[1])

/** The full 3d grid for mode analysis */
extern int mode_grid_3d[3];
/** Integer labels for grid axes compared to real axes*/
extern int xdir;
extern int ydir;
extern int zdir;

/** Enumerated constant indicating a Lipid in the top leaflet*/
#define LIPID_UP 0
/** Enumerated constant indicating a Lipid in the bottom leaflet*/
#define LIPID_DOWN 1
/** Enumerated constant indicating a Lipid that has left the bilayer
    but may have become incorporated into a periodic image bilayer */
#define LIPID_STRAY 2
/** Enumerated constant indicating a Lipid that has left the bilayer
    and truly floating in space */
#define REAL_LIPID_STRAY 3
/** The atom type corresponding to a lipid head group */
#define LIPID_HEAD_TYPE 0

/** Flag to indicate when the mode_grid is changed */
extern int mode_grid_changed;

/** Parameter indicating distance beyond which a lipid is said to have
    left the membrane the default value is set in \ref modes.cpp */
extern double stray_cut_off;

/* Exported Functions */
/* switch_fluc == 1 comuptes height grid
               == 0 thickness.
*/
int modes2d(fftw_complex* result, int switch_fluc );
void map_to_2dgrid();
/** 
    This routine performs a simple check to see whether a lipid is
    oriented up or down or if it has escaped the bilayer.  In order
    for this routine to work it is essential that the lipids head
    groups are of atom type LIPID_HEAD_TYPE

    \param id The particle identifier
    \param partCfg An array of sorted particles

    \param zref The average z position of all particles. This is used
    to check for stray lipids

    \param director director
    \param refdir is a vector indicating the direction indicating
    up. This is usually the zaxis. If it is not the z axis then lipids
    will not be returned as stray.
 */
int lipid_orientation( int id, Particle* partCfg , double zref, double director[3],double refdir[3]);

/**
   This routine calculates the orientational order parameter for a
   lipid bilayer.  It also calculates the direction of orientation for
   all lipids and places them in stored_dirs.
*/
int orient_order(double* result, double* stored_dirs);

/* Get the list of lipid orientations */
int get_lipid_orients(IntList* l_orient);

/**
   Calculate the average height of lipids or average thickness for each value in a grid
   switch_fluc == 1: height_grid
   switch_fluc == 0: thickness
   The output is written in *height_grid no matter what the argument is.
*/
int calc_fluctuations ( double* height_grid, int switch_fluc );

/** 
    Calculate a  vertical density profile for each of the specified beadtypes

    \param beadids The list of bead types for which profiles will be generated
    \param hrange The vertical range from the bilayer midplane over which the profiles are calculated
    \param density_profile The pre-allocated density profile into which data is written
    \param usegrid switch to determine whether grid should be used
 */
int bilayer_density_profile ( IntList *beadids, double hrange , DoubleList *density_profile, int usegrid );
int bilayer_density_profile_sphere (IntList *beadids, double rrange , DoubleList *density_profile, double radius, double center[3]);


#endif

#endif
