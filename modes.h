// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2005; all rights reserved unless otherwise stated.
#ifndef MODES_H
#define MODES_H

/** \file modes.h

    PLEASE INSERT DESCRIPTION

    <b>Responsible:</b>
    <a href="mailto:cooke@mpip-mainz.mpg.de">Ira Cooke</a>
*/

#include "utils.h"
#include "statistics.h"
#include "parser.h"

#ifdef USEFFTW3
#ifdef _Complex_I
#  warning the complex data type is predefined on your system, hoping it is compatible to a double[2]
#endif
#  include <fftw3.h>
#  define FFTW_REAL(x) (((double *)(x))[0])
#  define FFTW_IMAG(x) (((double *)(x))[1])
#else
#  include <rfftw.h>
#  define FFTW_REAL(x) ((x).re)
#  define FFTW_IMAG(x) ((x).im)
#endif

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
    left the membrane the default value is set in \ref modes.c */
extern double stray_cut_off;

/* Exported Functions */
int modes2d(fftw_complex* result);
void map_to_2dgrid();
/** 
    This routine performs a simple check to see whether a lipid is
    oriented up or down or if it has escaped the bilayer.  In order
    for this routine to work it is essential that the lipids head
    groups are of atom type LIPID_HEAD_TYPE

    \param id The particle identifier
    \param partCfg An array of sorted particles
    \param zref The average z position of all particles


 */
int lipid_orientation( int id, Particle* partCfg , double zref, double director[3]);

/**
   This routine calculates the orientational order parameter for a
   lipid bilayer.  It also calculates the direction of orientation for
   all lipids and places them in stored_dirs.
*/
int orient_order(double* result, double* stored_dirs);

/* Get the list of lipid orientations */
int get_lipid_orients(IntList* l_orient);

/**
   Calculate the average height of lipids for each value in a grid
*/
int calc_height_grid ( double* height_grid );

/** 
    Calculate a  vertical density profile for each of the specified beadtypes

    \param beadids The list of bead types for which profiles will be generated
    \param hrange The vertical range from the bilayer midplane over which the profiles are calculated
    \param density_profile The pre-allocated density profile into which data is written
 */
int bilayer_density_profile ( IntList *beadids, double hrange , DoubleList *density_profile, int usegrid );

#endif



