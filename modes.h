// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
#ifndef MODES_H
#define MODES_H

/** \file modes.h

    PLEASE INSERT DESCRIPTION

    <b>Responsible:</b>
    <a href="mailto:cooke@mpip-mainz.mpg.de">Ira Cooke</a>
*/

#include "statistics.h"
#include "parser.h"
#include "debug.h"
#include "utils.h"
//#include <fftw.h>
#include <rfftw.h>

/** The full 3d grid for mode analysis */
extern int mode_grid_3d[3];
/** Integer labels for grid axes compared to real axes*/
extern int xdir;
extern int ydir;
extern int zdir;
/** Flag to indicate when the mode_grid is changed */
extern int mode_grid_changed;

/** Parameter indicating distance beyond which a lipid is said to have
    left the membrane the default value is set in \ref modes.c */
extern double stray_cut_off;

/* Exported Functions */
int modes2d(fftw_complex* result);
void map_to_2dgrid();

#endif


