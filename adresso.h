// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.

#ifndef ADRESSO_H
#define ADRESSO_H
/** \file adresso.h
    This is the place for adaptive resolution scheme (adress)
    Implementation of adresso.h
    <b>Responsible:</b>
    <a href="mailto:junghans@mpip-mainz.mpg.de">Axel</a>
*/

#include <tcl.h>
#include "particle_data.h"

/** \name Exported Variables */
/************************************************************/
/*@{*/
extern double adress_vars[7];
/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/
/** Implements the Tcl command \ref tcl_adress. This allows for seetings for adress
*/
int adress_tcl(ClientData data, Tcl_Interp *interp, int argc, char **argv);

#ifdef ADRESS
/** Calc adress weight function of a vector
    @param x[3] vector
    @return weight of the vector
*/
double adress_wf_vector(double x[3]);

/** Calc adress weight function of a particle
    @param x[3] vector
    @return weight of the particle
*/
MDINLINE double adress_wf_particle(Particle *p){
   return adress_wf_vector(p->r.p);
}
#endif
/*@}*/

#endif
