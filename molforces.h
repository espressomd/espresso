// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.

#include <mpi.h>
#include "utils.h"
#include "grid.h"

/** \file molforces.h
 *  Routines for calculating and applying trap forces upon molecules.
 *  This trap force can be set to
 *  - a harmonic potential with a restlength of zero on the molecular centre of mass
 *  - a drag on the molecular velocity
 *  - a cancelation of the total force on the molecule (including thermostat forces)
 *  The centre of mass can be fixed to an absolute position or to a relative position in the
 *  simulation box.
 *  The molecular trap forces is distributed evenly upon all particles in a molecule.
 *  (see \file topology.c and file \molforces.c)  
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:cooke@mpip-mainz.mpg.de">Ira</a>
 *  <a href="mailto:reynolds@mpip-mainz.mpg.de">Ben</a>
 */

#ifdef MOLFORCES

extern int     IsTrapped;

/**
   Checks if there are any molecules trapped (IsTrapped=1) and if so calls calc_mol_info
   and apply_mol_constaints */
void calc_and_apply_mol_constraints();

#endif
