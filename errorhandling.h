// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
/** \file errorhandling.h
    This file contains the errorhandling code for severe errors, like a broken bond or illegal parameter
    combinations. See section \ref errors "Errorhandling" for details on the error format and
    \ref errorhandling "Errorhandling for developers" for implementation hints.

    <b>Responsible:</b>
    <a href="mailto:arnolda@mpip-mainz.mpg.de">Axel</a>
*/
#ifndef ERRORHANDLING_H
#define ERRORHANDLING_H

#include <tcl.h>

/** buffer for error messages during the integration process. */
extern char *error_msg;
extern int n_error_msg;


/* request space for leaving an error message to be passed to the master node.
   Also takes care of the error counter.
   @param errlen maximal length of the error message. If you use sprintf to create the error
   message, remember to use TCL_INTEGER/DOUBLE_SPACE as usual
   @return where to put the (null-terminated) string */
char *runtime_error(int errlen);

/** check for runtime errors on all nodes. This has to be called on all nodes synchronously.
    @return the number of characters in the error messages of all nodes together. */
int check_runtime_errors();

#endif
