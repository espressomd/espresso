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
/** \file errorhandling.hpp 
    This file contains the errorhandling code for severe errors, like
    a broken bond or illegal parameter combinations. See section
    "Errorhandling for developers" for details on the error format and
    how to use this.
*/
#ifndef _ERRORHANDLING_HPP
#define _ERRORHANDLING_HPP

/** buffer for error messages during the integration process. */
extern char *error_msg;
extern int n_error_msg;

/* request space for leaving an error message to be passed to the master node.
   Also takes care of the error counter.
   @param errlen maximal length of the error message. If you use sprintf to create the error
   message, remember to use ES_INTEGER/DOUBLE_SPACE as usual
   @return where to put the (null-terminated) string */
char *runtime_error(int errlen);

#define ERROR_SPRINTF sprintf

/** check for runtime errors on all nodes. This has to be called on all nodes synchronously.
    @return the number of characters in the error messages of all nodes together. */
int check_runtime_errors();

/** exit ungracefully, core dump if switched on. */
void errexit();

/** register a handler for sigint that translates it into an background error. */
void register_sigint_handler();

#endif
