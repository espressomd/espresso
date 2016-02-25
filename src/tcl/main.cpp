/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
/** \file main.cpp
    Main file of Espresso. Initialization of tcl interpreter.
*/
/* first, since we need the TK define */
#include "utils.hpp"
#ifdef TK
#include <tk.h>
#endif
#include "initialize_interpreter.hpp"
#include "initialize.hpp"
#include "communication.hpp"

int main(int argc, char **argv)
{
  /* first thing to do: fire up MPI */
  mpi_init(&argc, &argv);

  on_program_start();

  if (this_node == 0) {
    /* master node */
#ifdef TK
    Tk_Main(argc, argv, tcl_appinit);
#else
    Tcl_Main(argc, argv, tcl_appinit);
#endif
  }
  else {
    /* slave node */
    mpi_loop();
  }

  return 0;
}
