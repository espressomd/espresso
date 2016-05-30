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

#include "mpiio_tcl.hpp"
#include "communication.hpp"

#define LEN(x) (sizeof(x) / sizeof(*(x)))

/** Array describing the mapping of possible output descriptors to
 *  internal output specifiers.
 */
static const struct {
  const char *name;
  const MPIIOOutputFields val;
} out_desc[] = {
  { "pos", MPIIO_OUT_POS },
  { "v", MPIIO_OUT_VEL },
  { "type", MPIIO_OUT_TYP },
  { "bond", MPIIO_OUT_BND },
};


int tclcommand_mpiio(ClientData data, Tcl_Interp *interp,
                     int argc, char *argv[])
{
  char *filename;
  char *mode;
  unsigned output_fields = 0;
  bool correct_field;
  int write;

#ifdef MULTI_TIMESTEP
  Tcl_AppendResult(interp, "MPI-IO does not support multi-timestep velocity rescaling.",
                   (char *) NULL);
  return (TCL_ERROR);
#endif
  
  if (argc < 3) {
    Tcl_AppendResult(interp, "wrong # args:  should be \"",
                     argv[0], " <filename> read|write ?pos|v|types|bonds?* ...\"",
                     (char *) NULL);
    return (TCL_ERROR);
  }
  filename = argv[1];
  // Operation mode
  if (!strncmp(argv[2], "write", strlen(argv[2]))) {
    write = 1;
  } else if (!strncmp(argv[2], "read", strlen(argv[2]))) {
    write = 0;
  } else {
    Tcl_AppendResult(interp, "wrong operation mode: \"", argv[2],
                     "\"  should be \"read\" or \"write\"",
                     (char *) NULL);
    return (TCL_ERROR);
  }
  
  argv += 3;
  argc -= 3;

  // The remaining arguments are field output specifiers.
  // Parse them according to out_desc.
  for (; argc > 0; argc--, argv++) {
    correct_field = false;
    for (size_t i = 0; i < LEN(out_desc); ++i) {
      if (!strncmp(*argv, out_desc[i].name, strlen(*argv))) {
        output_fields |= out_desc[i].val;
        correct_field = true;
        break;
      }
    }
    if (!correct_field) {
      Tcl_AppendResult(interp, "unknown output field: \"", *argv,
                       "\"  see the user's guide for all possibilities",
                       (char *) NULL);
      return (TCL_ERROR);
    }
  }

  
  mpi_mpiio(filename, output_fields, write);

  return (TCL_OK);
}
