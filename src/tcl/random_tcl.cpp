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
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "utils.hpp"
#include "random.hpp"
#include "parser.hpp"
#include "tcl.h"
#include "communication.hpp"


/*----------------------------------------------------------------------*/

/**  Implementation of the tcl-command
     t_random [{ int \<n\> | seed [\<seed(0)\> ... \<seed(n_nodes-1)\>] | stat [status-list] }]
     <ul>
     <li> Without further arguments, it returns a random double between 0 and 1.
     <li> If 'int \<n\>' is given, it returns a random integer between 0 and n-1.
     <li> If 'seed'/'stat' is given without further arguments, it returns a tcl-list with
          the current seeds/status of the n_nodes active nodes; otherwise it issues the 
	  given parameters as the new seeds/status to the respective nodes.     
     </ul>
 */
int tclcommand_t_random (ClientData data, Tcl_Interp *interp, int argc, char **argv) {
  char buffer[100 + TCL_DOUBLE_SPACE + 3*TCL_INTEGER_SPACE];
  int i,j,cnt, i_out, temp; 
  double d_out;

  if (argc == 1) {                  /* 't_random' */
    d_out = d_random();
    sprintf(buffer, "%f", d_out); Tcl_AppendResult(interp, buffer, (char *) NULL); return (TCL_OK);
  }

  /* argc > 1 */
  argc--; argv++;

  if ( ARG_IS_S(0,"int") )          /* 't_random int <n>' */
  {
    if(argc < 2)
    { 
      Tcl_AppendResult(interp, "\nWrong # of args: Usage: 't_random int <n>'", (char *) NULL); 
      return (TCL_ERROR); 
    }
    else 
    {
      if( !ARG_IS_I(1,i_out) )
      { 
        Tcl_AppendResult(interp, "\nWrong type: Usage: 't_random int <n>'", (char *) NULL); 
        return (TCL_ERROR); 
      }

      i_out = i_random(i_out);

      sprintf(buffer, "%d", i_out);
      Tcl_AppendResult(interp, buffer, (char *) NULL); 
      return (TCL_OK);
    }
  }
  else if ( ARG_IS_S(0,"seed") )    /* 't_random seed [<seed(0)> ... <seed(n_nodes-1)>]' */
  {
    int *seed = (int *) Utils::malloc(n_nodes*sizeof(int));

    if (argc <= 1)                  /* ESPResSo generates a seed */
    {
      mpi_random_seed(0,seed);

      for (i=0; i < n_nodes; i++) 
      {
	      sprintf(buffer, "%d ", seed[i]); 
        Tcl_AppendResult(interp, buffer, (char *) NULL); 
      }
    }
    else if (argc < n_nodes+1)      /* Fewer seeds than nodes */
    { 
      sprintf(buffer, "Wrong # of args (%d)! Usage: 't_random seed [<seed(0)> ... <seed(%d)>]'", argc,n_nodes-1);
      Tcl_AppendResult(interp, buffer, (char *)NULL); 
      return (TCL_ERROR); 
    }
    else                            /* Get seeds for different nodes */
    {
      for (i=0; i < n_nodes; i++) 
      {
        if( !ARG_IS_I(i+1,temp) )
        { 
          sprintf(buffer, "\nWrong type for seed %d:\nUsage: 't_random seed [<seed(0)> ... <seed(%d)>]'", i+1 ,n_nodes-1);
          Tcl_AppendResult(interp, buffer, (char *)NULL);  
          return (TCL_ERROR); 
        }
        else
        {
          seed[i] = (int)temp;
        }
      }

      RANDOM_TRACE(
        printf("Got "); 

        for(i=0;i<n_nodes;i++) 
          printf("%ld ",seed[i]);

        printf("as new seeds.\n")
      );

      mpi_random_seed(n_nodes,seed);
    }

    free(seed); 
    return(TCL_OK);
  }

  /* else */
  sprintf(buffer, "Usage: 't_random [{ int <n> | seed [<seed(0)> ... <seed(%d)>] | stat [status-list] }]'",n_nodes-1);
  Tcl_AppendResult(interp, "Unknown job '",argv[0],"' requested!\n",buffer, (char *)NULL);
  return (TCL_ERROR); 
}

