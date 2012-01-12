/*
  Copyright (C) 2010,2011 The ESPResSo project
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
//#include "global.h"
#include "../random.h"
#include "tcl.h"
#include "communication.h"


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
  int i,j,cnt, i_out; double d_out;

  if (argc == 1) {                                          /* 't_random' */
    d_out = d_random();
    sprintf(buffer, "%f", d_out); Tcl_AppendResult(interp, buffer, (char *) NULL); return (TCL_OK);
  }
  /* argc > 1 */
  argc--; argv++;
  if (!strncmp(argv[0], "int", strlen(argv[0]))) {        /* 't_random int <n>' */
    if(argc < 2){ Tcl_AppendResult(interp, "Wrong # of args: Usage: 't_random int <n>'", (char *) NULL); return (TCL_ERROR); }
    else {
      Tcl_GetInt(interp, argv[1], &i_out);
      i_out = i_random(i_out);
      sprintf(buffer, "%d", i_out); Tcl_AppendResult(interp, buffer, (char *) NULL); return (TCL_OK);
    }
  }
  else if (!strncmp(argv[0], "seed", strlen(argv[0]))) {  /* 't_random seed [<seed(0)> ... <seed(n_nodes-1)>]' */
    long *seed = malloc(n_nodes*sizeof(long));
    if (argc <= 1) {
      mpi_random_seed(0,seed);
      for (i=0; i < n_nodes; i++) { 
	sprintf(buffer, "%ld ", seed[i]); Tcl_AppendResult(interp, buffer, (char *) NULL); 
      }
    }
    else if (argc < n_nodes+1) { 
      sprintf(buffer, "Wrong # of args (%d)! Usage: 't_random seed [<seed(0)> ... <seed(%d)>]'", argc,n_nodes-1);
      Tcl_AppendResult(interp, buffer, (char *)NULL); return (TCL_ERROR); }
    else {
      for (i=0; i < n_nodes; i++) { seed[i] = atol(argv[i+1]); }
      RANDOM_TRACE(printf("Got "); for(i=0;i<n_nodes;i++) printf("%ld ",seed[i]); printf("as new seeds.\n"));
      mpi_random_seed(n_nodes,seed);
    }
    free(seed); 
    return(TCL_OK);
  }
  else if (!strncmp(argv[0], "stat", strlen(argv[0]))) {  /* 't_random stat [status-list]' */
    RandomStatus *stat = malloc(n_nodes*sizeof(RandomStatus));
    if (argc <= 1) {
      mpi_random_stat(0,stat);
      for (i=0; i < n_nodes; i++) { 
	sprintf(buffer, "{"); Tcl_AppendResult(interp, buffer, (char *) NULL); 
	sprintf(buffer, "%ld %ld ", stat[i].idum,stat[i].iy); Tcl_AppendResult(interp, buffer, (char *) NULL);
	for (j=0; j < NTAB_RANDOM; j++) { 
	  sprintf(buffer, "%ld ", stat[i].iv[j]); Tcl_AppendResult(interp, buffer, (char *) NULL); }
	sprintf(buffer, "} "); Tcl_AppendResult(interp, buffer, (char *) NULL);
      }
    }
    else if (argc < n_nodes*(NTAB_RANDOM+2)+1) { 
      sprintf(buffer, "Wrong # of args (%d)! Usage: 't_random stat [<idum> <iy> <iv[0]> ... <iv[%d]>]^%d'", argc,NTAB_RANDOM-1,n_nodes);
      Tcl_AppendResult(interp, buffer, (char *)NULL); return (TCL_ERROR); }
    else {
      cnt = 1;
      for (i=0; i < n_nodes; i++) {
	stat[i].idum = atol(argv[cnt++]); stat[i].iy = atol(argv[cnt++]);
	for (j=0; j < NTAB_RANDOM; j++) stat[i].iv[j] = atol(argv[cnt++]);
      }
      RANDOM_TRACE(printf("Got "); for(i=0;i<n_nodes;i++) printf("%ld/%ld/... ",stat[i].idum,stat[i].iy); printf("as new status.\n"));
      mpi_random_stat(n_nodes,stat);
    }
    free(stat); 
    return(TCL_OK);
  }
  /* else */
  sprintf(buffer, "Usage: 't_random [{ int <n> | seed [<seed(0)> ... <seed(%d)>] | stat [status-list] }]'",n_nodes-1);
  Tcl_AppendResult(interp, "Unknown job '",argv[0],"' requested!\n",buffer, (char *)NULL);
  return (TCL_ERROR); 
}


/*----------------------------------------------------------------------*/

/**  Implementation of the tcl-command
     bit_random [{ seed [<seed(0)> ... <seed(n_nodes-1)>] | stat [status-list] }]
     <ul>
     <li> Without further arguments, it returns a random double between 0 and 1.
     <li> If 'seed'/'stat' is given without further arguments, it returns a tcl-list with
          the current seeds/status of the n_nodes active nodes; otherwise it issues the 
	  given parameters as the new seeds/status to the respective nodes. 
     </ul>
*/
int tclcommand_bit_random (ClientData data, Tcl_Interp *interp, int argc, char **argv) {
  char buffer[100 + TCL_DOUBLE_SPACE + 3*TCL_INTEGER_SPACE];
  int i,j,cnt; double d_out;

  if (argc == 1) {                                          /* 'bit_random' */
    d_out = bit_random_generator();
    sprintf(buffer, "%f", d_out); Tcl_AppendResult(interp, buffer, (char *) NULL); return (TCL_OK); }
  /* argc > 1 */
  argc--; argv++;
  if (!strncmp(argv[0], "seed", strlen(argv[0]))) {  /* 'bit_random seed [<seed(0)> ... <seed(n_nodes-1)>]' */
    int *seed = malloc(n_nodes*sizeof(int));
    if (argc <= 1) {
      mpi_bit_random_seed(0,seed);
      for (i=0; i < n_nodes; i++) { 
	  sprintf(buffer, "%d ", seed[i]); Tcl_AppendResult(interp, buffer, (char *) NULL); 
      }
    }
    else if (argc < n_nodes+1) { 
      sprintf(buffer, "Wrong # of args (%d)! Usage: 't_random seed [<seed(0)> ... <seed(%d)>]'", argc,n_nodes-1);
      Tcl_AppendResult(interp, buffer, (char *)NULL); return (TCL_ERROR); }
    else {
      for (i=0; i < n_nodes; i++) { seed[i] = atoi(argv[i+1]); }
      RANDOM_TRACE(printf("Got "); for(i=0;i<n_nodes;i++) printf("%d ",seed[i]); printf("as new seeds.\n"));
      mpi_bit_random_seed(n_nodes,seed);
    }
    free(seed); 
    return(TCL_OK);
  }
  else if (!strncmp(argv[0], "stat", strlen(argv[0]))) {  /* 'bit_random stat [status-list]' */
    BitRandomStatus *stat = malloc(n_nodes*sizeof(BitRandomStatus));
    if (argc <= 1) {
      mpi_bit_random_stat(0,stat);
      for (i=0; i < n_nodes; i++) { 
	sprintf(buffer, "{"); Tcl_AppendResult(interp, buffer, (char *) NULL); 
	sprintf(buffer, "%d %d ", stat[i].random_pointer_1,stat[i].random_pointer_2); Tcl_AppendResult(interp, buffer, (char *) NULL);
	for (j=0; j < MERS_BIT_RANDOM; j++) { 
	  sprintf(buffer, "%d ", stat[i].rand_w_array[j]); Tcl_AppendResult(interp, buffer, (char *) NULL); }
	sprintf(buffer, "} "); Tcl_AppendResult(interp, buffer, (char *) NULL);
      }
    }
    else if (argc < n_nodes*(MERS_BIT_RANDOM+2)+1) { 
      sprintf(buffer, "Wrong # of args (%d)! Usage: 't_random stat [<idum> <iy> <iv[0]> ... <iv[%d]>]^%d'", argc,MERS_BIT_RANDOM-1,n_nodes);
      Tcl_AppendResult(interp, buffer, (char *)NULL); return (TCL_ERROR); }
    else {
      cnt = 1;
      for (i=0; i < n_nodes; i++) {
	stat[i].random_pointer_1 = atoi(argv[cnt++]); stat[i].random_pointer_2 = atoi(argv[cnt++]);
	for (j=0; j < MERS_BIT_RANDOM; j++) stat[i].rand_w_array[j] = atoi(argv[cnt++]);
      }
      RANDOM_TRACE(printf("Got "); for(i=0;i<n_nodes;i++) printf("%d/%d/... as new status.\n",stat[i].random_pointer_1,stat[i].random_pointer_2));
      mpi_bit_random_stat(n_nodes,stat);
    }
    free(stat); 
    return(TCL_OK);
  }
  /* else */
  sprintf(buffer, "Usage: 'bit_random [{ seed [<seed(0)> ... <seed(%d)>] | stat [status-list] }]'",n_nodes-1);
  Tcl_AppendResult(interp, "Unknown job '",argv[0],"' requested!\n",buffer, (char *)NULL);
  return (TCL_ERROR); 
}



