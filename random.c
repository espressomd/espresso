#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"
#include "global.h"
#include "debug.h"
#include "random.h"
#include "tcl.h"
#include "communication.h"
#include "utils.h"

/** \file random.c A random generator. Be sure to run init_random() before
    you use any of the generators. */

const long     IA = 16807;
const long     IM = 2147483647;
const double   AM = (1.0/2147483647. );
const long     IQ = 127773;
const long     IR = 2836;
const double NDIV = (double) (1+(2147483647-1)/NTAB_RANDOM);
const double RNMX = (1.0-1.2e-7);


static long  idum = 1;
static long  idumInit = 1;
static long  iy=0;
static long  iv[NTAB_RANDOM];

/*----------------------------------------------------------------------*/

long l_random(void)
{
  /* 
   *    from Numerical Recipes in C by Press et al.,
   *    N O T E   T H A T   T H E R E   A R E   N O   S A F E T Y   C H E C K S  !!!
   */
  
  int    j;
  long   k;
  
  k = (idum) / IQ;
  idum = IA * (idum - k * IQ) - IR * k;
  if (idum < 0) idum += IM;
  j = iy / NDIV;
  iy = iv[j];
  iv[j] = idum;
  return iy;
}

/*----------------------------------------------------------------------*/

int i_random(int maxint)
{
  /* delivers an integer between 0 and maxint-1 */
  int temp;
  temp =  (int)( ( (double) maxint * l_random() )* AM );
  return temp;
}
  

/*----------------------------------------------------------------------*/

double d_random(void)
{
  /* delivers a uniform double between 0 and 1 */
  double temp;
  iy = l_random();
  if ((temp = AM * iy) > RNMX) 
    temp = RNMX;
  return temp;
}

/*----------------------------------------------------------------------*/

void init_random(void)
{
  /* initializes the random number generator. You MUST NOT FORGET THIS! */
  
  unsigned long seed;
  seed = (10*this_node+1)*1103515245 + 12345;
  seed = (seed/65536) % 32768;
  init_random_seed((long)seed);
}

/*----------------------------------------------------------------------*/

void init_random_seed(long seed)
{
  /* initializes the random number generator. You MUST NOT FORGET THIS! */

  int    j;
  long   k;

  /* This random generator is bad I know, thats why its only used
     for the seed (see Num. Rec. 7.1.) */
  idumInit = idum = seed;
  RANDOM_TRACE(fprintf(stderr, "%d: Init random with seed %ld in 'random.c'\n",this_node,idum));
  for (j = NTAB_RANDOM + 7;j >= 0; j--) {
    k = (idum) / IQ;
    idum = IA * (idum - k * IQ) - IR * k;
    if (idum < 0) idum += IM;
    if (j < NTAB_RANDOM) iv[j] = idum;
  }
  iy = iv[0];
}

/*----------------------------------------------------------------------*/

void init_random_stat(RandomStatus my_stat) {
  /* initializes the random number generator to a given status */
  int i;

  idum = my_stat.idum; iy = my_stat.iy;
  for (i=0; i < NTAB_RANDOM; i++) iv[i] = my_stat.iv[i];
}

/*----------------------------------------------------------------------*/

long print_random_idum(void) {
  /* returns current 'idum' */
  return(idum);
}

/*----------------------------------------------------------------------*/

long print_random_seed(void) {
  /* returns the seed originally used upon last initialize of the generator */
  return(idumInit);
}

/*----------------------------------------------------------------------*/

RandomStatus print_random_stat(void) {
  /* returns current status of random number generator */
  RandomStatus my_stat; int i;
  
  my_stat.idum = idum; my_stat.iy = iy;
  for (i=0; i < NTAB_RANDOM; i++) my_stat.iv[i] = iv[i];
  return(my_stat);
}

/*----------------------------------------------------------------------*/


/**  Implementation of the tcl-command
     t_random [{ int <n> | seed [<seed(0)> ... <seed(n_nodes-1)>] | stat [status-list] }]
     <li> Without further arguments, it returns a random double between 0 and 1.
     <li> If 'int <n>' is given, it returns a random integer between 0 and n-1.
     <li> If 'seed'/'stat' is given without further arguments, it returns a tcl-list with
          the current seeds/status of the n_nodes active nodes; otherwise it issues the 
	  given parameters as the new seeds/status to the respective nodes. */
int t_random (ClientData data, Tcl_Interp *interp, int argc, char **argv) {
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
  int i,j,cnt, i_out; double d_out;

  if (argc == 1) {                                          /* 't_random' */
    d_out = d_random();
    sprintf(buffer, "%f", d_out); Tcl_AppendResult(interp, buffer, (char *) NULL); return (TCL_OK); }
  else {
    argc--; argv++;
    if (!strncmp(argv[0], "int", strlen(argv[0]))) {        /* 't_random int <n>' */
      if(argc < 2){ Tcl_AppendResult(interp, "Wrong # of args: Usage: 't_random int <n>'", (char *) NULL); return (TCL_ERROR); }
      else {
	Tcl_GetInt(interp, argv[2], &i_out);
	i_out = i_random(i_out);
	sprintf(buffer, "%d", i_out); Tcl_AppendResult(interp, buffer, (char *) NULL); return (TCL_OK);
      }
    }
    else if (!strncmp(argv[0], "seed", strlen(argv[0]))) {  /* 't_random seed [<seed(0)> ... <seed(n_nodes-1)>]' */
      long seed[n_nodes];
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
      // free(seed); 
      return(TCL_OK);
    }
    else if (!strncmp(argv[0], "stat", strlen(argv[0]))) {  /* 't_random stat [status-list]' */
      RandomStatus stat[n_nodes];
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
      // free(stat); 
      return(TCL_OK);
    }
    else { 
      sprintf(buffer, "Usage: 't_random [{ int <n> | seed [<seed(0)> ... <seed(%d)>] | stat [status-list] }]'",n_nodes-1);
      Tcl_AppendResult(interp, "Unknown job '",argv[0],"' requested!\n",buffer, (char *)NULL); return (TCL_ERROR); 
    }
  }
  printf("Unknown error in tcl_seed(): This should have never been shown!\n");
  return(TCL_ERROR);
}
