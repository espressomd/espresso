#ifndef RANDOM_H
#define RANDOM_H

/** \file random.h a random generator
    <b>Responsible:</b>
    <a href="mailto:muehlbac@mpip-mainz.mpg.de">Frank</a>
*/

/*----------------------------------------------------------*/

# define NTAB_RANDOM  32

typedef struct {
  long  idum;
  long  iy;
  long  iv[NTAB_RANDOM];
} RandomStatus;

/*----------------------------------------------------------*/

extern long   l_random(void);
extern int    i_random(int maxint);
extern double d_random(void);
extern void   init_random(void);
extern void   init_random_seed(long seed);
extern void   init_random_stat(RandomStatus my_stat);
extern long   print_random_idum(void);
extern long   print_random_seed(void);
RandomStatus  print_random_stat(void);

/**  Implementation of the tcl-command
     t_random [{ int <n> | seed [<seed(0)> ... <seed(n_nodes-1)>] | stat [status-list] }]
     <li> Without further arguments, it returns a random double between 0 and 1.
     <li> If 'int <n>' is given, it returns a random integer between 0 and n-1.
     <li> If 'seed'/'stat' is given without further arguments, it returns a tcl-list with
          the current seeds/status of the n_nodes active nodes; otherwise it issues the 
	  given parameters as the new seeds/status to the respective nodes. */
int t_random(ClientData data, Tcl_Interp *interp, int argc, char **argv);

/*----------------------------------------------------------*/

#endif






