#ifndef RANDOM_H
#define RANDOM_H

/** \file random.h a random generator
    <b>Responsible:</b>
    <a href="mailto:muehlbac@mpip-mainz.mpg.de">Frank</a>
    \todo Put any debug output in a DEBUG statement!!!
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

/**  A random generator for tcl.
     <li> tcl_rand() for uniform double in ]0;1[
     <li> tcl_rand(i <n>) for integer between 0 and n-1 */
int tcl_rand(ClientData data, Tcl_Interp *interp,int argc, char **argv);

/**  Implementation of the tcl-command
     setmd_random { seed [<seed(0)> ... <seed(n_nodes)>] | stat [status-list] }
     Without further arguments, it returns the current seeds/status of the nodes as a tcl-list;
     otherwise it issues the parameters as the new seeds/status to the respective nodes. */
int setmd_random(ClientData data, Tcl_Interp *interp, int argc, char **argv);

/*----------------------------------------------------------*/

#endif






