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

/**  Implementation of the tcl command \ref tcl_t_random. Access to the
     parallel random number generator.
*/
int t_random(ClientData data, Tcl_Interp *interp, int argc, char **argv);

/*----------------------------------------------------------*/

#endif






