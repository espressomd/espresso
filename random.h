// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
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
/*----------------------------------------------------------*/
/*----------------------------------------------------------*/

#define MERS_BIT_RANDOM 250

typedef struct {
  int random_pointer_1;
  int random_pointer_2;
  int rand_w_array[MERS_BIT_RANDOM];
} BitRandomStatus;

double bit_random_generator(void);
void   init_bit_random(void);
void   init_bit_random_generator(int iseed);
void   init_bit_random_stat(BitRandomStatus my_stat);
int    print_bit_random_seed(void);
BitRandomStatus print_bit_random_stat(void);

/**  Implementation of the tcl command \ref tcl_bit_random. 
     Access to the parallel bit random number generator.
*/
int bit_random(ClientData data, Tcl_Interp *interp, int argc, char **argv);

/*----------------------------------------------------------*/

#endif






