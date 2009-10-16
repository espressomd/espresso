// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2009; all rights reserved unless otherwise stated.
#ifndef RANDOM_H
#define RANDOM_H

/** \file random.h 

    A random generator
*/

/*----------------------------------------------------------*/

/* Stuff for Franks ran1-generator */
/*@{*/
#define IA 16807
#define IM 2147483647
#define AM (1.0/2147483647.)
#define IQ 127773
#define IR 2836
#define NDIV ((double) (1+(2147483647-1)/NTAB_RANDOM))
#define RNMX (1.0-1.2e-7)
# define NTAB_RANDOM  32

extern long  idum;
extern long  idumInit;
extern long  iy;
extern long  iv[NTAB_RANDOM];
/*@}*/

/** Stuff for Burkhards r250-generator */
/*@{*/
#define MERS1 147
#define MERS_BIT_RANDOM 250
#define NBIT 32
#define BIGINTEGER 2147483647
#define BIGFLOAT 2147483647.
#define FACTOR 4.6566128752457969e-10
#define MULTIPLY 16807.
#define NWARM 10000

extern int bit_seed;
extern int rand_w_array[MERS_BIT_RANDOM];
extern int random_pointer_1;
extern int random_pointer_2;
/*@}*/

typedef struct {
  long  idum;
  long  iy;
  long  iv[NTAB_RANDOM];
} RandomStatus;

void   init_random(void);
void   init_random_seed(long seed);
void   init_random_stat(RandomStatus my_stat);
long   print_random_idum(void);
long   print_random_seed(void);
RandomStatus  print_random_stat(void);

/** classical RAN1 random number generator */
MDINLINE long l_random(void)
{
  /* 
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

/** same as l_random, but for integer */
MDINLINE int i_random(int maxint)
{
  /* delivers an integer between 0 and maxint-1 */
  int temp;
  temp =  (int)( ( (double) maxint * l_random() )* AM );
  return temp;
}
  

/*----------------------------------------------------------------------*/

MDINLINE double d_random(void)
{
  /* delivers a uniform double between 0 and 1 */
  double temp;
  iy = l_random();
  if ((temp = AM * iy) > RNMX) 
    temp = RNMX;
  return temp;
}


/**  Implementation of the tcl command \ref tcl_t_random. Access to the
     parallel random number generator.
*/
int t_random(ClientData data, Tcl_Interp *interp, int argc, char **argv);

/*----------------------------------------------------------*/
/*----------------------------------------------------------*/
/*----------------------------------------------------------*/

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






