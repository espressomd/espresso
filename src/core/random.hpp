/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
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
#ifndef RANDOM_H
#define RANDOM_H

/** \file random.hpp 

    A random generator
*/

#include "utils.hpp"

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
inline long l_random(void)
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
inline int i_random(int maxint)
{
  /* delivers an integer between 0 and maxint-1 */
  int temp;
  temp =  (int)( ( (double) maxint * l_random() )* AM );
  return temp;
}
  

/*----------------------------------------------------------------------*/

inline double d_random(void)
{
  /* delivers a uniform double between 0 and 1 */
  double temp;
  iy = l_random();
  if ((temp = AM * iy) > RNMX) 
    temp = RNMX;
  return temp;
}

/*----------------------------------------------------------------------*/

/** Generator for Gaussian random numbers. Uses the Box-Muller
 * transformation to generate two Gaussian random numbers from two
 * uniform random numbers.
 *
 * @return Gaussian random number.
 *
 */
inline double gaussian_random(void) {
  double x1, x2, r2, fac;
  static int calc_new = 1;
  static double save;

  /* On every second call two gaussian random numbers are calculated
     via the Box-Muller transformation. One is returned as the result
     and the second one is stored for use on the next call.
  */

  if (calc_new) {

    /* draw two uniform random numbers in the unit circle */
    do {      
      x1 = 2.0*d_random()-1.0;
      x2 = 2.0*d_random()-1.0;
      r2 = x1*x1 + x2*x2;
    } while (r2 >= 1.0 || r2 == 0.0);

    /* perform Box-Muller transformation */
    fac = sqrt(-2.0*log(r2)/r2);

    /* save one number for later use */
    save = x1*fac;
    calc_new = 0;

    /* return the second number */
    return x2*fac;

  } else {

    calc_new = 1;

    /* return the stored gaussian random number */
    return save;

  }

}

/** Generator for Gaussian random numbers. Uses the Box-Muller
 * transformation to generate two Gaussian random numbers from two
 * uniform random numbers. which generates numbers between -2 sigma and 2 sigma in the form of a Gaussian with standard deviation sigma=1.118591404 resulting in 
 * an actual standard deviation of 1.
 *
 * @return Gaussian random number.
 *
 */
inline double gaussian_random_cut(void) {
  double x1, x2, r2, fac;
  static int calc_new = 1;
  static double save, curr;

  /* On every second call two gaussian random numbers are calculated
     via the Box-Muller transformation. One is returned as the result
     and the second one is stored for use on the next call.
  */

  if (calc_new) {

    /* draw two uniform random numbers in the unit circle */
    do {      
      x1 = 2.0*d_random()-1.0;
      x2 = 2.0*d_random()-1.0;
      r2 = x1*x1 + x2*x2;
    } while (r2 >= 1.0 || r2 == 0.0);

    /* perform Box-Muller transformation */
    fac = sqrt(-2.0*log(r2)/r2);

    // save one number for later use 
    save = x1*fac*1.042267973;
    if ( fabs(save) > 2*1.042267973 ) {
      if ( save > 0 ) save = 2*1.042267973;
      else save = -2*1.042267973;
    }
    calc_new = 0;

    // return the second number 
    curr = x2*fac*1.042267973;
    if ( fabs(curr) > 2*1.042267973) {
      if ( curr > 0 ) curr = 2*1.042267973;
      else curr = -2*1.042267973;
    }
    return curr;
    
    /* save one number for later use */
    /*
    save = x1*fac*1.118591404;
    if ( fabs(save) > 2*1.118591404 ) {
      save = (2.0*d_random()-1.0)*2*1.118591404;
    }
    calc_new = 0;

    // return the second number 
    curr = x2*fac*1.118591404;
    if ( fabs(curr) > 2*1.118591404) {
      curr = (2.0*d_random()-1.0)*2*1.118591404;
    }
    return curr;
    */

  } else {

    calc_new = 1;

    /* return the stored gaussian random number */
    return save;

  }

}


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


/*----------------------------------------------------------*/

#endif






