#include "random.h"
#include "tcl.h"
#include "communication.h"

/** \file random.c A random generator. Be sure to run init_random() before
    you use any of the generators. */

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)


static long  idum = 1;
static long  iy=0;
static long  iv[NTAB];

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
  temp =  (int)( ( maxint * l_random() )* AM );
  return temp;
}
  

/*----------------------------------------------------------------------*/

double random(void)
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

  int    j;
  long   k;
  unsigned long seed;

  /* This random generator is bad I know, thats why its only used
     for the seed (see Num. Rec. 7.1.) */
  seed = (10*this_node+1)*1103515245 + 12345;
  seed = (seed/65536) % 32768;
  idum = (long) seed;
  printf("%d init random with seed %ld\n",this_node,idum);
  for (j = NTAB + 7;j >= 0; j--) {
    k = (idum) / IQ;
    idum = IA * (idum - k * IQ) - IR * k;
    if (idum < 0) idum += IM;
    if (j < NTAB) iv[j] = idum;
  }
  iy = iv[0];
}

/*----------------------------------------------------------------------*/
