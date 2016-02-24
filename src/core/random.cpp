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
#include "global.hpp"
#include "random.hpp"
#include "communication.hpp"

/** \file random.cpp A random generator. 
    Be sure to run init_random() before you use any of the generators. */

/* Stuff for Franks ran1-generator */
long  idum = -1;
long  idumInit = -1;
long  iy=0;
long  iv[NTAB_RANDOM];

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

  /* This random generator is bad I know {why, Frank? It's the same as the
     one in l_random!}, thats why its only {no, in l_random as well!} used
     for the seed (see Num. Rec. 7.1.) */
  if (seed < 1) {
    fprintf(stderr,"The initial seed of the random number generator must be a positive integer!\n");
    fprintf(stderr,"Using 0 will result in a plain 0-sequence, hence it's forbidden (you used: %ld)!\n",seed);
    fflush(NULL); 
    errexit();
  }
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
