/*
  Copyright (C) 2010,2012,2013,2014 The ESPResSo project
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
#include <random>

void init_random(void);
void init_random_seed(int seed);

extern std::mt19937 generator;
extern std::normal_distribution<double> normal_distribution;
extern std::uniform_real_distribution<double> uniform_real_distribution;

inline double d_random() {
	/* draws a random real number from the uniform distribution in the range [0,1) */
	return uniform_real_distribution(generator); 
}

inline int i_random(int maxint){
	/* draws a random integer from the uniform distribution in the range [0,maxint-1] */
	std::uniform_int_distribution<int> uniform_int_dist(0, maxint-1);
	return uniform_int_dist(generator);
}

inline double gaussian_random(void){
	/* draws a random number from the normal distribution with mean 0 and variance 1 */
	return normal_distribution(generator);
}


/** Generator for Gaussian random numbers. Generates a Gaussian random number and generates a number between -2 sigma and 2 sigma in the form of a Gaussian with standard deviation sigma=1.118591404 resulting in 
 * an actual standard deviation of 1.
 *
 * @return Gaussian random number.
 *
 */
inline double gaussian_random_cut(void){
	double random_number=normal_distribution(generator);
	random_number=random_number*1.042267973;
	if ( fabs(random_number) > 2*1.042267973 ) {
		if ( random_number > 0 ) random_number = 2*1.042267973;
		else random_number = -2*1.042267973;
	}
	return random_number;
}

#endif
