/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
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

// ############## usage of standard C++ <random> objects
std::mt19937 generator;
std::normal_distribution<double> normal_distribution(0,1);
std::uniform_real_distribution<double> uniform_real_distribution(0,1);

void init_random(void)
{
  /* initializes the random number generator with a seed that depends on the node. You MUST NOT FORGET THIS! */
  unsigned int seed;
  seed = (10*this_node+1)*1103515245 + 12345;
  seed = (seed/65536) % 32768;
  init_random_seed((int)seed);
}

void init_random_seed(int seed)
{
  /* seed the random number generator. You MUST NOT FORGET THIS! */
  extern std::mt19937 generator;
  generator.seed(seed);
}

//checkpointing for random number generators is not yet implemented
