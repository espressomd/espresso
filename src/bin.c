/*
  Copyright (C) 2010,2011 The ESPResSo project
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
#include "bin.h"
#include "utils.h"

void setup_linear_bins(DoubleList *dl, double min_bin, double max_bin, int bins)
{
  int i;
  realloc_doublelist(dl, dl->n = bins + 1);
  for (i = 0; i <= bins; i++)
    dl->e[i] = min_bin + ((max_bin - min_bin)/bins)*i;
}

void setup_log_bins(DoubleList *dl, double min_bin, double max_bin, int bins)
{
  int i;
  realloc_doublelist(dl, dl->n = bins + 1);
  for (i = 0; i <= bins; i++)
    dl->e[i] = min_bin*pow(max_bin/min_bin, ((double)i)/bins);
}


