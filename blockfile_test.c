/*
  Copyright (C) 2010 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
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
/** \file blockfile_test.c
    Reference implementation for the blockfile C-interface.
*/

#include "blockfile.h"

/*
  Demonstration of the use of the blockfile interface.
  Build using gcc -Wall -o test test.c -LLinux -lEspresso.
*/
int main()
{
  double epsilon = 1;
  double p3m_mesh_offset[3] = {.5, .5, .5};
  int    node_grid[3] = {2, 2, 2};

  FILE *f = fopen("demofile","w");
  block_writestart(f, "file");
  fprintf(f, "{Demonstration of the block format}\n");

  /* variable epsilon */
  block_writestart(f, "variable");
  fprintf(f, "epsilon ");
  block_write_data(f, TYPE_DOUBLE, 1, &epsilon);
  block_writeend(f);
  fprintf(f, "\n");

  /* variable p3m_mesh_offset */
  block_writestart(f, "variable");
  fprintf(f, "p3m_mesh_offset ");
  block_write_data(f, TYPE_DOUBLE, 1, p3m_mesh_offset);
  block_writeend(f);
  fprintf(f, "\n");

  /* variable node_grid */
  block_writestart(f, "variable");
  fprintf(f, "node_grid ");
  block_write_data(f, TYPE_DOUBLE, 1, node_grid);
  block_writeend(f);
  fprintf(f, "\n");

  /* end tag */
  block_writestart(f, "end");
  block_writeend(f);

  block_writeend(f);
  fclose(f);
  return 0;
}
