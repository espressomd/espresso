// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2009; all rights reserved unless otherwise stated.
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
