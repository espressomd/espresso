#include "blockfile.h"

/*
  Demonstration of the use of the blockfile interface.
  Build using gcc -Wall -o test test.c -LLinux -ltcl_md.
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
