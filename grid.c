#include "grid.h"
#include <mpi.h>
#include <stdio.h>
#include <math.h>

int setup_processor_grid()
{
  if (processor_grid[0] < 0) {
    fprintf(stderr, "not implemented: setup_processor_grid()\n");
    processor_grid[0] = nprocs;
    processor_grid[1] = processor_grid[2] = 0;
    return 1;
  }
  if (processor_grid[0]*processor_grid[1]*processor_grid[2] != nprocs)
    return 0;
  return 1;
}

int processor_grid_is_set()
{
  return (processor_grid[1] > 0);
}

int find_node(double pos[3])
{
  return map_array_node((int)floor(processor_grid[0]*pos[0]/box_l[0]),
			(int)floor(processor_grid[1]*pos[1]/box_l[1]),
			(int)floor(processor_grid[2]*pos[2]/box_l[2]));
}

void map_node_array(int node, int *a, int *b, int *c)
{
  *a = node % processor_grid[0];
  node /= processor_grid[0];
  *b = node % processor_grid[1];
  node /= processor_grid[1];
  *c = node;
}

int map_array_node(int a, int b, int c) {
  return (a + processor_grid[0]*(b + processor_grid[1]*c));
}
