#include <stdlib.h>
#include "statistics.h"
#include "forces.h"
#include "communication.h"
#include "grid.h"

int mindist(ClientData data, Tcl_Interp *interp, int argc, char **argv)
{
  char buffer[TCL_DOUBLE_SPACE];
  double *buf;
  double mindist;
  int i;

  if (argc != 1) {
    Tcl_AppendResult(interp, "mindist takes no arguments", (char *)NULL);
    return (TCL_ERROR);
  }

  if (minimum_part_dist == -1) {
    /* here we could use a double loop (over ALL particles) to calculate... */
    Tcl_AppendResult(interp, "(not yet set)", (char *)NULL);
    return (TCL_OK);
  }

  buf = malloc(n_nodes*sizeof(double));
  mpi_gather_stats(0, buf);
  mindist = buf[0];
  for (i = 1; i < n_nodes; i++) {
    if (buf[i] < mindist)
      mindist = buf[i];
  }
  free(buf);

  if (mindist >= box_l[0] + box_l[1] + box_l[2]) {
    /* here we could use a double loop (over ALL particles) to calculate... */
    Tcl_AppendResult(interp, "(> ramp cutoffs)", (char *)NULL);
    return (TCL_OK);
  }

  Tcl_PrintDouble(interp, mindist, buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return (TCL_OK);
}
