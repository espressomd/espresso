#include "chain.h"
#include <stdio.h>
#include <stdlib.h>

/* change chains */
int mdchain(ClientData data, Tcl_Interp *interp,
	    int argc, char **argv);

/********** variables *********/
int Nchains;
Chain *chains;

void chain_init(Tcl_Interp *interp)
{
  chains = NULL;
  Tcl_CreateCommand(interp, "mdchain", mdchain, 0, NULL);
}


void chain_set_nchains(int newVal)
{
  if (chains)
    free(chains);

  chains = (Chain *)malloc(newVal*sizeof(Chain));
}

int mdchain(ClientData data, Tcl_Interp *interp,
	    int argc, char **argv)
{
  int index;
  char buffer[256];

  if (!((argc == 2) || (argc == 4))) {
    Tcl_AppendResult(interp, "wrong # args:  should be \"",
		     argv[0], " <chainindex> ?range identifier?\"",
		     (char *) NULL);
    return (TCL_ERROR);
  }

  index = atol(argv[1]);
  if (argc == 2) {
    sprintf(buffer, "%d-%d %d", chains[index].start,
	    chains[index].end, chains[index].ident);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    return (TCL_OK);
  }

  if (sscanf(argv[2], "%d-%d", &chains[index].start,
	     &chains[index].end) != 2) {
    Tcl_AppendResult(interp, "cannot parse range \"", argv[2], "\"",
		     (char *)NULL);
    return (TCL_ERROR);
  }

  chains[index].ident = atol(argv[3]);
  return (TCL_OK);
}
