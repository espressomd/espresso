#include "parser.h"

int parse_int_list(Tcl_Interp *interp, char *list, IntList *il)
{
  int i, tmp_argc, res = 1;
  char  **tmp_argv;
  Tcl_SplitList(interp, list, &tmp_argc, &tmp_argv);
  realloc_intlist(il, tmp_argc);
  for(i = 0 ; i < tmp_argc; i++) if (Tcl_GetInt(interp, tmp_argv[i], &(il->e[i])) == TCL_ERROR) { res = 0; break; }
  Tcl_Free((char *)tmp_argv);
  return res;
}
