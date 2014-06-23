

#include "external_potential_tcl.hpp"

int tclcommand_external_potential_tabulated(Tcl_Interp* interp, int argc, char **argv, ExternalPotential* e);

int tclcommand_external_potential(ClientData _data, Tcl_Interp *interp,
    int argc, char **argv) {
  ExternalPotential* e;
  int error = generate_external_potential(&e);
  if (error == ES_ERROR)
    return TCL_ERROR;

  if (argc<2) {
    Tcl_AppendResult(interp, "Usage: external_potential <tabulate|rod|...>\n" , (char *)NULL);
    return TCL_ERROR;
  }
  if (ARG1_IS_S("tabulated")) {
    return tclcommand_external_potential_tabulated(interp, argc-2, argv+2, e);
  }
  Tcl_AppendResult(interp, "Usage: external_potential <tabulate|rod|...>\n" , (char *)NULL);
  return TCL_ERROR;

}



int tclcommand_external_potential_tabulated(Tcl_Interp* interp, int argc, char **argv, ExternalPotential* e) 
{
  char* filename =0;

  DoubleList scalelist;

  init_doublelist(&scalelist);

  while (argc>0) {
    if (ARG0_IS_S("file") ) {
      if (argc>1) {
        filename = argv[1];
        argc-=2;
        argv+=2;
      } else {
        Tcl_AppendResult(interp, "Usage: external_potential file <filename>\n" , (char *)NULL);
        return TCL_ERROR;
      }
    } else if (ARG0_IS_S("scale")) {
      if (argc>1  && ARG_IS_DOUBLELIST(1, scalelist)) {
        argc-=2;
        argv+=2;
      } else {
        Tcl_AppendResult(interp, "Usage: external_potential tabulated scale <float>\n" , (char *)NULL);
        return TCL_ERROR;
      }
    } else {
      Tcl_AppendResult(interp, "Unknown argument to external_potential: " , argv[0], "\n", (char *)NULL);
      return TCL_ERROR;
    }
  }
  if (filename == 0) {
    Tcl_AppendResult(interp, "No filename given to external_potential tabulated\n" , (char *)NULL);
    return TCL_ERROR;
  }
  if (external_potential_tabulated_init(n_external_potentials-1, filename, scalelist.n, scalelist.e)==ES_ERROR) {
    Tcl_AppendResult(interp, "Error creating external potential\n" , (char *)NULL);
    return TCL_ERROR;
  }
  return gather_runtime_errors(interp, TCL_OK);
}




