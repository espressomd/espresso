#ifndef _REACTION_ENSEMBLE_TCL_H
#define _REACTION_ENSEMBLE_TCL_H
#include "parser.hpp"

#ifdef REACTION_ENSEMBLE

/*********************************
 * functions
 *********************************/
/** Implementation of the Tcl command \ref tclcommand_reaction_ensemble. This function
 *  allows to change the parameters of reaction_ensemble */
int tclcommand_reaction_ensemble(ClientData data, Tcl_Interp *interp, int argc, char **argv);

#endif

#endif
