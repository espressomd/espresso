#include "triel_tcl.hpp"

#ifdef TRIELASTIC
#include "triel.hpp"
#include "communication.hpp"

/// parse parameters for the triel potential
int tclcommand_inter_parse_triel(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  int ind1, ind2, ind3;
  double max, ks, ka;
  double lo, lpo, sinpo, cospo, Area0;

  if (argc != 7 && argc !=9) {
    Tcl_AppendResult(interp, "The elastic triangle initially needs 4 parameters: "
		     "<ind1> <ind2> <ind3> <maxdist> <ks> <ka>", (char *) NULL);
    return TCL_ERROR;
  } else if (argc == 7) {

      if (! ARG_IS_I(1, ind1)) {
	Tcl_AppendResult(interp, "Triangle index must be int: "
		     "<ind1> ", (char *) NULL);
	return TCL_ERROR;
      }

      if (! ARG_IS_I(2, ind2)) {
	Tcl_AppendResult(interp, "Triangle index must be int: "
		     "<ind2> ", (char *) NULL);
	return TCL_ERROR;
      }
  
      if (! ARG_IS_I(3, ind3)) {
	Tcl_AppendResult(interp, "Triangle index must be int: "
		     "<ind3> ", (char *) NULL);
	return TCL_ERROR;
      }
  
      if (! ARG_IS_D(4, max)) {
	Tcl_AppendResult(interp, "Triangle needs max stretch : "
		     "<maxdist> ", (char *) NULL);
	return TCL_ERROR;
      }
  
      if (! ARG_IS_D(5, ks)) {
	Tcl_AppendResult(interp, "Triangle needs shear modulus (double) : "
		     "<ks> ", (char *) NULL);
	return TCL_ERROR;
      }
  
      if (! ARG_IS_D(6, ka)) {
	Tcl_AppendResult(interp, "Triangle needs dilation modulus (double) : "
		     "<ka> ", (char *) NULL);
	return TCL_ERROR;
      }
  
    CHECK_VALUE(triel_set_params(bond_type, ind1, ind2, ind3, max,ks,ka), "bond type must be nonnegative");
    
  } else {
    
      if (! ARG_IS_D(1, lo)) {
	Tcl_AppendResult(interp, "Triangle lo must be double: "
		     "<ind1> ", (char *) NULL);
	return TCL_ERROR;
      }

      if (! ARG_IS_D(2, lpo)) {
	Tcl_AppendResult(interp, "Triangle lpo must be double: "
		     "<ind2> ", (char *) NULL);
	return TCL_ERROR;
      }
  
      if (! ARG_IS_D(3, cospo)) {
	Tcl_AppendResult(interp, "Triangle cospo must be double: "
		     "<ind3> ", (char *) NULL);
	return TCL_ERROR;
      }
      
      if (! ARG_IS_D(4, sinpo)) {
	Tcl_AppendResult(interp, "Triangle sinpo must be double: "
		     "<ind3> ", (char *) NULL);
	return TCL_ERROR;
      }
      
      if (! ARG_IS_D(5, Area0)) {
	Tcl_AppendResult(interp, "Triangle Area0 must be double: "
		     "<ind3> ", (char *) NULL);
	return TCL_ERROR;
      }
  
      if (! ARG_IS_D(6, max)) {
	Tcl_AppendResult(interp, "Triangle needs max stretch : "
		     "<maxdist> ", (char *) NULL);
	return TCL_ERROR;
      }
  
      if (! ARG_IS_D(7, ks)) {
	Tcl_AppendResult(interp, "Triangle needs shear modulus (double) : "
		     "<ks> ", (char *) NULL);
	return TCL_ERROR;
      }
  
      if (! ARG_IS_D(8, ka)) {
	Tcl_AppendResult(interp, "Triangle needs dilation modulus (double) : "
		     "<ka> ", (char *) NULL);
	return TCL_ERROR;
      }
    
      CHECK_VALUE(triel_reset_params(bond_type, lo, lpo, cospo, sinpo, Area0, max,ks,ka), "bond type must be nonnegative");
    
  }
}

int tclprint_to_result_trielIA(Tcl_Interp *interp, Bonded_ia_parameters *params)
{
  char buffer[TCL_DOUBLE_SPACE];

   Tcl_PrintDouble(interp, params->p.triel.lo, buffer);
  Tcl_AppendResult(interp, "triel ", buffer," ", (char *) NULL);
   Tcl_PrintDouble(interp, params->p.triel.lpo, buffer);
  Tcl_AppendResult(interp," ",buffer, " ", (char *) NULL);
   Tcl_PrintDouble(interp, params->p.triel.cospo, buffer);
  Tcl_AppendResult(interp," ", buffer, " ", (char *) NULL);
   Tcl_PrintDouble(interp, params->p.triel.sinpo, buffer);
    Tcl_AppendResult(interp," ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, params->p.triel.Area0, buffer);
  Tcl_AppendResult(interp," ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, params->p.triel.maxdist, buffer);
  Tcl_AppendResult(interp," ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, params->p.triel.ks, buffer);
  Tcl_AppendResult(interp," ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, params->p.triel.ka, buffer);
  Tcl_AppendResult(interp," ", buffer, " ", (char *) NULL);
  return TCL_OK;
}

#endif
