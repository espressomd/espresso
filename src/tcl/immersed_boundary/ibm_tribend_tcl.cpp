
#include "config.hpp"

#ifdef IMMERSED_BOUNDARY

#include "parser.hpp"
#include "interaction_data.hpp"
#include "immersed_boundary/ibm_tribend.hpp"

/****************
   tclcommand_inter_parse_ibm_tribend
*****************/

int tclcommand_inter_parse_ibm_tribend(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  if (!ARG_IS_S_EXACT(1, "cpt") )
  {
    // ****** Setting up a new interaction *********
    
    // Number of arguments
    if ( argc != 8  )
    {
      printf("ibm_tribend has %d arguments.\n", argc-1);
      Tcl_AppendResult(interp, "Tribend for IBM needs seven parameters:\n "
                       "<ind1> <ind2> <ind3> <ind4> <method> <kb> <flat/initial>", (char *) NULL);
      return TCL_ERROR;
    }
    
    // PARTNER LIST ONLY REQUIRED FOR KRUEGER BENDING
    
    int ind1, ind2, ind3, ind4;
    double kb;
    
    if (! ARG_IS_I(1, ind1)) {
      Tcl_AppendResult(interp, "Node index must be int.", (char *) NULL);
      return TCL_ERROR;
    }
    
    if (! ARG_IS_I(2, ind2)) {
      Tcl_AppendResult(interp, "Node index must be int.", (char *) NULL);
      return TCL_ERROR;
    }
    
    if (! ARG_IS_I(3, ind3)) {
      Tcl_AppendResult(interp, "Node index must be int.", (char *) NULL);
      return TCL_ERROR;
    }
    
    if (! ARG_IS_I(4, ind4)) {
      Tcl_AppendResult(interp, "Node index must be int.", (char *) NULL);
      return TCL_ERROR;
    }
    
    // Bending method
    bool done = false;
    tBendingMethod method;
    if (strcmp(argv[5], "TriangleNormals") == 0) { method = TriangleNormals; done = true; }
    if (strcmp(argv[5], "NodeNeighbors") == 0) { method = NodeNeighbors; done = true; }
    if ( !done )
    {
      Tcl_AppendResult(interp, "Wrong bending method: TriangleNormals or NodeNeighbors ", (char *) NULL);
      return TCL_ERROR;
    }
    
    if (! ARG_IS_D(6, kb))
    {
      Tcl_AppendResult(interp, "Triangle needs bending modulus (double)", (char *) NULL);
      return TCL_ERROR;
    }
    
    bool flat = false;
    if ( strcmp(argv[7], "flat") == 0 ) flat = true;
    else if ( strcmp(argv[7], "initial") == 0 ) flat = false;
    else { Tcl_AppendResult(interp, "Equilibrium angle must be flat or initial", (char *) NULL); return TCL_ERROR; }
    
    // Set interaction
    CHECK_VALUE(IBM_Tribend_SetParams(bond_type, ind1, ind2, ind3, ind4, method, kb, flat), "failed to setup bonding interaction");
    
    // Everyting ok
    return TCL_OK;
  }
  else
  {
    // ****** Continue from checkpoint *********
    if ( argc != 3 )
    {
      Tcl_AppendResult(interp, "Wrong number of arguments reading checkpoint for ibm_tribend ", (char *) NULL);
      return TCL_ERROR;
    }
    
    double k;
    if (! ARG_IS_D(2, k))
    {
      Tcl_AppendResult(interp, "Checkpoint reading failed for tribend ", (char *) NULL);
      return TCL_ERROR;
    }
    CHECK_VALUE(IBM_Tribend_ResetParams(bond_type, k), "failed to reset bonding interaction");

    return TCL_OK;

  }
}

/***************
   tclprint_to_result_ibm_tribend
****************/

int tclprint_to_result_ibm_tribend(Tcl_Interp *interp, Bonded_ia_parameters *params)
{
  char buffer[TCL_DOUBLE_SPACE];
  
  // Reference angle
  Tcl_PrintDouble(interp, params->p.ibm_tribend.kb, buffer);
  Tcl_AppendResult(interp, "ibm_tribend cpt ", buffer," ", (char *) NULL);
  
  return TCL_OK;
}

#endif