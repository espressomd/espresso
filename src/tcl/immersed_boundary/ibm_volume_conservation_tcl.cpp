

#include "config.hpp"

#ifdef IMMERSED_BOUNDARY

#include "parser.hpp"
#include "interaction_data.hpp"
#include "immersed_boundary/ibm_volume_conservation_tcl.hpp"
#include "immersed_boundary/ibm_volume_conservation.hpp"

/****************
   tclcommand_inter_parse_ibm_volcons
*****************/

int tclcommand_inter_parse_ibm_volume_conservation(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  if (!ARG_IS_S_EXACT(1, "cpt") )
  {
    // ****** Setting up a new interaction *********
    
    // Number of arguments
    if ( argc != 3  )
    {
      printf("ibm_volcons has %d arguments.\n", argc-1);
      Tcl_AppendResult(interp, "VolCons for IBM needs two parameters:\n "
                       "<softID> <kappaV>", (char *) NULL);
      return TCL_ERROR;
    }
    
    int softID;
    double kappaV;
    
    // Check correctness of numeric parameters
    if (! ARG_IS_I(1, softID)) {
      Tcl_AppendResult(interp, "softID must be int\n", (char *) NULL);
      return TCL_ERROR;
    }
    
    if (! ARG_IS_D(2, kappaV)) {
      Tcl_AppendResult(interp, "kappaV must be double\n", (char *) NULL);
      return TCL_ERROR;
    }
    
/*    bool writeCom = false;
    if ( argc == 4 )
    {
      if ( ARG_IS_S_EXACT(3, "writeCOM") ) writeCom = true;
      else
      {
        Tcl_AppendResult(interp, "illegal optional argument: must be writeCOM or nothing", (char *) NULL);
        return TCL_ERROR;
      }
    }*/
    
    // Setup interaction
    const int status = IBM_VolumeConservation_SetParams(bond_type, softID, kappaV);
    if (status == ES_ERROR)
    {
      Tcl_AppendResult(interp, "Unspecified error reading volcons", (char *)NULL);
      return TCL_ERROR;
    }
  }
  else
  {
    // ****** Continue from checkpoint *********
    if ( argc != 3 )
    {
      Tcl_AppendResult(interp, "Wrong number of arguments reading checkpoint for ibm_volcons ", (char *) NULL);
      return TCL_ERROR;
    }
    
    // ****** Continue from checkpoint *********
    double volRef;
    if (! ARG_IS_D(2, volRef))
    {
      Tcl_AppendResult(interp, "volRef must be double ", (char *) NULL);
      return TCL_ERROR;
    }
    
    // Reset interaction
    const int status = IBM_VolumeConservation_ResetParams(bond_type, volRef);
    if (status == ES_ERROR)
    {
      Tcl_AppendResult(interp, "Unspecified error reading volcons checkpoint", (char *)NULL);
      return TCL_ERROR;
    }
  }
  
  return TCL_OK;
}

/***************
   tclprint_to_result_ibm_volcons
****************/
    
int tclprint_to_result_ibm_volume_conservation(Tcl_Interp *interp, Bonded_ia_parameters *params)
{
  // Writes out interaction data to a text file
  // Needed for checkpointing
  
  char buffer[TCL_DOUBLE_SPACE];
  
  Tcl_PrintDouble(interp, params->p.ibmVolConsParameters.volRef, buffer);
  Tcl_AppendResult(interp, "ibm_volcons cpt ", buffer," ", (char *) NULL);
  
  return TCL_OK;
}

#endif
