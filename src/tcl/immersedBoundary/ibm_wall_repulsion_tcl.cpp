
#include "config.hpp"

#ifdef IMMERSED_BOUNDARY

#include "parser.hpp"
#include "interaction_data.hpp"
#include "immersedBoundary/ibm_wall_repulsion.hpp"

/****************
   tclcommand_inter_parse_ibm_wall_repulsion
*****************/

int tclcommand_inter_parse_ibm_wall_repulsion(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  if (!ARG_IS_S_EXACT(1, "cpt") )
  {
    // ****** Setting up a new interaction *********
    
    // Number of arguments
    if ( argc != 2  )
    {
      printf("ibm_wallRep has %d arguments.\n", argc-1);
      Tcl_AppendResult(interp, "Wall repulsion for IBM needs one parameter:\n "
                       "<kappaWall>", (char *) NULL);
      return TCL_ERROR;
    }
    
    // Check correctness of numeric parameters
    double kappaWall;
    if (! ARG_IS_D(1, kappaWall)) {
      Tcl_AppendResult(interp, "kappaWall must be double.", (char *) NULL);
      return TCL_ERROR;
    }
    
    // Setup interaction
    const int status = IBM_WallRepulsion_SetParams(bond_type, kappaWall);
    if (status == ES_ERROR)
    {
      Tcl_AppendResult(interp, "Unspecified error reading ibm_wallRep", (char *)NULL);
      return TCL_ERROR;
    }
  }
  else
  {
    // ****** Continue from checkpoint *********
    if ( argc != 2 )
    {
      Tcl_AppendResult(interp, "Wrong number of arguments reading checkpoint for ibm_wallRep ", (char *) NULL);
      return TCL_ERROR;
    }
    // Set interaction
    const int status = IBM_WallRepulsion_ResetParams(bond_type);
    if (status == ES_ERROR)
    {
      Tcl_AppendResult(interp, "Unspecified error reading ibm_wallRep checkpoint", (char *)NULL);
      return TCL_ERROR;
    }

  }
  
  return TCL_OK;
}

/***************
   tclprint_to_result_ibm_wall_repulsion
****************/

int tclprint_to_result_ibm_wall_repulsion(Tcl_Interp *interp, Bonded_ia_parameters *params)
{
    // Writes out interaction data to a text file
    // Needed for checkpointing
  
    Tcl_AppendResult(interp, "ibm_wallRep cpt", (char *) NULL);
    return TCL_OK;
}

#endif