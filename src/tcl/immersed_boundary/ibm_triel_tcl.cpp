

#include "config.hpp"

#ifdef IMMERSED_BOUNDARY

#include "parser.hpp"
#include "interaction_data.hpp"
#include "immersed_boundary/ibm_triel.hpp"

/****************
   tclcommand_inter_parse_ibm_triel
*****************/

int tclcommand_inter_parse_ibm_triel(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  if (!ARG_IS_S_EXACT(1, "cpt") )
  {
    // ****** Setting up a new interaction *********
    
    // Number of arguments
    if ( argc != 7 && argc != 8 )
    {
      printf("ibm_triel has %d arguments.\n", argc-1);
      Tcl_AppendResult(interp, "Triel for IBM needs six or seven parameters:\n "
                       "<ind1> <ind2> <ind3> <maxdist> <elastic_law> <elastic constant(s)>", (char *) NULL);
      return TCL_ERROR;
    }
    
    
    // ******* This is the syntax for setting up a new triel form the tcl **********
    
    int ind1, ind2, ind3;
    double max;
    
    // Check correctness of numeric parameters
    if (! ARG_IS_I(1, ind1)) {
      Tcl_AppendResult(interp, "Triangle index 1 must be int", (char *) NULL);
      return TCL_ERROR;
    }
    
    if (! ARG_IS_I(2, ind2)) {
      Tcl_AppendResult(interp, "Triangle index 2 must be int.", (char *) NULL);
      return TCL_ERROR;
    }
    
    if (! ARG_IS_I(3, ind3)) {
      Tcl_AppendResult(interp, "Triangle index 3 must be int.", (char *) NULL);
      return TCL_ERROR;
    }
    
    if (! ARG_IS_D(4, max)) {
      Tcl_AppendResult(interp, "Triangle needs max stretch", (char *) NULL);
      return TCL_ERROR;
    }
   
    // Check elastic law
    tElasticLaw law = NeoHookean;
    bool done = false;
    if ( strcmp(argv[5], "NeoHookean") == 0) { law = NeoHookean; done = true; }
    if ( strcmp(argv[5], "Skalak") == 0) { law = Skalak; done = true; }
    if ( !done )
    {
      Tcl_AppendResult(interp, "Law must be NeoHookean or Skalak", (char *) NULL);
      return TCL_ERROR;
    }
    
    // Elastic constants
    double k1=0, k2 = 0;
    if ( law == NeoHookean )
    {
      if ( argc != 7 )
      {
        Tcl_AppendResult(interp, "NeoHookean needs one elastic constant", (char *) NULL);
        return TCL_ERROR;
      }
      else
      {
        if (! ARG_IS_D(6, k1))
        {
          Tcl_AppendResult(interp, "Spring constant must be a double", (char *) NULL);
          return TCL_ERROR;
        }
      }
    }
    
    if ( law == Skalak )
    {
      if ( argc != 8 )
      {
        Tcl_AppendResult(interp, "Skalak needs two elastic constants", (char *) NULL);
        return TCL_ERROR;
      }
      else
      {
        if (! ARG_IS_D(6, k1))
        {
          Tcl_AppendResult(interp, "Spring constant must be a double", (char *) NULL);
          return TCL_ERROR;
        }
        if (! ARG_IS_D(7, k2))
        {
          Tcl_AppendResult(interp, "Spring constant must be a double", (char *) NULL);
          return TCL_ERROR;
        }
      }
    }
    
    // Set values
    const int status = IBM_Triel_SetParams(bond_type, ind1, ind2, ind3, max, law, k1,k2);
    if (status == ES_ERROR)
    {
      Tcl_AppendResult(interp, "Unspecified error reading ibm_triel", (char *)NULL);
      return TCL_ERROR;
    }
    
    // Everyting ok
    return TCL_OK;

    

    
    
    // Argument 5 is passed as a single string to allow flexibility
    // Copy and tokenize argument 5
    // Needs to be copied, otherwise the original tcl strign is modified
/*    char elasticLawString[400];
    strcpy(elasticLawString, argv[5]);
    char *tok = strtok(elasticLawString, " ");
    
    // Check elastic law
    tElasticLaw law = NeoHookean;
    bool done = false;
    if ( tok == NULL )
    {
      Tcl_AppendResult(interp, "No elastic law given", (char *) NULL);
      return TCL_ERROR;
    }
    if ( strcmp(tok, "NeoHookean") == 0) { law = NeoHookean; done = true; }
    if ( strcmp(tok, "Skalak") == 0) { law = Skalak; done = true; }
    if ( !done )
    {
      Tcl_AppendResult(interp, "Law must be NeoHookean or Skalak", (char *) NULL);
      return TCL_ERROR;
    }
    
    tok = strtok(NULL, " ");
    if ( tok == NULL )
    {
      Tcl_AppendResult(interp, "Elastic law needs parameters: k for NeoHookean and k1, k2 for Skalak", (char *) NULL);
      return TCL_ERROR;
    }
    const double k1 = atof(tok);
    
    tok = strtok(NULL, " ");
    if ( law == NeoHookean && tok != NULL )
    {
      printf("tok = %s\n", tok);
      Tcl_AppendResult(interp, "NeoHookean law has only one elastic constant", (char *) NULL);
      return TCL_ERROR;
    }
    if ( law == Skalak && tok == NULL )
    {
      Tcl_AppendResult(interp, "Skalak law needs two elastic constants", (char *) NULL);
      return TCL_ERROR;
    }
    // Read 2nd k only if Skalak
    double k2 = 0;
    if ( law == Skalak )
      k2 = atof(tok);
    
    // Now everything is there. Check if there are not further arguments
    tok = strtok(NULL, " ");
    
    if ( tok == NULL )
    {
      // No further parameters --> compute triangle data automatically at start of simulation
      // Set interaction
      const int status = IBM_Triel_SetParams(bond_type, ind1, ind2, ind3, max, law, k1,k2);
      if (status == ES_ERROR)
      {
        Tcl_AppendResult(interp, "Unspecified error reading ibm_triel", (char *)NULL);
        return TCL_ERROR;
      }
      
      // Everyting ok
      return TCL_OK;
    }
    else
    {
      Tcl_AppendResult(interp, "Too many parameters reading triel", (char *)NULL);
      return TCL_ERROR;
    }*/
  }
  else
  {
    // ****** Continue from checkpoint *********
    if ( argc != 4 )
    {
      Tcl_AppendResult(interp, "Wrong number of arguments reading checkpoint for ibm_triel ", (char *) NULL);
      return TCL_ERROR;
    }
    
    double k1, l0;

    // Reference shape
    if (! ARG_IS_D(2, k1))
    {
      Tcl_AppendResult(interp, "Error reading triel checkpoint: k1 must be double", (char *) NULL);
      return TCL_ERROR;
    }
    
    if (! ARG_IS_D(3, l0))
    {
      Tcl_AppendResult(interp, "Error reading triel checkpoint: l0 must be a double", (char *) NULL);
      return TCL_ERROR;
    }
    
   
    // Set interaction
    const int status = IBM_Triel_ResetParams(bond_type, k1, l0);
    if (status == ES_ERROR)
    {
      Tcl_AppendResult(interp, "Unspecified error reading ibm_triel checkpoint", (char *)NULL);
      return TCL_ERROR;
    }

  }
  
  return TCL_OK;
}

/***************
   tclprint_to_result_ibm_triel
****************/

int tclprint_to_result_ibm_triel(Tcl_Interp *interp, Bonded_ia_parameters *params)
{
  char buffer[TCL_DOUBLE_SPACE];
  
  // Reference shape
  Tcl_PrintDouble(interp, params->p.ibm_triel.k1, buffer);
  Tcl_AppendResult(interp, "ibm_triel cpt ", buffer," ", (char *) NULL);
  
  Tcl_PrintDouble(interp, params->p.ibm_triel.l0, buffer);
  Tcl_AppendResult(interp," ",buffer, " ", (char *) NULL);
  
  return TCL_OK;
}

#endif