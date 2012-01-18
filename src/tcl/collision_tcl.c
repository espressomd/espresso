#include "particle_data.h"
#include "interaction_data.h"
#include "virtual_sites_relative.h"
#include "collision.h"
#include "virtual_sites.h"
#include "integrate.h"
#include "cells.h"
#include "communication.h" 
#include "parser.h" 


#ifdef COLLISION_DETECTION

int tclcommand_on_collision(ClientData data, Tcl_Interp *interp, int argc, char **argv) 
{
 // If no argumens are given, print status
 if (argc==1)
 {
  char s[200];
  if (collision_detection_mode==0)
  {
   sprintf(s,"off");
  }
  
  if (collision_detection_mode==1)
  {
   sprintf(s,"bind_centers %f %d",collision_distance,collision_detection_bond_centers);
  }
  
  if (collision_detection_mode==2)
  {
   sprintf(s,"bind_at_point_of_collision %f %d %d %d",collision_distance,collision_detection_bond_centers, collision_detection_bond_vs,collision_vs_particle_type);
  }
  Tcl_AppendResult(interp, s, (char*) NULL);
  return TCL_OK;
 }

 // Otherwise, we set parametes

 int res;

 if (ARG1_IS_S("off"))
 {
  collision_detection_set_params(0,0.,0,0,0);
 }
 else if (ARG1_IS_S("bind_centers"))
 {
  if (argc!=4)
  {
   Tcl_AppendResult(interp, "Need a ditance and a bond type as args.", (char*) NULL);
   return TCL_ERROR;
  }
  double d;
  if (!ARG_IS_D(2,d))
  {
   Tcl_AppendResult(interp, "Need a ditance as 1st arg.", (char*) NULL);
   return TCL_ERROR;
  }
  int bond1;
  if (!ARG_IS_I(3,bond1))
  {
   Tcl_AppendResult(interp, "Need a bond type as 2nd argument.", (char*) NULL);
   return TCL_ERROR;
  }
  res=collision_detection_set_params(1,d,bond1,0,0);
  if (res==2)
  {
   Tcl_AppendResult(interp, "Collision detection only works on a single cpu.", (char*) NULL);
   return TCL_ERROR;
  }
  if (res==3)
  {
   Tcl_AppendResult(interp, "A bond of the specified type does not exist.", (char*) NULL);
   return TCL_ERROR;
  }
 }
 else if (ARG1_IS_S("bind_at_point_of_collision"))
 {
  if (argc!=6)
  {
   Tcl_AppendResult(interp, "Need a ditance, two bond types, and a particle type as args.", (char*) NULL);
   return TCL_ERROR;
  }
  double d;
  if (!ARG_IS_D(2,d))
  {
   Tcl_AppendResult(interp, "Need a ditance as 1st arg.", (char*) NULL);
   return TCL_ERROR;
  }
  int bond1,bond2,t;
  if ((!ARG_IS_I(3,bond1)) || (!ARG_IS_I(4,bond2)) || (!ARG_IS_I(5,t)))
  {
   Tcl_AppendResult(interp, "Need two bond types as 2nd and 3rd and a particle type as 4th argument.", (char*) NULL);
   return TCL_ERROR;
  }
  res=collision_detection_set_params(2,d,bond1,bond2,t);
  if (res==1)
  {
   Tcl_AppendResult(interp, "This mode requires the VIRTUAL_SITES_RELATIVE feature to be comiled in.", (char*) NULL);
   return TCL_ERROR;
  }
  if (res==2)
  {
   Tcl_AppendResult(interp, "Collision detection only works on a single cpu.", (char*) NULL);
   return TCL_ERROR;
  }
  if (res==3)
  {
   Tcl_AppendResult(interp, "A bond of the specified type does not exist.", (char*) NULL);
   return TCL_ERROR;
  }
 }
 else
 {
  Tcl_AppendResult(interp,"Unknown mode.",(char*)NULL);
  return TCL_ERROR;
 }
 return TCL_OK;
}

#endif
