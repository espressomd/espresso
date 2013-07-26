/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

/** \file reaction_tcl.c
 *
*/

#include "galilei.hpp"
#include "parser.hpp"
#include "utils.hpp"
#include "grid.hpp"
#include "communication.hpp"
#include "integrate.hpp"

/* ############### */

/* Stop motion of the particles */
int tcl_command_kill_particle_motion_print_usage(Tcl_Interp * interp){
  char buffer[256];
  sprintf(buffer, "Usage: kill_particle_motion [rotation]\n");
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return TCL_ERROR;
}

/* Set the forces on the particles to zero */
int tcl_command_kill_particle_forces_print_usage(Tcl_Interp * interp){
  char buffer[256];
  sprintf(buffer, "Usage: kill_particle_forces [torques]\n");
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return TCL_ERROR;
}

/* Calculate the CMS of the system */
int tcl_command_system_CMS_print_usage(Tcl_Interp * interp){
  char buffer[256];
  sprintf(buffer, "Usage: system_CMS [folded]\n");
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return TCL_ERROR;
}

/* Calculate the CMS velocity of the system */
int tcl_command_system_CMS_velocity_print_usage(Tcl_Interp * interp){
  char buffer[256];
  sprintf(buffer, "Usage: system_CMS_velocity\n");
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return TCL_ERROR;
}

/* Remove the CMS velocity of the system */
int tcl_command_galilei_transform_print_usage(Tcl_Interp * interp){
  char buffer[256];
  sprintf(buffer, "Usage: galilei_transform\n");
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return TCL_ERROR;
}

/* ############### */

/* Stop motion of the particles */
int tclcommand_kill_particle_motion(ClientData data, Tcl_Interp * interp, int argc, char ** argv){
  if (argc != 1 && argc != 2 ) { 
    return tcl_command_kill_particle_motion_print_usage(interp);
  } else {
    if (ARG_IS_S_EXACT(0,"kill_particle_motion")) {
      if ( argc == 2 ) {
        if (ARG_IS_S_EXACT(1,"rotation")) {
          mpi_kill_particle_motion( 1 );
          return TCL_OK;
        } else {
          return tcl_command_kill_particle_motion_print_usage(interp);
        }
      } else {
        mpi_kill_particle_motion( 0 );
        return TCL_OK;
      }
    } else {
      return tcl_command_kill_particle_motion_print_usage(interp);
    }
  }
}

/* Set the forces on the particles to zero */
int tclcommand_kill_particle_forces(ClientData data, Tcl_Interp * interp, int argc, char ** argv){
  if (argc != 1 && argc != 2 ) { 
    return tcl_command_kill_particle_forces_print_usage(interp);
  } else {
    if (ARG_IS_S_EXACT(0,"kill_particle_forces")) {
      if ( argc == 2 ) {
        if (ARG_IS_S_EXACT(1,"torques")) {
          mpi_kill_particle_forces( 1 );
          return TCL_OK;
        } else {
          return tcl_command_kill_particle_forces_print_usage(interp);
        }
      } else {
        mpi_kill_particle_forces( 0 );
        return TCL_OK;
      }
    } else {
      return tcl_command_kill_particle_forces_print_usage(interp);
    }
  }
}

/* Calculate the CMS of the system */
int tclcommand_system_CMS(ClientData data, Tcl_Interp * interp, int argc, char ** argv){
  char buffer[256];
  double cmspos[3];
  int box[3];

  if (argc != 1 && argc != 2  ) { 
    return tcl_command_system_CMS_print_usage(interp);
  } else {
    if (ARG_IS_S_EXACT(0,"system_CMS")) {
      if ( argc == 2 ) {
        if (ARG_IS_S_EXACT(1,"folded")) {
          mpi_system_CMS();

          memcpy(cmspos, gal.cms, 3*sizeof(double));
          box[0] = 0; box[1] = 0; box[2] = 0;
          fold_position(cmspos, box);

          Tcl_PrintDouble(interp, cmspos[0], buffer);
          Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
          Tcl_PrintDouble(interp, cmspos[1], buffer);
          Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
          Tcl_PrintDouble(interp, cmspos[2], buffer);
          Tcl_AppendResult(interp, buffer, (char *)NULL);

          return TCL_OK;
        } else {
          return tcl_command_system_CMS_print_usage(interp);
        }
      } else {
        mpi_system_CMS();

        Tcl_PrintDouble(interp, gal.cms[0], buffer);
        Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
        Tcl_PrintDouble(interp, gal.cms[1], buffer);
        Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
        Tcl_PrintDouble(interp, gal.cms[2], buffer);
        Tcl_AppendResult(interp, buffer, (char *)NULL);


        return TCL_OK;     
      }
    } else {
      return tcl_command_system_CMS_print_usage(interp);
    }
  }
}

/* Calculate the CMS velocity of the system */
int tclcommand_system_CMS_velocity(ClientData data, Tcl_Interp * interp, int argc, char ** argv){
  char buffer[256];

  if (argc != 1  ) { 
    return tcl_command_system_CMS_velocity_print_usage(interp);
  } else {
    if (ARG_IS_S_EXACT(0,"system_CMS_velocity")) {
        mpi_system_CMS_velocity();

        Tcl_PrintDouble(interp, gal.cms_vel[0]/time_step, buffer);
        Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
        Tcl_PrintDouble(interp, gal.cms_vel[1]/time_step, buffer);
        Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
        Tcl_PrintDouble(interp, gal.cms_vel[2]/time_step, buffer);
        Tcl_AppendResult(interp, buffer, (char *)NULL);

        return TCL_OK;
    } else {
      return tcl_command_system_CMS_velocity_print_usage(interp);
    }
  }
}

/* Remove the CMS velocity of the system */
int tclcommand_galilei_transform(ClientData data, Tcl_Interp * interp, int argc, char ** argv){
  if (argc != 1  ) { 
    return tcl_command_galilei_transform_print_usage(interp);
  } else {
    if (ARG_IS_S_EXACT(0,"galilei_transform")) {
        mpi_galilei_transform();
        return TCL_OK;
    } else {
      return tcl_command_galilei_transform_print_usage(interp);
    }
  }
}
