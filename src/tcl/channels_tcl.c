/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
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

#include <unistd.h>
#include <string.h>
#include <tcl.h>
#include "utils.h"

int tclcommand_replacestdchannel(ClientData clientData, Tcl_Interp *interp, int argc, char **argv)
{
  Tcl_Channel channel=NULL;

  if (argc != 3){
    Tcl_AppendResult(interp,"Wrong # of args! Usage: ",argv[0]," [stdout|stderr|stdin] <pipename>",(char *)NULL);
    return TCL_ERROR;
  }

  if(access(argv[2],F_OK)<0){
    Tcl_AppendResult(interp,"File ",argv[2]," does not exist!",(char *) NULL);
    return TCL_ERROR;
  }

  if(strcmp(argv[1],"stdout")==0){
    if(access(argv[2],W_OK)<0){
      Tcl_AppendResult(interp,"You do not have permission to access ",argv[2],(char *) NULL);
      return TCL_ERROR;
    }
    Tcl_UnregisterChannel(interp,Tcl_GetStdChannel(TCL_STDOUT));
    channel = Tcl_OpenFileChannel(interp, argv[2], "WRONLY",0666);
    Tcl_RegisterChannel(interp,channel);
    Tcl_SetStdChannel(channel,TCL_STDOUT);
  }
  else if(strcmp(argv[1],"stderr")==0){
    if(access(argv[2],W_OK)<0){
      Tcl_AppendResult(interp,"You do not have permission to access ",argv[2],(char *) NULL);
      return TCL_ERROR;
    }
    Tcl_UnregisterChannel(interp,Tcl_GetStdChannel(TCL_STDERR));
    channel = Tcl_OpenFileChannel(interp, argv[2], "WRONLY",0666);
    Tcl_RegisterChannel(interp,channel);
    Tcl_SetStdChannel(channel,TCL_STDERR);
  }
  else if(strcmp(argv[1],"stdin")==0){
    if(access(argv[2],R_OK)<0){
      Tcl_AppendResult(interp,"You do not have permission to access ",argv[2],(char *) NULL);
      return TCL_ERROR;
    }
    Tcl_UnregisterChannel(interp,Tcl_GetStdChannel(TCL_STDIN));
    channel = Tcl_OpenFileChannel(interp, argv[2], "RDONLY",0666);
    Tcl_RegisterChannel(interp,channel);
    Tcl_SetStdChannel(channel,TCL_STDIN);
  }
  else{
    Tcl_AppendResult(interp,"invalid first argument (got: ",argv[1]," )",(char *) NULL);
    return TCL_ERROR;
  }

  if(channel == NULL)
    return TCL_ERROR;

  return TCL_OK;
}
