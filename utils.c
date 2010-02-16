// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2009; all rights reserved unless otherwise stated.

#include <unistd.h>
#include <string.h>
#include <tcl.h>
#include "utils.h"



int replacestdchannel(ClientData clientData, Tcl_Interp *interp, int argc, char **argv)
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
