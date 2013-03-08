/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
  Copyright (C) 2008,2009,2010 
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
/** \file adresso.c
    This is the place for adaptive resolution scheme
    Implementation of adresso.h
*/

#include "adresso_tcl.h"
#include "communication.h"
#include "parser.h"
#include "cells.h"
#include "adresso.h"
#include <tcl.h>
#include "virtual_sites.h"
#include "interaction_data.h"


/** \name Privat Functions */
/************************************************************/
/*@{*/
#ifdef ADRESS
/** prints adress settings */
int tclcommand_adress_parse_print(Tcl_Interp *interp,int argc, char **argv);

/** prints adress settings */
int tclcommand_adress_parse_set(Tcl_Interp *interp,int argc, char **argv);

#endif

/*@}*/

int tclcommand_adress(ClientData data, Tcl_Interp *interp, int argc, char **argv){
   int err = TCL_OK;
#ifndef ADRESS
   Tcl_ResetResult(interp);
   Tcl_AppendResult(interp, "Adress is not compiled in (change config.h).", (char *)NULL);
   err = (TCL_ERROR);
#else
   if (argc < 2) {
      Tcl_AppendResult(interp, "Wrong # of args! Usage: adress (set|print)", (char *)NULL);
      err = (TCL_ERROR);
   }
   else{
      if (ARG1_IS_S("print")) err=tclcommand_adress_parse_print(interp,argc,argv);
      else if (ARG1_IS_S("set")) err=tclcommand_adress_parse_set(interp,argc,argv);
      else {
         Tcl_ResetResult(interp);
         Tcl_AppendResult(interp, "The operation \"", argv[1],"\" you requested is not implemented.", (char *)NULL);
         err = (TCL_ERROR);
      }
   }
#endif
   return gather_runtime_errors(interp, err);
}

#ifdef ADRESS
int tclcommand_adress_parse_print(Tcl_Interp *interp,int argc, char **argv){
   int topo=(int)adress_vars[0],dim;
   char buffer[3*TCL_DOUBLE_SPACE];
   argv+=2;argc-=2;
   Tcl_ResetResult(interp);
   if (topo == 0) {
      Tcl_AppendResult(interp,"adress topo 0", (char *)NULL);
      return TCL_OK;
   }
   else if (topo == 1) {
      Tcl_PrintDouble(interp, adress_vars[1], buffer);
      Tcl_AppendResult(interp,"adress topo 1 width ",buffer, (char *)NULL);
      return TCL_OK;
   }
   //topo 2 and 3
   sprintf(buffer,"%i",topo);
   Tcl_AppendResult(interp,"adress topo ",buffer," width ",(char *)NULL);
   Tcl_PrintDouble(interp, adress_vars[1], buffer);
   Tcl_AppendResult(interp,buffer, " ", (char *)NULL);
   Tcl_PrintDouble(interp, adress_vars[2], buffer);
   Tcl_AppendResult(interp,buffer, " center ", (char *)NULL);
   
   if (topo==2) {
      dim=(int)adress_vars[3];
      if (dim==0) sprintf(buffer,"x");
      else if (dim==1) sprintf(buffer,"y");
      else sprintf(buffer,"z");
      Tcl_AppendResult(interp,buffer," ", (char *)NULL);
      Tcl_PrintDouble(interp, adress_vars[4], buffer);
   }
   else{ // topo == 3
      Tcl_PrintDouble(interp, adress_vars[3], buffer);
      Tcl_AppendResult(interp,buffer," ", (char *)NULL);
      Tcl_PrintDouble(interp, adress_vars[4], buffer);
      Tcl_AppendResult(interp,buffer," ", (char *)NULL);
      Tcl_PrintDouble(interp, adress_vars[5], buffer);
   }
   Tcl_AppendResult(interp,buffer, " wf ", (char *)NULL);
   sprintf(buffer,"%i",(int)adress_vars[6]);
   Tcl_AppendResult(interp,buffer, (char *)NULL);

   return TCL_OK;
}

int tclcommand_adress_parse_set(Tcl_Interp *interp,int argc, char **argv){
   int topo=-1,i,wf=0,set_center=0;
   double width[2],center[3];
   char buffer[3*TCL_DOUBLE_SPACE];
   argv+=2;argc-=2;

   for(i=0;i<3;i++) center[i]=box_l[i]/2;

   if (argc < 2) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp, "Wrong # of args! adress set needs at least 2 arguments\n", (char *)NULL);
      Tcl_AppendResult(interp, "Usage: adress set topo [0|1|2|3] width X.X Y.Y (center X.X Y.Y Z.Z) (wf [0|1])\n", (char *)NULL);
      Tcl_AppendResult(interp, "topo:   0 - switched off (no more values needed)\n", (char *)NULL);
      Tcl_AppendResult(interp, "        1 - constant (weight will be first value of width)\n", (char *)NULL);
      Tcl_AppendResult(interp, "        2 - divided in one direction (default x, or give a negative center coordinate\n", (char *)NULL);
      Tcl_AppendResult(interp, "        3 - spherical topology\n", (char *)NULL);
      Tcl_AppendResult(interp, "width:  X.X  - half of size of ex zone(r0/2 in the papers)\n", (char *)NULL);
      Tcl_AppendResult(interp, "        Y.Y  - size of hybrid zone (d in the papers)\n", (char *)NULL);
      Tcl_AppendResult(interp, "        Note: Only one value need for topo 1 \n", (char *)NULL);
      Tcl_AppendResult(interp, "center: center of the ex zone (default middle of the box) \n", (char *)NULL);
      Tcl_AppendResult(interp, "        Note: x|y|x X.X for topo 2  \n", (char *)NULL);
      Tcl_AppendResult(interp, "        Note: X.X Y.Y Z.Z for topo 3  \n", (char *)NULL);
      Tcl_AppendResult(interp, "wf:     0 - cos weighting function (default)\n", (char *)NULL);
      Tcl_AppendResult(interp, "        1 - polynom weighting function\n", (char *)NULL);
      Tcl_AppendResult(interp, "ALWAYS set box_l first !!!", (char *)NULL);
      return (TCL_ERROR);
   }

   //parse topo
   if ( (argc<2) || (!ARG0_IS_S("topo"))  || (!ARG1_IS_I(topo)) || (topo < 0) || (topo > 3) ) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp, "expected \'topo 0|1|2|3\'\n", (char *)NULL);
      return (TCL_ERROR);
   }
   argv+=2;argc-=2;
   
   //stop if topo is 0
   if (topo==0) {
      adress_vars[0]=0.0;
      mpi_bcast_parameter(FIELD_ADRESS);
      return TCL_OK;
   }

   //parse width
   if ( (argc>1) && (ARG0_IS_S("width")) ) {
      if (topo==1) {
         if ( (!ARG1_IS_D(width[0])) || (width[0]<0) ){
            Tcl_ResetResult(interp);
            Tcl_AppendResult(interp, "expected \'width X.X (X.X non-negative)\'", (char *)NULL);
            return (TCL_ERROR);
         }
         if ((width[0]> 1.0) || (width[0]< 0.0)) {
            Tcl_ResetResult(interp);
            Tcl_AppendResult(interp, "for constant topo, first width must be between 0 and 1", (char *)NULL);
            return (TCL_ERROR);
         }
         //stop if topo is 1
         adress_vars[0]=1;
         adress_vars[1]=width[0];
         mpi_bcast_parameter(FIELD_ADRESS);
         return TCL_OK;
      }
      else {//topo 2 and 3 are left over
         if ( (argc<3) || (!ARG1_IS_D(width[0])) || (width[0]<0) ||(!ARG_IS_D(2,width[1])) || (width[1]<0) ){
            Tcl_ResetResult(interp);
            Tcl_AppendResult(interp, "expected \'width X.X Y.Y (both non-negative)\'", (char *)NULL);
            return (TCL_ERROR);
         }
         argv+=3;argc-=3;
      }
   }
   else{
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp, "expected \'width\'", (char *)NULL);
      return (TCL_ERROR);
   }

   while (argc!=0){
      if (ARG0_IS_S("wf")){
         if ( (argc<2) || (!ARG1_IS_I(wf)) || (wf < 0) || (wf > 1) ){
            Tcl_ResetResult(interp);
            Tcl_AppendResult(interp, "expected \'wf 0|1\'", (char *)NULL);
            return (TCL_ERROR);
         }
         else{
            argv+=2;argc-=2;
         }
      }
      else if (ARG0_IS_S("center")){
         if (topo == 2) {
            if ( (argc<3) || ( (!ARG1_IS_S("x"))&&(!ARG1_IS_S("y"))&&(!ARG1_IS_S("z")) ) || (!ARG_IS_D(2,center[1])) ){
               Tcl_ResetResult(interp);
               Tcl_AppendResult(interp, "expected \'center x|y|z X.X\'", (char *)NULL);
               return (TCL_ERROR);
            }
            if (ARG1_IS_S("x")) center[0]=0;
            else if  (ARG1_IS_S("y")) center[0]=1;
            else center[0]=2;
            if ( (center[1]<0) || (center[1]>box_l[(int)center[0]]) ) {
               Tcl_ResetResult(interp);
               Tcl_AppendResult(interp, "The center component is outside the box", (char *)NULL);
               return (TCL_ERROR);
            }
            set_center=1;
            argv+=3;argc-=3;
         }
         else  { //topo 3
            if ( (argc<4) || (!ARG_IS_D(1,center[0])) || (!ARG_IS_D(2,center[1])) || (!ARG_IS_D(3,center[2])) ){
               Tcl_ResetResult(interp);
               Tcl_AppendResult(interp, "expected \'center X.X Y.Y Z.Z\'", (char *)NULL);
               return (TCL_ERROR);
            }
            argv+=4;argc-=4;
            //check components of center
            for (i=0;i<3;i++){
               if ( (center[i]<0)||(center[i]>box_l[i]) ){
                  Tcl_ResetResult(interp);
                  sprintf(buffer,"%i",i);
                  Tcl_AppendResult(interp, "The ",buffer," th component of center is outside the box\n", (char *)NULL);
                  return (TCL_ERROR);
               }
            }
         }
      }
      else{
         Tcl_ResetResult(interp);
         Tcl_AppendResult(interp, "The unknown operation \"", argv[0],"\".", (char *)NULL);
         return (TCL_ERROR);
      }
   }

   //set standard center value for topo 2
   if ((topo==2) && (set_center==0) ) center[0]=0;

   //width check
   if (topo==2){
      if (width[0]+width[1]>box_l[(int)center[0]]/2){
         Tcl_ResetResult(interp);
         Tcl_AppendResult(interp, "The width of ex+hy must smaller than box_l/2\n", (char *)NULL);
         return (TCL_ERROR);
      }
   }
   else if (topo==3){
      for (i=0;i<3;i++){
         if (width[0]+width[1]>box_l[i]/2){
            Tcl_ResetResult(interp);
            sprintf(buffer,"%i",i);
            Tcl_AppendResult(interp, "The width of ex+hy must smaller than box_l/2 in dim " ,buffer,"\n", (char *)NULL);
            return (TCL_ERROR);
         }
      }
   }

   adress_vars[0]=topo;
   adress_vars[1]=width[0];
   adress_vars[2]=width[1];
   adress_vars[3]=center[0];
   adress_vars[4]=center[1];
   adress_vars[5]=center[2];
   adress_vars[6]=wf;

   mpi_bcast_parameter(FIELD_ADRESS);

   return TCL_OK;
}

/* #ifdef THERMODYNAMIC_FORCE */
int tclcommand_thermodynamic_force(ClientData _data, Tcl_Interp * interp, int argc, char ** argv)
{
  int i, part_type, err_code;
  double j, prefactor;
  
  Tcl_ResetResult(interp);
  
  if(argc != 4){
    Tcl_AppendResult(interp, "wrong # args:  should be \"",
		     "thermodynamic_force <type> <filename> <prefactor>\"",
		     (char *) NULL);
    err_code = TCL_ERROR;
  }
  else {
    i=ARG_IS_I(1, part_type);
    j=ARG_IS_D(3,prefactor);
    if(i && j)
      err_code =  tclcommand_thermodynamic_force_parse_opt(interp, part_type, prefactor, argc-2, argv+2);
    else 
      err_code = TCL_ERROR;
  }
  
  return err_code;
}


int tclcommand_thermodynamic_force_parse_opt(Tcl_Interp * interp, int type, double prefactor, int argc, char ** argv){
  char * filename = NULL;
  
  filename = argv[0];
  
  switch(tf_set_params(type, prefactor, filename)){
    case 1:
    Tcl_AppendResult(interp, "particle type must be non-negative", (char *) NULL);
    return 0;
  case 2:
    Tcl_AppendResult(interp, "the length of the filename must be less than 256 characters,"
		     "but is \"", filename, "\"", (char *)NULL);
    return 0;
  case 3:
    Tcl_AppendResult(interp, "cannot open \"", filename, "\"", (char *)NULL);
    return 0;
  case 4:
    Tcl_AppendResult(interp, "attempt to read file \"", filename,
		     "\" failed, could not find start the start token <#>", (char *)NULL);
    return 0;
  case 5:
    Tcl_AppendResult(interp, "number of data points does not match the existing table", (char *)NULL);
    return 0;
  }
  return TCL_OK;
}
/* #endif */

int tclcommand_update_adress_weights(ClientData _data, Tcl_Interp * interp, int argc, char ** argv)
{
  int  err_code = TCL_OK;
  
  adress_update_weights();
  
  return err_code;
}

#ifdef INTERFACE_CORRECTION

int adress_tab_parser(Tcl_Interp * interp,
		      int part_type_a, int part_type_b,
		      int argc, char ** argv)
{
    char *filename = NULL;
    
    /* adress_tab interactions should supply a file name for a file containing
     both force and energy profiles as well as number of points, max
     values etc.
     */
    if (argc < 2) {
        Tcl_AppendResult(interp, "tabulated potentials require a filename: "
                         "<filename>",
                         (char *) NULL);
        return 0;
    }
    
    /* copy tabulated parameters */
    filename = argv[1];
    
    switch (adress_tab_set_params(part_type_a, part_type_b, filename)) {
        case 1:
            Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
            return 0;
        case 2:
            Tcl_AppendResult(interp, "the length of the filename must be less than 256 characters,"
                             "but is \"", filename, "\"", (char *)NULL);
            return 0;
        case 3:
            Tcl_AppendResult(interp, "cannot open \"", filename, "\"", (char *)NULL);
            return 0;
        case 4:
            Tcl_AppendResult(interp, "attempt to read file \"", filename,
                             "\" failed, could not find start the start token <#>", (char *)NULL);
            return 0;
    }
    return 2;
}

#endif

#endif
