// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2009; all rights reserved unless otherwise stated.
/** \file adresso.c
    This is the place for adaptive resolution scheme
    Implementation of adresso.h
*/

#include "adresso.h"
#include "communication.h"
#include "parser.h"
#include "cells.h"

/** \name Privat Functions */
/************************************************************/
/*@{*/
#ifdef ADRESS
/** prints adress settings */
int adress_print(Tcl_Interp *interp,int argc, char **argv);

/** prints adress settings */
int adress_set(Tcl_Interp *interp,int argc, char **argv);

/** calc weighting function of a distance
    @param dist distance
    @return weight of the distance
*/
double adress_wf(double dist);

#endif

/*@}*/

double adress_vars[7]       = {0, 0, 0, 0, 0, 0, 0};

int adress_tcl(ClientData data, Tcl_Interp *interp, int argc, char **argv){
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
      if (ARG1_IS_S("print")) err=adress_print(interp,argc,argv);
      else if (ARG1_IS_S("set")) err=adress_set(interp,argc,argv);
      else {
         Tcl_ResetResult(interp);
         Tcl_AppendResult(interp, "The operation \"", argv[1],"\" you requested is not implemented.", (char *)NULL);
         err = (TCL_ERROR);
      }
   }
#endif
   return mpi_gather_runtime_errors(interp, err);
}

#ifdef ADRESS
int adress_print(Tcl_Interp *interp,int argc, char **argv){
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

int adress_set(Tcl_Interp *interp,int argc, char **argv){
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

double adress_wf_vector(double x[3]){
  int topo=(int)adress_vars[0];
  double dist;
  int dim;
  
  int img_box[3];
  double temp_pos[3];
  
  switch (topo) {
  case 0:
    return 0.0;
    break;
  case 1:
    return adress_vars[1];
    break;
  case 2:
    dim=(int)adress_vars[3];
    //dist=fabs(x[dim]-adress_vars[4]);
    dist = x[dim]-adress_vars[4];
    if(dist>0)
      while(dist>box_l[dim]/2.0)
	dist = dist - box_l[dim];
    else if(dist < 0)
      while(dist< -box_l[dim]/2.0)
	dist = dist + box_l[dim];
    dist = fabs(dist);
    return adress_wf(dist);
    break;
  case 3:
    for(dim=0;dim<3;dim++){
      img_box[dim]=0;
      temp_pos[dim]=x[dim];
    }
    fold_position(temp_pos,img_box);
    dist=distance(temp_pos,&(adress_vars[3]));
    //printf("%f %f \n", dist, adress_wf(dist));
    return adress_wf(dist);
    break;
  default:
    return 0.0;
    break;
  }
}

double adress_wf(double dist){
   int wf;
   double tmp;
   
   //explicit region
   if (dist < adress_vars[1]) return 1;
   //cg regime
   else if (dist> adress_vars[1]+adress_vars[2]) return 0;
   else {
      wf=(int)adress_vars[6];
      if (wf == 0){ //cos
         tmp=PI/2/adress_vars[2]*(dist-adress_vars[1]);
         return cos(tmp)*cos(tmp);
      }
      else{ //wf == 1
         tmp=(dist-adress_vars[1]);
         return 1+2*tmp*tmp-3*tmp*tmp*tmp;
      }
   }
}

void adress_update_weights(){
  Particle *p;
  int i, np, c;
  Cell *cell;
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      if (ifParticleIsVirtual(&p[i])) {
         p[i].p.adress_weight=adress_wf_vector((&p[i])->r.p);
	 //printf("LOCAL %f %f\n", p[i].r.p[0], p[i].p.adress_weight);
      }
    }
  }
  for (c = 0; c < local_cells.n; c++) {
    cell = ghost_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      if (ifParticleIsVirtual(&p[i])) {
         p[i].p.adress_weight=adress_wf_vector((&p[i])->r.p);
	 //printf("GHOST %f %f\n", p[i].r.p[0], p[i].p.adress_weight);
      }
    }
  }
}

/** #ifdef THERMODYNAMIC_FORCE */
int tf_tcl(ClientData _data, Tcl_Interp * interp, int argc, char ** argv)
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
      err_code =  tf_parse(interp, part_type, prefactor, argc-2, argv+2);
    else 
      err_code = TCL_ERROR;
  }
  
  return err_code;
}


int tf_parse(Tcl_Interp * interp, int type, double prefactor, int argc, char ** argv){
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
/** #endif */

int manual_update_weights(ClientData _data, Tcl_Interp * interp, int argc, char ** argv)
{
  int  err_code = TCL_OK;
  
  adress_update_weights();
  
  return err_code;
}
#endif
