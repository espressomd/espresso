/*
  Copyright (C) 2010,2011,2012,2012,2013 The ESPResSo project
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
#include "virtual_sites_com.h"
#include "virtual_sites_com_tcl.h"

#ifdef VIRTUAL_SITES_COM

#include "parser.h"
#include "virtual_sites.h"
#include "topology.h"
#include "grid.h"
#include "errorhandling.h"

int tclcommand_analyze_parse_and_print_pressure_mol(Tcl_Interp *interp,int argc, char **argv)
{
   char buffer[TCL_DOUBLE_SPACE];
   int type1, type2;
   double psum;
   #ifdef ELECTROSTATICS
   #ifndef INTER_RF
   Tcl_ResetResult(interp);
   Tcl_AppendResult(interp, "parse_and_print_pressure_mol is only possible with INTER_RF ", (char *)NULL);
   return (TCL_ERROR);
   #endif
   #endif
   updatePartCfg(WITHOUT_BONDS);
   if (!sortPartCfg()) {
      char *errtxt = runtime_error(128);
      ERROR_SPRINTF(errtxt, "{059 parse_and_print_pressure_mol: could not sort particle config, particle ids not consecutive?} ");
      return TCL_ERROR;
   }
   if (argc < 2) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp, "usage: analyze pressure_mol <type1> <type2>", (char *)NULL);
      return (TCL_ERROR);
   }

   if (!ARG0_IS_I(type1)) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp, "usage: analyze pressure_mol <type1> <type2>", (char *)NULL);
      return (TCL_ERROR);
   }
   if (!ARG1_IS_I(type2)) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp, "usage: analyze pressure_mol <type1> <type2>", (char *)NULL);
      return (TCL_ERROR);
   }
   argc-=2; argv+=2;

   if (n_molecules==0) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp, "No molecules defined !", (char *)NULL);
      return (TCL_ERROR);
   }
   psum=calc_pressure_mol(type1,type2);
   //sprintf(buffer,"%i",type1);
   //Tcl_AppendResult(interp,"{ analyze pressure_mol ",buffer," ",(char *)NULL);   
   //sprintf(buffer,"%i",type2);
   //Tcl_AppendResult(interp,buffer," ",(char *)NULL);   
   sprintf(buffer,"%e",psum);
   Tcl_AppendResult(interp, buffer,(char *)NULL);
   return TCL_OK;
}

int tclcommand_analyze_parse_and_print_energy_kinetic_mol(Tcl_Interp *interp,int argc, char **argv)
{
   char buffer[TCL_DOUBLE_SPACE];
   int type;
   double Ekin;
   updatePartCfg(WITHOUT_BONDS);
   if (!sortPartCfg()) {
      char *errtxt = runtime_error(128);
      ERROR_SPRINTF(errtxt, "{059 parse_and_print_energy_kinetic_mol: could not sort particle config, particle ids not consecutive?} ");
      return TCL_ERROR;
   }
   if (argc < 1) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp, "usage: analyze energy_kinetic <type>", (char *)NULL);
      return (TCL_ERROR);
   }

   if (!ARG0_IS_I(type)) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp, "usage: analyze energy_kinetic <type>", (char *)NULL);
      return (TCL_ERROR);
   }
   argc-=1; argv+=1;

   if (n_molecules==0) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp, "No molecules defined !", (char *)NULL);
      return (TCL_ERROR);
   }
   Ekin=calc_energy_kinetic_mol(type);
   if (Ekin < 0.0) {
     Tcl_ResetResult(interp);
      sprintf(buffer,"%i",-(int)Ekin);
      Tcl_AppendResult(interp,"Could not fetch com in calc_energy_kinetic_mol! From mol_id",buffer, (char *)NULL);
      return (TCL_ERROR);
   }
   //sprintf(buffer,"%i",type);
   //Tcl_AppendResult(interp,"{ analyze pressure_mol ",buffer," ",(char *)NULL);   
   sprintf(buffer,"%e",Ekin);
   Tcl_AppendResult(interp, buffer,(char *)NULL);
   return TCL_OK;
}

int tclcommand_analyze_parse_and_print_dipmom_mol(Tcl_Interp *interp,int argc, char **argv)
{
#ifndef ELECTROSTATICS
   Tcl_ResetResult(interp);
   Tcl_AppendResult(interp, "calc_dipole_mol is not possible without ELECTROSTATICS", (char *)NULL);
   return (TCL_ERROR);
#else
   int k,type;
   char buffer[TCL_DOUBLE_SPACE];
   double dipole[4];
   updatePartCfg(WITHOUT_BONDS);
   if (!sortPartCfg()) {
      char *errtxt = runtime_error(128);
      ERROR_SPRINTF(errtxt, "{059 parse_and_print_dipole: could not sort particle config, particle ids not consecutive?} ");
      return TCL_ERROR;
   }
   if (n_molecules==0) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp, "No molecules defined !", (char *)NULL);
      return (TCL_ERROR);
   }
   if (argc < 2) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp, "usage: analyze parse_and_print_dipole_mol <type>", (char *)NULL);
      return (TCL_ERROR);
   }

   if (!ARG1_IS_I(type)) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp, "usage: analyze parse_and_print_dipole_mol <type>", (char *)NULL);
      return (TCL_ERROR);
   }
   if (ARG0_IS_S("total")){
      calc_total_dipolmoment_mol(type,dipole);
      sprintf(buffer,"%i ",type);
      Tcl_AppendResult(interp,"{ dipolemoment_mol total ",buffer,(char *)NULL);
      for (k=0;k<3;k++)
      {
            sprintf(buffer,"%e ",dipole[k]);
            Tcl_AppendResult(interp, buffer,(char *)NULL);
      }
      sprintf(buffer,"%e",dipole[3]);
      Tcl_AppendResult(interp,buffer,"}",(char *)NULL);
   }
   else if (ARG0_IS_S("absolute")){
      calc_absolute_dipolmoment_mol(type,dipole);
      sprintf(buffer,"%i ",type);
      Tcl_AppendResult(interp,"{ dipolemoment_mol absolute ",buffer,(char *)NULL);
      sprintf(buffer,"%e ",dipole[0]);
      Tcl_AppendResult(interp, buffer,(char *)NULL);
      sprintf(buffer,"%e",dipole[1]);
      Tcl_AppendResult(interp,buffer,"}",(char *)NULL);
   }
   else
   {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp, "Feature not implemented", (char *)NULL);
      return (TCL_ERROR);
   }
   argc-=2; argv+=2;
   return TCL_OK;
#endif
}

int tclcommand_analyze_parse_and_print_check_mol(Tcl_Interp *interp,int argc, char **argv){
   int j,count=0;
   double dist;
   char buffer[TCL_DOUBLE_SPACE];
   Particle p;
   updatePartCfg(WITHOUT_BONDS);
   for(j=0; j<n_total_particles; j++){
      if (!ifParticleIsVirtual(&partCfg[j])) continue;
      get_particle_data(j,&p);
      //dist=min_distance(partCfg[j].r.p,p.r.p);
      unfold_position(p.r.p,p.l.i);
      dist=distance(partCfg[j].r.p,p.r.p);
      if (dist > 0.01){
         if (count==0) Tcl_AppendResult(interp,"BEGIN Particle Missmatch: \n", (char *)NULL);
         sprintf(buffer,"%i",j);
         Tcl_AppendResult(interp,"Particle ",buffer, (char *)NULL);
         Tcl_PrintDouble(interp,partCfg[j].r.p[0] , buffer);
         Tcl_AppendResult(interp," partCfg x ",buffer, (char *)NULL);
         Tcl_PrintDouble(interp,partCfg[j].r.p[1] , buffer);
         Tcl_AppendResult(interp," y ",buffer, (char *)NULL);
         Tcl_PrintDouble(interp,partCfg[j].r.p[2] , buffer);
         Tcl_AppendResult(interp," z ",buffer, (char *)NULL);
         Tcl_PrintDouble(interp,p.r.p[0] , buffer);
         Tcl_AppendResult(interp," my_partCfg x ",buffer, (char *)NULL);
         Tcl_PrintDouble(interp,p.r.p[1] , buffer);
         Tcl_AppendResult(interp," y ",buffer, (char *)NULL);
         Tcl_PrintDouble(interp,p.r.p[2] , buffer);
         Tcl_AppendResult(interp," z ",buffer, (char *)NULL);
         Tcl_PrintDouble(interp, dist, buffer);
         Tcl_AppendResult(interp," dist ",buffer,"\n", (char *)NULL);
         count++;
      }
   }
   if (count!=0){
      Tcl_AppendResult(interp,"END Particle Missmatch\n", (char *)NULL);
      return(TCL_ERROR);
   }
   return(TCL_OK);
}

#endif
