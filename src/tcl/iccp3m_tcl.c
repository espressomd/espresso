/*
  Copyright (C) 2010,2011 The ESPResSo project
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

/** \file iccp3m.c
    Detailed Information about the method is included in the corresponding header file \ref iccp3m.h.

 */

#include "iccp3m.h"

#ifdef ELECTROSTATICS
enum { ICCP3M_AREA , ICCP3M_EPSILON, ICCP3M_NORMAL, ICCP3M_EXTFIELD } ;
static int tclcommand_iccp3m_parse_params(Tcl_Interp *interp,int normal_args, char *string, int flag);

/** Parses the ICCP3M command.
 */
int tclcommand_iccp3m(ClientData data, Tcl_Interp *interp, int argc, char **argv) {
  int last_ind_id,num_iteration,normal_args,area_args;
  char buffer[TCL_DOUBLE_SPACE];
  double e1,convergence,relax;

  Tcl_AppendResult(interp, "The ICCP3M algorithm is still experimental. Function can not be guaranteed, therefore it is still disabled.\n", (char *)NULL);
  return (TCL_ERROR);

  if(iccp3m_initialized==0){
      iccp3m_init();
      iccp3m_initialized=1;
  }

  if(argc != 9 && argc != 2 && argc != 10) { 
         Tcl_AppendResult(interp, "Wrong # of args! Usage: iccp3m { iterate | <last_ind_id> <e1> <num_iteration> <convergence> <relaxation> <area> <normal_components> <e_in/e_out>  [<ext_field>] }", (char *)NULL); 
         return (TCL_ERROR); 
   }
   if (argc == 2 ){
      if(ARG_IS_S(1,"iterate")) { 
           if (iccp3m_cfg.set_flag==0) {
                 Tcl_AppendResult(interp, "iccp3m parameters not set!", (char *)NULL);
                 return (TCL_ERROR);
           }
           else{ 
              Tcl_PrintDouble(interp,mpi_iccp3m_iteration(0),buffer); 
              Tcl_AppendResult(interp, buffer, (char *) NULL);
              return TCL_OK;
	   }
      }
   }
   else {
      if(!ARG_IS_I(1, last_ind_id)) {
          Tcl_ResetResult(interp);
          Tcl_AppendResult(interp, "Last induced id must be integer (got: ", argv[1],")!", (char *)NULL); return (TCL_ERROR);
          } else if(last_ind_id < 2) { 
          Tcl_ResetResult(interp);
          Tcl_AppendResult(interp, "Last induced id can not be smaller then 2 (got: ", argv[1],")!", (char *)NULL); return (TCL_ERROR);
       }
       if(!ARG_IS_D(2, e1)) {
          Tcl_ResetResult(interp);
          Tcl_AppendResult(interp, "Dielectric constant e1(inner) must be double(got: ", argv[2],")!", (char *)NULL); return (TCL_ERROR);
       }
       if(!ARG_IS_I(3, num_iteration)) {
          Tcl_ResetResult(interp);
          Tcl_AppendResult(interp, "Number of maximum iterations must be integer (got: ", argv[4],")!", (char *)NULL); return (TCL_ERROR);
       }
       if(!ARG_IS_D(4, convergence)) {
          Tcl_ResetResult(interp);
          Tcl_AppendResult(interp, "Convergence criterion must be double (got: ", argv[5],")!", (char *)NULL); return (TCL_ERROR);
          } else if( (convergence < 1e-10) || (convergence > 1e-1) ) { 
            Tcl_ResetResult(interp);
            Tcl_AppendResult(interp, "Convergence criterion can not be smaller then 1e-10 and greater then 1e-2(got: ", argv[5],")!", (char *)NULL); return (TCL_ERROR);
         }
          if(!ARG_IS_D(5, relax)) {
             Tcl_ResetResult(interp);
             Tcl_AppendResult(interp, "Relaxation parameter must be double (got: ", argv[6],")!", (char *)NULL); return (TCL_ERROR);
       }
   
       iccp3m_cfg.last_ind_id = last_ind_id;      /* Assign needed options */
       iccp3m_cfg.num_iteration = num_iteration; /* maximum number of iterations to go */
       iccp3m_cfg.eout = e1;
       iccp3m_cfg.convergence = convergence;
       iccp3m_cfg.relax = relax;
       iccp3m_cfg.update = 0;
       iccp3m_cfg.set_flag = 1;
       
       normal_args = (iccp3m_cfg.last_ind_id+1)*3;
       /* Now get Normal Vectors components consecutively */
       
       if( tclcommand_iccp3m_parse_params(interp,normal_args,argv[7],ICCP3M_NORMAL) == 1) { 
              Tcl_ResetResult(interp);
              Tcl_AppendResult(interp, "ICCP3M: Error in following normal vectors\n", argv[7],"\nICCP3M: Error in previous normal vectors\n", (char *)NULL); 
              return (TCL_ERROR);
       }      

       area_args=(iccp3m_cfg.last_ind_id+1);
       /* Now get area of the boundary elements */
       
       if ( tclcommand_iccp3m_parse_params(interp,area_args,argv[6],ICCP3M_AREA) == 1 ){
             Tcl_ResetResult(interp);
             Tcl_AppendResult(interp, "ICCP3M: Error in following areas\n", argv[6],"\nICCP3M: Error in previous areas\n", (char *)NULL); return (TCL_ERROR);
       }
       /* Now get the ration ein/eout for each element. 
          It's the user's duty to make sure that only disconnected 
          regions have different ratios */
      if ( tclcommand_iccp3m_parse_params(interp,area_args,argv[8],ICCP3M_EPSILON) == 1 ) {
             Tcl_ResetResult(interp);
             Tcl_AppendResult(interp, "ICCP3M: Error in following dielectric constants\n", argv[8],"\nICCP3M:  Error in previous dielectric constants\n", (char *)NULL); return (TCL_ERROR);
       } 

       if( argc == 10 ) {
         if( tclcommand_iccp3m_parse_params(interp,normal_args,argv[9],ICCP3M_EXTFIELD) == 1) { 
              Tcl_ResetResult(interp);
              Tcl_AppendResult(interp, "ICCP3M: Error in following external field vectors\n", argv[9],"\nICCP3M: Error in previous external field vectors\n", (char *)NULL); return (TCL_ERROR);
         }      
       }
       else {
         printf("allocating zeroes for external field \n");
         iccp3m_cfg.extx = (double*)calloc((last_ind_id +1), sizeof(double));
         iccp3m_cfg.exty = (double*)calloc((last_ind_id +1), sizeof(double));
         iccp3m_cfg.extz = (double*)calloc((last_ind_id +1), sizeof(double));
       }
      
       mpi_iccp3m_init(0);
       Tcl_PrintDouble(interp,mpi_iccp3m_iteration(0),buffer); 
       Tcl_AppendResult(interp, buffer, (char *) NULL);
       return TCL_OK;
   } /* else (argc==10) */
   return TCL_OK;
}

static int tclcommand_iccp3m_parse_params(Tcl_Interp *interp,int normal_args, char *string, int flag) {
  /* This function parses a vector give as C-String */
  /* It currently is not too elegant. But works. */
  int size,i,k=0;
  int scan_succes;
  static double *numbers=NULL;
  const char delimiters[] = " ";
  char *token,*cp;
  float temp;

  size= normal_args;
  numbers = malloc((size)*sizeof(double));

  cp = strdup(string);                /* Make writable copy.  */
  token = strtok (cp, delimiters);
  scan_succes = sscanf(token,"%f",&temp);
  if (scan_succes < 1) 
    return 1;

  numbers[0]=temp;
    
  for(i=1;i<size;i++) {
    token = strtok (NULL, delimiters);
    if(token == NULL)
      return 1;
    scan_succes = sscanf(token,"%lf",&(numbers[i]));
    if (scan_succes < 1 ) 
      return 1;
  }

  switch(flag) {
    case ICCP3M_AREA: 
      iccp3m_cfg.areas = (double*) realloc(iccp3m_cfg.areas, (size+1)*sizeof(double)); 
      for( i = 0 ; i < size ; i++ )  
        iccp3m_cfg.areas[i]=numbers[i];
      break;
    case ICCP3M_EPSILON:
      iccp3m_cfg.ein = (double*) realloc(iccp3m_cfg.ein,(size+1)*sizeof(double));
      for( i = 0 ; i < size; i++)  
        iccp3m_cfg.ein[i]=numbers[i];
    break;
    case ICCP3M_NORMAL:
      iccp3m_cfg.nvectorx = (double*) realloc(iccp3m_cfg.nvectorx,sizeof(double)*(iccp3m_cfg.last_ind_id+1));
      iccp3m_cfg.nvectory = (double*) realloc(iccp3m_cfg.nvectory,sizeof(double)*(iccp3m_cfg.last_ind_id+1));
      iccp3m_cfg.nvectorz = (double*) realloc(iccp3m_cfg.nvectorz,sizeof(double)*(iccp3m_cfg.last_ind_id+1));
      for(i=0;i<size;i++) {
        if( i%3 == 0 ) { iccp3m_cfg.nvectorx[k] = numbers[i]; } 
        if( i%3 == 1 ) { iccp3m_cfg.nvectory[k] = numbers[i]; }
        if( i%3 == 2 ) { iccp3m_cfg.nvectorz[k] = numbers[i];  k++; } 
       }
    break;

    case ICCP3M_EXTFIELD:
      iccp3m_cfg.extx = (double*) realloc(iccp3m_cfg.extx,sizeof(double)*(iccp3m_cfg.last_ind_id+1));
      iccp3m_cfg.exty = (double*) realloc(iccp3m_cfg.exty,sizeof(double)*(iccp3m_cfg.last_ind_id+1));
      iccp3m_cfg.extz = (double*) realloc(iccp3m_cfg.extz,sizeof(double)*(iccp3m_cfg.last_ind_id+1));
      for(i=0;i<size;i++) {
        if( i%3 == 0 ) { iccp3m_cfg.extx[k] = numbers[i]; } 
        if( i%3 == 1 ) { iccp3m_cfg.exty[k] = numbers[i]; }
        if( i%3 == 2 ) { iccp3m_cfg.extz[k] = numbers[i];  k++; } 
      }
    break;
  }

  free(numbers);
  return (0);
}

#endif

