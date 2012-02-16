/*
  Copyright (C) 2010,2011,2012 The ESPResSo project
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
#include "parser.h"

#ifdef ELECTROSTATICS
enum { ICCP3M_AREA , ICCP3M_EPSILON, ICCP3M_NORMAL, ICCP3M_SIGMA, ICCP3M_EXTFIELD } ;
int tclcommand_iccp3m_parse_normals(Tcl_Interp *interp,int n_ic, char *string);
int tclcommand_iccp3m_parse_double_list(Tcl_Interp *interp, int n_ic, char *string, int flag); 

/** Parses the ICCP3M command.
 */
int tclcommand_iccp3m(ClientData data, Tcl_Interp *interp, int argc, char **argv) {
  char buffer[TCL_DOUBLE_SPACE];

  if(iccp3m_initialized==0){
      iccp3m_init();
      iccp3m_initialized=1;
  }

    iccp3m_cfg.num_iteration=30;
  iccp3m_cfg.convergence=1e-2;
  iccp3m_cfg.relax=0.7;
  iccp3m_cfg.eout=1;


  if(argc < 2 ) { 
         Tcl_AppendResult(interp, "Usage of ICCP3M: RTFM", (char *)NULL); 
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
     if(ARG_IS_I(1, iccp3m_cfg.n_ic)) {
       argc-=2;
       argv+=2;
     } else {
       Tcl_AppendResult(interp, "ICCP3M: First argument has to be the number of induced charges", (char *)NULL); 
       return (TCL_ERROR);
     }
     while (argc > 0) {
       if (ARG0_IS_S("convergence")) {
         if (argc>1 && ARG1_IS_D(iccp3m_cfg.convergence)) {
           argc-=2;
           argv+=2;
         } else {
           Tcl_AppendResult(interp, "ICCP3M Usage: convergence <convergence>", (char *)NULL); 
           return (TCL_ERROR);
         }
       } else if (ARG0_IS_S("relaxation")) {
         if (argc>1 && ARG1_IS_D(iccp3m_cfg.relax)) {
           argc-=2;
           argv+=2;
         } else {
           Tcl_AppendResult(interp, "ICCP3M Usage: convergence <convergence>", (char *)NULL); 
           return (TCL_ERROR);
         }
       } else if (ARG0_IS_S("eps_out")) {
         if (argc>1 && ARG1_IS_D(iccp3m_cfg.eout)) {
           argc-=2;
           argv+=2;
         } else {
           Tcl_AppendResult(interp, "ICCP3M Usage: eps_out <eps_out>", (char *)NULL); 
           return (TCL_ERROR);
         }
       } else if (ARG0_IS_S("max_iterations")) {
         if (argc>1 && ARG1_IS_I(iccp3m_cfg.num_iteration)) {
           argc-=2;
           argv+=2;
         } else {
           Tcl_AppendResult(interp, "ICCP3M Usage: max_iterations <max_iterations>", (char *)NULL); 
           return (TCL_ERROR);
         }
       } else if (ARG0_IS_S("normals")) {
         if (argc>1) {
           if (tclcommand_iccp3m_parse_normals(interp, iccp3m_cfg.n_ic, argv[1]) != TCL_OK) {
             return TCL_ERROR;
           }
           argc-=2;
           argv+=2;
         } else {
           Tcl_AppendResult(interp, "ICCP3M Usage: normals <List of normal vectors>", (char *)NULL); 
           return (TCL_ERROR);
         }
       } else if (ARG0_IS_S("areas")) {
         if (argc>1) {
           if (tclcommand_iccp3m_parse_double_list(interp, iccp3m_cfg.n_ic, argv[1], ICCP3M_AREA)!=TCL_OK) {
             return TCL_ERROR;
           }
           argc-=2;
           argv+=2;
         } else {
           Tcl_AppendResult(interp, "ICCP3M Usage: areas <list of areas>", (char *)NULL); 
           return (TCL_ERROR);
         }
       } else if (ARG0_IS_S("sigmas")) {
         if (argc>1) {
           if (tclcommand_iccp3m_parse_double_list(interp, iccp3m_cfg.n_ic, argv[1], ICCP3M_SIGMA)!=TCL_OK) {
             return TCL_ERROR;
           }
           argc-=2;
           argv+=2;
         } else {
           Tcl_AppendResult(interp, "ICCP3M Usage: sigmas <list of sigmas>", (char *)NULL); 
           return (TCL_ERROR);
         }
       } else if (ARG0_IS_S("epsilons")) {
         if (argc>1) {
           if (tclcommand_iccp3m_parse_double_list(interp, iccp3m_cfg.n_ic, argv[1], ICCP3M_EPSILON) != TCL_OK) {
             return TCL_ERROR;
           }
           argc-=2;
           argv+=2;
         } else {
           Tcl_AppendResult(interp, "ICCP3M Usage: epsilons <list of epsilons>", (char *)NULL); 
           return (TCL_ERROR);
         }
       } else {
         Tcl_AppendResult(interp, "Unknown Argument to ICCP3M ", argv[0], (char *)NULL); 
         return (TCL_ERROR);
       }
     }
   }
   iccp3m_initialized=1;
   iccp3m_cfg.set_flag = 1;
      
      
   mpi_iccp3m_init(0);

   return TCL_OK;
}

int tclcommand_iccp3m_parse_normals(Tcl_Interp *interp,int n_ic, char *string) {
  char *arg, *token;
  int scan_succes;
  arg=strdup(string);
  iccp3m_cfg.nvectorx = (double*) realloc(iccp3m_cfg.nvectorx,sizeof(double)*(iccp3m_cfg.n_ic));
  iccp3m_cfg.nvectory = (double*) realloc(iccp3m_cfg.nvectory,sizeof(double)*(iccp3m_cfg.n_ic));
  iccp3m_cfg.nvectorz = (double*) realloc(iccp3m_cfg.nvectorz,sizeof(double)*(iccp3m_cfg.n_ic));
  const char opening_bracket[] = "{";
  const char closing_bracket[] = "}";
  const char space[] = " ";
  

  // Searching for first opening bracket
  for (int i = 0; i<n_ic; i++) {
    if (i==0) {
      token=strtok(arg, space);
      token++;
    } else {
      token=strtok(NULL, space);
    }
    if (token==0) {
      return TCL_ERROR;
      Tcl_AppendResult(interp, "Unexpected argument ", token, (char *)NULL); 
    }
    // convert to float
    scan_succes = sscanf(token,"%lf",&(iccp3m_cfg.nvectorx[i]));
    if (!scan_succes) {
      Tcl_AppendResult(interp, "Unexpected argument ", token, (char *)NULL); 
      return TCL_ERROR;
    } 
    token=strtok(NULL, space);
    if (token==0) {
      Tcl_AppendResult(interp, "Unexpected argument ", token, (char *)NULL); 
      return TCL_ERROR;
    }
    // convert to float
    scan_succes = sscanf(token,"%lf",&(iccp3m_cfg.nvectory[i]));
    if (!scan_succes) {
      Tcl_AppendResult(interp, "Unexpected argument ", token, (char *)NULL); 
      return TCL_ERROR;
    } 
    token=strtok(NULL, closing_bracket);
    if (token==0) {
      Tcl_AppendResult(interp, "Unexpected argument ", token, (char *)NULL); 
      return TCL_ERROR;
    }
    // convert to float
    scan_succes = sscanf(token,"%lf",&(iccp3m_cfg.nvectorz[i]));
    if (!scan_succes) {
      Tcl_AppendResult(interp, "Unexpected argument ", token, (char *)NULL); 
      return TCL_ERROR;
    } 

    token=strtok(NULL, opening_bracket);

  }
  free(arg);
  return TCL_OK;
}
  

int tclcommand_iccp3m_parse_double_list(Tcl_Interp *interp, int n_ic, char *string, int flag) {
  /* This function parses a vector give as C-String */
  /* It currently is not too elegant. But works. */
  int size,i=0;
  int scan_succes;
  static double *numbers=NULL;
  const char delimiters[] = " ";
  char *token,*cp;
  float temp;

  size= n_ic;
  numbers = malloc((size)*sizeof(double));

  cp = strdup(string);                /* Make writable copy.  */
  token = strtok (cp, delimiters);
  scan_succes = sscanf(token,"%f",&temp);
  if (scan_succes < 1) 
    return 1;

  numbers[0]=temp;
    
  for(i=1;i<size;i++) {
    token = strtok (NULL, delimiters);
    if(token == NULL) {
      Tcl_AppendResult(interp, "Unexpected argument ", token, (char *)NULL); 
      return TCL_ERROR;
    }
    scan_succes = sscanf(token,"%lf",&(numbers[i]));
    if (scan_succes < 1 ) {
      Tcl_AppendResult(interp, "Unexpected argument ", token, (char *)NULL); 
      return TCL_ERROR;
    }
  }

  switch(flag) {
    case ICCP3M_AREA: 
      iccp3m_cfg.areas = (double*) realloc(iccp3m_cfg.areas, (size)*sizeof(double)); 
      for( i = 0 ; i < size ; i++ )  
        iccp3m_cfg.areas[i]=numbers[i];
      break;
    case ICCP3M_EPSILON:
      iccp3m_cfg.ein = (double*) realloc(iccp3m_cfg.ein,(size)*sizeof(double));
      for( i = 0 ; i < size; i++)  
        iccp3m_cfg.ein[i]=numbers[i];
      break;
    case ICCP3M_SIGMA:
      iccp3m_cfg.sigma = (double*) realloc(iccp3m_cfg.sigma,(size)*sizeof(double));
      for( i = 0 ; i < size; i++)  {
        iccp3m_cfg.sigma[i]=numbers[i];
      }

    break;
  }

  free(numbers);
  return (0);
}
#endif

