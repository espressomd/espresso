/*
  Copyright (C) 2010 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
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

/** \file magnetic_non_p3m__methods.h   Header of all 3d non P3M methods to deal with the magnetic dipoles
 *   
 *  DAWAANR => DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA
 *   Handling of a system of dipoles where no replicas exist (i.e. it is not replicated) and we DO NOT 
 *    assume minimum image convention, i.e.,  we work with the unfolded coordinates 
 *
 *   Specific application for which it was intendeed:  handling of a single magnetic filament, where all particles are bonded and hence, it 
 *   it has sense to work with unfolded coordinate, to avoid minimum image convention and replicating the chain when computing the
 *   interactions
 *
 * 
 *  DS => Direct sum , compute the interactions via direct sum, 
 *
 *  For more information about the DAWAANR => DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA
 *  see \ref magnetic_non_p3m__methods.c "magnetic_non_p3m__methods.c"
 */




/** \file magnetic_non_p3m__methods.h  header for the method for  handling  a system of dipoles where no replicas exist 
 *                            (i.e. it is not replicated) and we DO NOT assume minimum image convention,
 *                            i.e.,  we work with the unfolded coordinates 
 *
 *                                    Specific application for which it was intendeed:  handling of a single magnetic filament, where all particles are bonded and hence, it 
 *                                    it has sense to work with unfolded coordinate, to avoid minimum image convention and replicating the chain when computing the
 *                                    interactions
 *
 *  For more information about the DAWAANR => DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA
 *  see \ref magnetic_non_p3m__methods.c "magnetic_non_p3m__methods.c"
 */


#ifndef  MAG_NON_P3M_H
#define MAG_NON_P3M_H



#ifdef MAGNETOSTATICS




/* =============================================================================
                  DAWAANR => DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA                
   =============================================================================
*/
#ifdef DAWAANR
     /*  Information about the status of the method */
     int tclprint_to_result_DAWAANR(Tcl_Interp *interp);
         
  /* Parsing function for the dawaanr method*/
  int tclcommand_inter_magnetic_parse_dawaanr(Tcl_Interp * interp, int argc, char ** argv);

   /* Sanity checks for the dawaanr*/
  int DAWAANR_sanity_checks();

  /* Core of the method: here you compute all the magnetic forces, torques and the magnetic energy for the whole system*/
   double dawaanr_calculations(int force_flag, int energy_flag) ;
#endif/*of  DAWAANR */



/* =============================================================================
                  DIRECT SUM FOR MAGNETIC SYSTEMS               
   =============================================================================
*/
#ifdef MAGNETIC_DIPOLAR_DIRECT_SUM

    /*  Information about the status of the method */
     int tclprint_to_result_Magnetic_dipolar_direct_sum_(Tcl_Interp *interp);


  /* Parsing function for the magnetic dipolar direct sum method*/
  int tclcommand_inter_magnetic_parse_mdds(Tcl_Interp * interp, int argc, char ** argv);

  /* Sanity checks for the magnetic dipolar direct sum*/
  int magnetic_dipolar_direct_sum_sanity_checks();

 /* Core of the method: here you compute all the magnetic forces,torques and the energy for the whole system using direct sum*/
double  magnetic_dipolar_direct_sum_calculations(int force_flag, int energy_flag);

#endif /*of ifdef MAGNETIC_DIPOLAR_DIRECT_SUM */






#endif /*of ifdef MAGNETOSTATICS  */
#endif /* of ifndef  MAG_NON_P3M_H */
