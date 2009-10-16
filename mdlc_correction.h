// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2009; all rights reserved unless otherwise stated.
/** \file p3m.h   main header-file for MDLC (magnetic dipolar layer correction).
 *
 *  Developer: Joan J. Cerda.
 *  Purpose:   get the corrections for dipolar 3D algorithms 
 *             when applied to a slab geometry and dipolar
 *	      particles. DLC & co
 *  Article:   A. Brodka, Chemical Physics Letters 400, 62-67 (2004).
 * 	      
 *	      We also include a tuning function that returns the
 *	      cut-off necessary to attend a certain accuracy.
 *	      
 *  Restrictions: the slab must be such that the z is the short 
 *                direction. Othewise we get trash.    	      
 * 
 *  Limitations:  at this moment it is restricted to work with 1 cpu
 */

#ifndef DLC_DIPOLAR_H
#define DLC_DIPOLAR_H


#ifdef MDLC
#ifdef MAGNETOSTATICS
#ifdef DIPOLES

/** parameters for the MDLC method */
typedef struct {
  /** maximal pairwise error of the potential and force */
  double maxPWerror;
  
  /** the cutoff of the exponential sum. Since in all other MMM methods this is
      the far formula, we call it here the same, although in the ELC context it
      does not make much sense. */
  double far_cut;
  
  /** size of the empty gap. Note that ELC relies on the user to make sure that
      this condition is fulfilled. */
   double gap_size;
   
  /** whether the cutoff was set by the user, or calculated by Espresso. In the latter case, the
      cutoff will be adapted if important parameters, such as the box dimensions, change. */
  int far_calculated;
  
  /** up to where particles can be found */
  double h;

} DLC_struct;
extern DLC_struct dlc_params;

   int         mdlc_sanity_checks(); 
   void       add_mdlc_force_corrections();
   double   add_mdlc_energy_corrections();
   int         inter_parse_mdlc_params(Tcl_Interp * interp, int argc, char ** argv) ; 
   int printMDLCToResult(Tcl_Interp *interp);
   
   
#endif  /*of DIPOLES */
#endif   
#endif  /* of ifdef MDLC */


#endif

