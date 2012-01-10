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
//

/** \file iccp3m.h 

    ICCP3M is a method that allows to take into account the influence
    of arbitrarliy shaped dielectric interfaces.  The dielectric
    properties of a dielectric medium in the bulk of the simulation
    box are taken into account by reproducing the jump in the electric
    field at the inface with charge surface segments. The charge
    density of the surface segments have to be determined
    self-consistently using an iterative scheme.  It can at presently
    - despite its name - be used with P3M, ELCP3M, MMM2D and MMM1D.
    For details see:<br> S. Tyagi, M. Suzen, M. Sega, C. Holm,
    M. Barbosa: A linear-scaling method for computing induced charges
    on arbitrary dielectric boundaries in large system simulations
    (Preprint)

    To set up ICCP3M first the dielectric boundary has to be modelled
    by espresso particles 0..n where n has to be passed as a parameter
    to ICCP3M. This is still a bit inconvenient, as it forces the user
    to reserve the first n particle ids to wall charges, but as the
    other parts of espresso do not suffer from a limitation like this,
    it can be tolerated.
    
    For the determination of the induced charges only the forces
    acting on the induced charges has to be determined. As P3M an the
    other coulomb solvers calculate all mutual forces, the force
    calculation was modified to avoid the calculation of the short
    range part of the source-source force calculation.  For different
    particle data organisation schemes this is performed differently.
    */

#ifndef _ICCP3M_TCL_H
#define _ICCP3M_TCL_H

#if defined(ELECTROSTATICS)

/** Implementation of the tcl-command <br>
    iccp3m  { \<last_ind_id\> \<e1\> \<num_iteration\> \<convergence\> \<relaxation\> \<area\> \<normal_components\> \<e_in/e_out\>  [\<ext_field\>] | iterate } 
    ICC sets up and calculates induced charges on dielectric surfaces. At the beginning of every simulation run particles on the surface boundary 
    have to be set up (before any real particle) together with the list of areas, normal vectors and dielectric constant associated with them. 
    After that the iterate flag can be used during the simulation to update the value of the induced charges.
    
    Parameters: <br>
                 \<last_ind_id\> ID of the last surface charge. Note that the IDs of the surface charges must range from 0 to \<last_ind_id\>
                 \<e1\>          = Dielectric Constant of the Bulk accessible to free particles 
                 \<num_iteration\> = Maximum number of ICCP3M iterations calculating the induced charges. 
                 \<relaxation\> = Relaxaxion parameter \f$omega\f$ for the successive over-relaxation scheme. 
                 \<area\>       = List of the areas of each surface element.
                 \<normal_components\> = List of normal vectors of each surface element. 3n list entries. Do not have to be normalized.
                 \<e_in/e_out\> = Ratio of dielectric co


                 iterate         = Indicates that a previous surface discretization shall be used. T
*/
int tclcommand_iccp3m(ClientData data, Tcl_Interp *interp, int argc, char **argv);

#endif /* ELECTROSTATICS */

#endif /* ICCP3M_H */
