/*
  Copyright (C) 2014,2015,2016 The ESPResSo project
  
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
#include "Mmm1dgpu_tcl.hpp"
#include "forces.hpp"
#include "energy.hpp"
#include "actor/Mmm1dgpuForce.hpp"
#include "mmm1d.hpp"
#include "EspressoSystemInterface.hpp"

#ifdef MMM1D_GPU

Mmm1dgpuForce *mmm1dgpuForce;

int tclprint_to_result_MMM1DGPU (Tcl_Interp * interp)
{
  char buffer[TCL_DOUBLE_SPACE];

  Tcl_PrintDouble (interp, sqrt (mmm1d_params.far_switch_radius_2), buffer);
  Tcl_AppendResult (interp, "mmm1dgpu ", buffer, " ", (char *) NULL);
  sprintf (buffer, "%d", mmm1d_params.bessel_cutoff);
  Tcl_AppendResult (interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble (interp, mmm1d_params.maxPWerror, buffer);
  Tcl_AppendResult (interp, buffer, (char *) NULL);

  return TCL_OK;
}

//void Mmm1dgpuForce::disable()
//{
//	if (!mmm1dgpuForce) {
//		forceActors.remove(mmm1dgpuForce);
//		if (coulomb.method == COULOMB_MMM1D_GPU) {
//			coulomb.method = COULOMB_NONE;
//			mpi_bcast_coulomb_params();
//		}
//	}
//}

int tclcommand_inter_coulomb_parse_mmm1dgpu (Tcl_Interp * interp, int argc, char **argv)
{
  double switch_rad, maxPWerror;
  int bessel_cutoff;

  if (argc < 2)
  {
    Tcl_AppendResult (interp, "wrong # arguments: inter coulomb mmm1dgpu <switch radius> " 
      "{<bessel cutoff>} <maximal error for near formula> | tune  <maximal pairwise error>", (char *) NULL);
    return TCL_ERROR;
  }

  if (ARG0_IS_S ("tune"))
  {
    /* autodetermine bessel cutoff AND switching radius */
    if (!ARG_IS_D (1, maxPWerror))
      return TCL_ERROR;
    bessel_cutoff = -1;
    switch_rad = -1;
  }
  else
  {
    if (argc == 2)
  	{
  	  /* autodetermine bessel cutoff */
  	  if ((!ARG_IS_D (0, switch_rad)) || (!ARG_IS_D (1, maxPWerror)))
  	    return TCL_ERROR;
  	  bessel_cutoff = -1;
  	}
    else if (argc == 3)
    {
      /* fully manual */
      if ((!ARG_IS_D (0, switch_rad)) || (!ARG_IS_I (1, bessel_cutoff)) || (!ARG_IS_D (2, maxPWerror)))
        return TCL_ERROR;

      if (bessel_cutoff <= 0)
      {
        Tcl_AppendResult (interp, "bessel cutoff too small", (char *) NULL);
        return TCL_ERROR;
      }
	  }
    else
  	{
  	  Tcl_AppendResult (interp, "wrong # arguments: inter coulomb mmm1dgpu <switch radius> " "{<bessel cutoff>} <maximal error for near formula> | tune  <maximal pairwise error>", (char *) NULL);
  	  return TCL_ERROR;
  	}

    if (switch_rad <= 0 || switch_rad > box_l[2])
  	{
  	  Tcl_AppendResult (interp, "switching radius is not between 0 and box_l[2]", (char *) NULL);
  	  return TCL_ERROR;
  	}
  }

  // turn on MMM1DGPU

  /* coulomb.prefactor apparently is not available yet at this point */
  if (!mmm1dgpuForce) // inter coulomb mmm1dgpu was never called before
  {
	  mmm1dgpuForce = new Mmm1dgpuForce(espressoSystemInterface, 0, maxPWerror,
			  switch_rad, bessel_cutoff);
	  forceActors.add(mmm1dgpuForce);
	  energyActors.add(mmm1dgpuForce);
  }
  /* set_params needs to be called both upon create and upon update because it is responsible for writing
		to the struct from which the TCL command "inter coulomb" retrieves the current parameter set */
  mmm1dgpuForce->set_params(0, 0, maxPWerror, switch_rad, bessel_cutoff, true);

  coulomb.method = COULOMB_MMM1D_GPU;
  mpi_bcast_coulomb_params();

  return TCL_OK;
}

#endif
