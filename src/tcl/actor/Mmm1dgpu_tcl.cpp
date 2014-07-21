#include "Mmm1dgpu_tcl.hpp"
#include "forces.hpp"
#include "actor/Mmm1dgpuForce.hpp"
#include "mmm1d.hpp"
#include "EspressoSystemInterface.hpp"

#ifdef MMM1D_GPU

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

  addMmm1dgpuForce (maxPWerror, switch_rad, bessel_cutoff);	// if this command is run multiple times, subsequent runs will only update parameters

  return TCL_OK;
}

#endif
