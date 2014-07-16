#include "ewaldgpu_tcl.hpp"
#include "forces.hpp"
#include "EwaldgpuForce.hpp"

#ifdef EWALD_GPU

int tclprint_to_result_ewaldgpu(Tcl_Interp *interp)
{
  char buffer[TCL_DOUBLE_SPACE];

  Tcl_PrintDouble(interp, ewaldgpu_params.rcut, buffer);
  Tcl_AppendResult(interp, "ewaldgpu ", buffer, " ", (char *) NULL);
  sprintf(buffer,"%d",ewaldgpu_params.num_kx);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  sprintf(buffer,"%d",ewaldgpu_params.num_ky);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  sprintf(buffer,"%d",ewaldgpu_params.num_kz);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, ewaldgpu_params.alpha, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, ewaldgpu_params.accuracy, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);

  return TCL_OK;
}
int tclcommand_inter_coulomb_parse_ewaldgpu(Tcl_Interp * interp, int argc, char ** argv)
{
	double r_cut;
	int num_kx;
	int num_ky;
	int num_kz;
	double accuracy;
	double alpha=-1;
	IntList il;
	init_intlist(&il);

  if (argc < 1)
  {
  	Tcl_AppendResult(interp, "expected: inter coulomb <bjerrum> ewaldgpu <r_cut> (<K_cut> | {<K_cut,x> <K_cut,y><K_cut,z>}) <alpha> \nexpected: inter coulomb <bjerrum> ewaldgpu tune <accuracy> [<precision>] \nexpected: inter coulomb <bjerrum> ewaldgpu tunealpha <r_cut> <K_cut> [<precision>]",(char *) NULL);
    return TCL_ERROR;
  }

  if (ARG0_IS_S("tune"))
  {
    int status = tclcommand_inter_coulomb_parse_ewaldgpu_tune(interp, argc-1, argv+1, 0);
    if(status==TCL_OK) return TCL_OK;
    if(status==TCL_ERROR)
    {
    	Tcl_AppendResult(interp, "Accuracy could not been reached. Choose higher K_max or lower accuracy",(char *) NULL);
    	    return TCL_ERROR;
    }

  }

  if (ARG0_IS_S("tunealpha"))
    return tclcommand_inter_coulomb_parse_ewaldgpu_tunealpha(interp, argc-1, argv+1);

  if(! ARG0_IS_D(r_cut))
    return TCL_ERROR;

  if(argc < 3 || argc > 5) {
  	Tcl_AppendResult(interp, "expected: inter coulomb <bjerrum> ewaldgpu <r_cut> (<K_cut> | {<K_cut,x> <K_cut,y><K_cut,z>}) <alpha> \nexpected: inter coulomb <bjerrum> ewaldgpu tune <accuracy> [<precision>] \nexpected: inter coulomb <bjerrum> ewaldgpu tunealpha <r_cut> <K_cut> [<precision>]",(char *) NULL);
    return TCL_ERROR;
  }

  if(! ARG_IS_I(1, num_kx))
  {
    if( ! ARG_IS_INTLIST(1, il) || !(il.n == 3) )
    {
      Tcl_AppendResult(interp, "integer or interger list of length 3 expected", (char *) NULL);
      return TCL_ERROR;
    }
    else
    {
      num_kx = il.e[0];
      num_ky = il.e[1];
      num_kz = il.e[2];
    }
  }
  else
  {
  	num_kz = num_ky = num_kx;
  }

  if(argc > 2)
  {
    if(! ARG_IS_D(2, alpha))
      return TCL_ERROR;
  }
  else
  {
    Tcl_AppendResult(interp, "Automatic ewaldgpu tuning not implemented.",
		     (char *) NULL);
    return TCL_ERROR;
  }

  if(argc > 3)
  {
    if(! ARG_IS_D(3, accuracy))
    {
      Tcl_AppendResult(interp, "accuracy double expected", (char *) NULL);
      return TCL_ERROR;
    }
  }

  // Create object
  EwaldgpuForce *A=new EwaldgpuForce(r_cut, num_kx, num_ky, num_kz, alpha);
  FI.addMethod(A);
  rebuild_verletlist = 1;
  ewaldgpu_params.ewaldgpu_is_running = true;
  ewaldgpu_params.isTuned = true;
	mpi_bcast_coulomb_params();
	mpi_bcast_event(INVALIDATE_SYSTEM);
  return TCL_OK;
}
int tclcommand_inter_coulomb_parse_ewaldgpu_tune(Tcl_Interp * interp, int argc, char ** argv, int adaptive)
{
	double r_cut;
	double alpha;
	int num_kx;
	int num_ky;
	int num_kz;
	int K_max = 30;
	int time_calc_steps = 0;
	double accuracy = 0.0001;
	double precision = 0.000001;

  while(argc > 0)
  {
    if(ARG0_IS_S("accuracy"))
    {
      if(! (argc > 1 && ARG1_IS_D(accuracy) && accuracy > 0))
      {
				Tcl_AppendResult(interp, "accuracy expects a positive double ",(char *) NULL);
				return TCL_ERROR;
      }
    }
    else if(ARG0_IS_S("precision"))
    {
      if(! (argc > 1 && ARG1_IS_D(precision) && precision > 0))
      {
				Tcl_AppendResult(interp, "precision expects a positive double ",(char *) NULL);
				return TCL_ERROR;
      }
    }
    else if(ARG0_IS_S("K_max"))
    {
      if(! (argc > 1 && ARG1_IS_I(K_max) && K_max > 0))
      {
				Tcl_AppendResult(interp, "K_max expects a positive integer ",(char *) NULL);
				return TCL_ERROR;
      }
    }
    else if(ARG0_IS_S("time_calc_steps"))
    {
      if(! (argc > 1 && ARG1_IS_I(time_calc_steps) && time_calc_steps > 0))
      {
				Tcl_AppendResult(interp, "time_calc_steps expects a positive integer ",(char *) NULL);
				return TCL_ERROR;
      }
    }
    /* unknown parameter. Probably one of the optionals */
    else break;

    argc -= 2;
    argv += 2;
  }

  ewaldgpu_set_params_tune(accuracy, precision, K_max, time_calc_steps);

  /* Create object */
  EwaldgpuForce *A=new EwaldgpuForce(r_cut, num_kx, num_ky, num_kz, alpha);
  FI.addMethod(A);
  rebuild_verletlist = 1;

  /* do the tuning */
  char *log = NULL;
  if (ewaldgpu_adaptive_tune(&log) == ES_ERROR)
  {
    Tcl_AppendResult(interp, log, "\nfailed to tune ewaldgpu parameters to required accuracy ", (char *) NULL);
    if (log) free(log);
    return TCL_ERROR;
  }

  /* Tell the user about the tuning outcome */
  Tcl_AppendResult(interp, log, (char *) NULL);
  if (log) free(log);

  rebuild_verletlist = 1;
	mpi_bcast_coulomb_params();
	mpi_bcast_event(INVALIDATE_SYSTEM);
  return TCL_OK;
}
int tclcommand_inter_coulomb_parse_ewaldgpu_tunealpha(Tcl_Interp * interp, int argc, char ** argv)
{
	double r_cut;
	double alpha;
	int num_kx;
	int num_ky;
	int num_kz;
	double precision=0.000001;

	IntList il;
	init_intlist(&il);

  if (argc < 3)
  {
    Tcl_AppendResult(interp, "wrong # arguments: <r_cut> <K_cut,x> <K_cut,y><K_cut,z> [<precision>]  ", (char *) NULL);
    return TCL_ERROR;
  }

  /* PARSE EWALD COMMAND LINE */
  /* epsilon */
  if (! ARG_IS_D(0, r_cut))
  {
    Tcl_AppendResult(interp, "<r_cut> should be a double",(char *) NULL);
  }
  /* k_cut */
  if(! ARG_IS_I(1, num_kx))
  {
    if( ! ARG_IS_INTLIST(1, il) || !(il.n == 3) )
    {
      Tcl_AppendResult(interp, "integer or integer list of length 3 expected", (char *) NULL);
      return TCL_ERROR;
    }
    else
    {
      num_kx = il.e[0];
      num_ky = il.e[1];
      num_kz = il.e[2];
    }
  }
  else
  {
  	num_kz = num_ky = num_kx;
  }
  /* precision */
  if (! ARG_IS_D(2, precision))
  {
    Tcl_AppendResult(interp, "<precision> should be a double",
		     (char *) NULL);
    return 0;
  }
  //Compute alpha
  double q_sqr = ewaldgpu_compute_q_sqare();
  alpha = ewaldgpu_compute_optimal_alpha(r_cut, num_kx, num_ky, num_kz, q_sqr, box_l, precision);
  ewaldgpu_params.isTuned = true;
  mpi_bcast_coulomb_params();
  mpi_bcast_event(INVALIDATE_SYSTEM);
  // Create object
  EwaldgpuForce *A=new EwaldgpuForce(r_cut, num_kx, num_ky, num_kz, alpha);
  FI.addMethod(A);
  rebuild_verletlist = 1;

  return TCL_OK;
}


#endif
