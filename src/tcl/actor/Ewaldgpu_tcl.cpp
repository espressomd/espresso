#include "Ewaldgpu_tcl.hpp"
#include "forces.hpp"
#include "energy.hpp"
#include "actor/EwaldGPU.hpp"
#include "EspressoSystemInterface.hpp"

#ifdef EWALD_GPU

EwaldgpuForce *ewaldgpuForce;

int tclprint_to_result_ewaldgpu(Tcl_Interp *interp)
{
	//Used in tcl-script output
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
  //PARSE EWALD COMMAND LINE
  if (argc < 1)
  {
  	Tcl_AppendResult(interp, "\nExpected: inter coulomb <bjerrum> ewaldgpu <r_cut> (<K_cut> | {<K_cut,x> <K_cut,y><K_cut,z>}) <alpha> \nExpected: inter coulomb <bjerrum> ewaldgpu tune <accuracy> [<precision>] \nExpected: inter coulomb <bjerrum> ewaldgpu tunealpha <r_cut> <K_cut> [<precision>]",(char *) NULL);
    return TCL_ERROR;
  }

  //Tune
  if (ARG0_IS_S("tune"))
  {
    return tclcommand_inter_coulomb_parse_ewaldgpu_tune(interp, argc-1, argv+1, 0);
  }
  //Tune alpha
  else if (ARG0_IS_S("tunealpha"))
  {
  	return tclcommand_inter_coulomb_parse_ewaldgpu_tunealpha(interp, argc-1, argv+1);
  }
  //No tuning
  else
  {
  	return tclcommand_inter_coulomb_parse_ewaldgpu_notune(interp, argc, argv);
  }
}

int tclcommand_inter_coulomb_parse_ewaldgpu_notune(Tcl_Interp * interp, int argc, char ** argv)
{
	double r_cut=-1;
	int num_kx=-1;
	int num_ky=-1;
	int num_kz=-1;
	double alpha=-1;
	IntList il;
	init_intlist(&il);

	if(argc < 3 || argc > 5) {
		Tcl_AppendResult(interp, "\nExpected: inter coulomb <bjerrum> ewaldgpu <r_cut> (<K_cut> | {<K_cut,x> <K_cut,y><K_cut,z>}) <alpha>",(char *) NULL);
		return TCL_ERROR;
	}
	if(! ARG_IS_D(0,r_cut))
	{
		Tcl_AppendResult(interp, "\n<r_cut> double expected", (char *) NULL);
		return TCL_ERROR;
	}
	if(! ARG_IS_I(1, num_kx))
	{
		if( ! ARG_IS_INTLIST(1, il) || !(il.n == 3) )
		{
			Tcl_AppendResult(interp, "\n(<K_cut> | {<K_cut,x> <K_cut,y><K_cut,z>}) integer or integer list of length 3 expected", (char *) NULL);
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
		{
			Tcl_AppendResult(interp, "\n<alpha> double expected", (char *) NULL);
			return TCL_ERROR;
		}
	}
	else
	{
		Tcl_AppendResult(interp, "\nAutomatic ewaldgpu tuning not implemented.",
				 (char *) NULL);
		return TCL_ERROR;
	}

	//Turn on ewaldgpu
	if (!ewaldgpuForce) // inter coulomb ewaldgpu was never called before
	{
		ewaldgpuForce = new EwaldgpuForce(espressoSystemInterface, r_cut, num_kx, num_ky, num_kz, alpha);
		forceActors.add(ewaldgpuForce);
		energyActors.add(ewaldgpuForce);
	}
	//Broadcast parameters
	ewaldgpuForce->set_params(r_cut, num_kx, num_ky, num_kz, alpha);
	coulomb.method = COULOMB_EWALD_GPU;
	ewaldgpu_params.isTunedFlag = false;
	ewaldgpu_params.isTuned = true;
	rebuild_verletlist = 1;
	mpi_bcast_coulomb_params();
	return TCL_OK;
}
int tclcommand_inter_coulomb_parse_ewaldgpu_tune(Tcl_Interp * interp, int argc, char ** argv, int adaptive)
{
	double r_cut=-1;
	double alpha=-1;
	int num_kx=-1;
	int num_ky=-1;
	int num_kz=-1;
	int K_max = 30;
	int time_calc_steps = 100;
	double accuracy = 0.0001;
	double precision = 0.000001;

	//PARSE EWALD COMMAND LINE
	while(argc > 0)
	{
		if(ARG0_IS_S("accuracy"))
		{
			if(! (argc > 1 && ARG1_IS_D(accuracy) && accuracy > 0))
			{
				Tcl_AppendResult(interp, "\n<accuracy> expects a positive double ",(char *) NULL);
				return TCL_ERROR;
			}
		}
		else if(ARG0_IS_S("precision"))
		{
			if(! (argc > 1 && ARG1_IS_D(precision) && precision > 0))
			{
				Tcl_AppendResult(interp, "\n<precision> expects a positive double ",(char *) NULL);
				return TCL_ERROR;
			}
		}
		else if(ARG0_IS_S("K_max"))
		{
			if(! (argc > 1 && ARG1_IS_I(K_max) && K_max > 0))
			{
				Tcl_AppendResult(interp, "\n<K_max> expects a positive integer ",(char *) NULL);
				return TCL_ERROR;
			}
		}
		else if(ARG0_IS_S("time_calc_steps"))
		{
			if(! (argc > 1 && ARG1_IS_I(time_calc_steps) && time_calc_steps > 0))
			{
				Tcl_AppendResult(interp, "\n<time_calc_steps> expects a positive integer ",(char *) NULL);
				return TCL_ERROR;
			}
		}
		else break;

		argc -= 2;
		argv += 2;
	}

  //Turn on ewaldgpu
  ewaldgpuForce->set_params_tune(accuracy, precision, K_max, time_calc_steps);
  if (!ewaldgpuForce) // inter coulomb ewaldgpu was never called before
  {
	  ewaldgpuForce = new EwaldgpuForce(espressoSystemInterface, r_cut, num_kx, num_ky, num_kz, alpha);
	  forceActors.add(ewaldgpuForce);
	  energyActors.add(ewaldgpuForce);
  }

  //Broadcast parameters
  coulomb.method = COULOMB_EWALD_GPU;
  ewaldgpu_params.isTunedFlag = false;
  rebuild_verletlist = 1;
  mpi_bcast_coulomb_params();
  //Tuning
  char *log = NULL;
  if (ewaldgpuForce->adaptive_tune(&log,espressoSystemInterface) == ES_ERROR)
  {
    Tcl_AppendResult(interp,  "\nAccuracy could not been reached. Choose higher K_max or lower accuracy", (char *) NULL);
    return TCL_ERROR;
  }
  //Tell the user about the tuning outcome
  Tcl_AppendResult(interp, log, (char *) NULL);
  if (log) free(log);

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

  //PARSE EWALD COMMAND LINE
  if (argc < 3)
  {
    Tcl_AppendResult(interp, "\nWrong # arguments: <r_cut> (<K_cut> | {<K_cut,x> <K_cut,y><K_cut,z>}) <precision>", (char *) NULL);
    return TCL_ERROR;
  }
  if (! ARG0_IS_D(r_cut))
  {
    Tcl_AppendResult(interp, "\n<r_cut> should be a double",(char *) NULL);
    return TCL_ERROR;
  }
  if(! ARG_IS_I(1, num_kx))
  {
    if( ! ARG_IS_INTLIST(1, il) || !(il.n == 3) )
    {
      Tcl_AppendResult(interp, "\n(<K_cut> | {<K_cut,x> <K_cut,y><K_cut,z>}) integer or integer list of length 3 expected", (char *) NULL);
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
  if (! ARG_IS_D(2, precision))
  {
    Tcl_AppendResult(interp, "\n<precision> should be a double", (char *) NULL);
    return TCL_ERROR;
  }

  //Compute alpha
	Particle *particle;
	particle = (Particle*)Utils::malloc(n_part*sizeof(Particle));
	mpi_get_particles(particle, NULL);
	double q_sqr = ewaldgpuForce->compute_q_sqare(particle);
  alpha = ewaldgpuForce->compute_optimal_alpha(r_cut, num_kx, num_ky, num_kz, q_sqr, box_l, precision);

  //Turn on ewaldgpu
  if (!ewaldgpuForce) // inter coulomb ewaldgpu was never called before
  {
	  ewaldgpuForce = new EwaldgpuForce(espressoSystemInterface, r_cut, num_kx, num_ky, num_kz, alpha);
	  forceActors.add(ewaldgpuForce);
	  energyActors.add(ewaldgpuForce);
  }
  //Broadcast parameters
  coulomb.method = COULOMB_EWALD_GPU;
  ewaldgpu_params.isTunedFlag = false;
  ewaldgpu_params.isTuned = true;
  rebuild_verletlist = 1;
  mpi_bcast_coulomb_params();
  return TCL_OK;
}

#endif

