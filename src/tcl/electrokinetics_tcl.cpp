#include "config.hpp"
#include "electrokinetics_tcl.hpp"
#include "electrokinetics.hpp"
#include "lb-boundaries.hpp"
#include "lb-boundaries_tcl.hpp"
#include "initialize.hpp"

int tclcommand_electrokinetics(ClientData data, Tcl_Interp *interp, int argc, char **argv) {
#ifndef ELECTROKINETICS
  Tcl_AppendResult(interp, "Feature ELECTROKINETICS required", NULL);
  return TCL_ERROR;
#else

  argc--;
  argv++;

  int err = TCL_OK;
  int species;
  double floatarg;
  double vectarg[3];
  int coord[3];
  char int_buffer[TCL_INTEGER_SPACE];
  char double_buffer[TCL_DOUBLE_SPACE];

#ifdef EK_REACTION
  int reactant,
      product0,
      product1;
  double ct_rate,
         rho_reactant_reservoir, 
         rho_product0_reservoir, 
         rho_product1_reservoir, 
         radius,
         fraction0,
         fraction1;
#endif

  if(argc < 2) {
    Tcl_AppendResult(interp, "Usage of \"electrokinetics\":", (char *)NULL);
    Tcl_AppendResult(interp, "electrokinetics [agrid #float] [viscosity #float] [friction #float]\n", (char *)NULL);
    Tcl_AppendResult(interp, "                [bulk_viscosity #float] [gamma_even #float] [gamma_odd #float]\n", (char *)NULL);
    Tcl_AppendResult(interp, "electrokinetics print <density|velocity|potential|pressure|lbforce> vtk #string]\n", (char *)NULL);
    Tcl_AppendResult(interp, "electrokinetics node i j k print velocity\n", (char *)NULL);
    Tcl_AppendResult(interp, "electrokinetics reaction [reactant_index #int]\n", (char *)NULL);
    Tcl_AppendResult(interp, "                         [product0_index #int]\n", (char *)NULL);
    Tcl_AppendResult(interp, "                         [product1_index #int]\n", (char *)NULL);
    Tcl_AppendResult(interp, "                         [reactant_resrv_density #float]\n", (char *)NULL);
    Tcl_AppendResult(interp, "                         [product0_resrv_density #float]\n", (char *)NULL);
    Tcl_AppendResult(interp, "                         [product1_resrv_density #float]\n", (char *)NULL);
    Tcl_AppendResult(interp, "                         [reaction_rate #float]\n", (char *)NULL);
    Tcl_AppendResult(interp, "                         [reaction_radius #float]\n", (char *)NULL);
    Tcl_AppendResult(interp, "                         [reaction_fraction_pr_0 #float]\n", (char *)NULL);
    Tcl_AppendResult(interp, "                         [reaction_fraction_pr_1 #float]\n", (char *)NULL);
    Tcl_AppendResult(interp, "electrokinetics accelerated_frame on\n", (char *)NULL);
    Tcl_AppendResult(interp, "                  boundary_mass_density #double\n", (char *)NULL);
    Tcl_AppendResult(interp, "                  ext_acceleration_force #double #double #double\n", (char *)NULL);
    Tcl_AppendResult(interp, "electrokinetics accelerated_frame off\n", (char *)NULL);
    Tcl_AppendResult(interp, "electrokinetics accelerated_frame print boundary_velocity\n", (char *)NULL);
    Tcl_AppendResult(interp, "electrokinetics boundary charge_density #float [wall ]\n", (char *)NULL); //TODO full description
    Tcl_AppendResult(interp, "                                               [sphere ]\n", (char *)NULL);
    Tcl_AppendResult(interp, "                                               [cylinder ]\n", (char *)NULL);
    Tcl_AppendResult(interp, "                                               [rhomboid ]\n", (char *)NULL);
    Tcl_AppendResult(interp, "                                               [pore ]\n", (char *)NULL);
    Tcl_AppendResult(interp, "                                               [stomatocyte ]\n", (char *)NULL);
    Tcl_AppendResult(interp, "electrokinetics #int [density #float] [D #float] [T #float] [valency #float]\n", (char *)NULL);
    Tcl_AppendResult(interp, "                     [ext_force #float #float #float]\n", (char *)NULL);
    Tcl_AppendResult(interp, "                     [print density vtk #string]\n", (char *)NULL);
    Tcl_AppendResult(interp, "                     [print flux vtk #string]\n", (char *)NULL);
    return TCL_ERROR;
  }
  else if(ARG0_IS_S("boundary")) {
#ifndef EK_BOUNDARIES
    Tcl_AppendResult(interp, "Feature EK_BOUNDARIES required", (char *) NULL);
    return (TCL_ERROR);
#else
    argc--;
    argv++;
    
    if(!ARG0_IS_S("charge_density") || !ARG1_IS_D(floatarg)) {
      Tcl_AppendResult(interp, "You need to specify the boundary charge density using\n", (char *) NULL);
      Tcl_AppendResult(interp, "electrokinetics boundary charge_density #float ...\n", (char *)NULL);
      return (TCL_ERROR);
    }
    
    argc -= 2;
    argv += 2;
    
    if(floatarg != 0.0) {
      species = -1;
      
      for(unsigned int i = 0; i < ek_parameters.number_of_species; i++)
        if(ek_parameters.valency[i] != 0.0) {
          species = i;
          break;
        }
      
      if(species == -1) {
        Tcl_AppendResult(interp, "You need to define at least one charged species in order to use charged walls\n", (char *) NULL);
        return (TCL_ERROR);
      }
    }
    
    int c_num;
    LB_Boundary *lbboundary_tmp;
    
    if(ARG0_IS_S("wall")) {
      lbboundary_tmp = generate_lbboundary();
      err = tclcommand_lbboundary_wall(lbboundary_tmp, interp, argc - 1, argv + 1);
      lbboundary_tmp->charge_density = floatarg;
    }
    else if(ARG0_IS_S("sphere")) {
      lbboundary_tmp = generate_lbboundary();
      err = tclcommand_lbboundary_sphere(lbboundary_tmp, interp, argc - 1, argv + 1);
      lbboundary_tmp->charge_density = floatarg;
    }
    else if(ARG0_IS_S("cylinder")) {
      lbboundary_tmp = generate_lbboundary();
      err = tclcommand_lbboundary_cylinder(lbboundary_tmp, interp, argc - 1, argv + 1);
      lbboundary_tmp->charge_density = floatarg;
    }
    else if(ARG0_IS_S("rhomboid")) {
      lbboundary_tmp = generate_lbboundary();
      err = tclcommand_lbboundary_rhomboid(lbboundary_tmp, interp, argc - 1, argv + 1);
      lbboundary_tmp->charge_density = floatarg;
    }
    else if(ARG0_IS_S("pore")) {
      lbboundary_tmp = generate_lbboundary();
      err = tclcommand_lbboundary_pore(lbboundary_tmp, interp, argc - 1, argv + 1);
      lbboundary_tmp->charge_density = floatarg;
    }
    else if(ARG0_IS_S("stomatocyte")) {
      lbboundary_tmp = generate_lbboundary();
      err = tclcommand_lbboundary_stomatocyte(lbboundary_tmp, interp, argc - 1, argv + 1);
      lbboundary_tmp->charge_density = floatarg;
    }
    else if(ARG0_IS_S("delete")) {
      if(argc < 3) {
        /* delete all */
        Tcl_AppendResult(interp, "Can only delete individual electrokinetics boundaries", (char *) NULL);
        err = TCL_ERROR;
      }
      else {
        if(Tcl_GetInt(interp, argv[2], &(c_num)) == TCL_ERROR)
          return (TCL_ERROR);
          
        if(c_num < 0 || c_num >= n_lb_boundaries) {
	        Tcl_AppendResult(interp, "Can not delete non existing electrokinetics boundary", (char *) NULL);
	        return (TCL_ERROR);
        }
      }
    }
    else {
      Tcl_AppendResult(interp, "possible electrokinetics boundary charge_density #float parameters: wall, sphere, cylinder, rhomboid, pore, stomatocyte, delete {c} to delete a boundary", (char *) NULL);
      return (TCL_ERROR);
    }
        
    on_lbboundary_change();
#endif /* EK_BOUNDARIES */
  }
  else if(ARG0_IS_I(species)) {
    argc--;
    argv++;
    
    if(species < 0 || species > MAX_NUMBER_OF_SPECIES) {
      sprintf(int_buffer, "%d", MAX_NUMBER_OF_SPECIES);
      Tcl_AppendResult(interp, "electrokinetics #int requires a number between 0 and", int_buffer, "denoting the species\n", (char *)NULL);
      return TCL_ERROR;
    }
    
    while(argc > 0) {
      if(ARG0_IS_S("D")) {
        if(argc < 2 || !ARG1_IS_D(floatarg)) {
          Tcl_AppendResult(interp, "electrokinetics #int D requires one floating point number as argument\n", (char *)NULL);
          return TCL_ERROR;
        }
        else if(floatarg < 0) {
          Tcl_AppendResult(interp, "electrokinetics #int D can not be negative\n", (char *)NULL);
          return TCL_ERROR;
        }
        else {
          if(ek_set_D(species, floatarg) == 0) {
            argc -= 2;
            argv += 2;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error setting electrokinetics #int D\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if(ARG0_IS_S("density")) {
        if(argc < 2 || !ARG1_IS_D(floatarg)) {
          Tcl_AppendResult(interp, "electrokinetics #int density requires one floating point number as argument\n", (char *)NULL);
          return TCL_ERROR;
        }
        else if(floatarg < 0) {
          Tcl_AppendResult(interp, "electrokinetics #int density must be positive\n", (char *)NULL);
          return TCL_ERROR;
        }
        else {
          if(ek_set_density(species, floatarg) == 0) {
            argc -= 2;
            argv += 2;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error setting electrokinetics #int density\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if(ARG0_IS_S("ext_force")) {
        if(argc < 4 || !ARG_IS_D(1, vectarg[0]) || !ARG_IS_D(2, vectarg[1]) || !ARG_IS_D(3, vectarg[2])) {
          Tcl_AppendResult(interp, "electrokinetics #int ext_force requires three floating point numbers as arguments\n", (char *)NULL);
          return TCL_ERROR;
        }
        else if(ek_set_ext_force(species, vectarg[0], vectarg[1], vectarg[2]) == 0) {
          argc -= 4;
          argv += 4;
        }
        else {
          Tcl_AppendResult(interp, "Unknown error setting electrokinetics #int ext_force\n", (char *)NULL);
          return TCL_ERROR;
        }
      }
      else if(ARG0_IS_S("valency")) {
        if(argc < 2 || !ARG1_IS_D(floatarg)) {
          Tcl_AppendResult(interp, "electrokinetics #int valency requires one floating point number as argument\n", (char *)NULL);
          return TCL_ERROR;
        }
        else {
          if(ek_set_valency(species, floatarg) == 0) {
            argc -= 2;
            argv += 2;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error setting electrokinetics #int valency\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if(ARG0_IS_S("print")) {
        argc--;
        argv++;
        
        if(argc != 3 || !ARG1_IS_S("vtk") || ( !ARG0_IS_S("density") && !ARG0_IS_S("flux") ) ) {
          Tcl_AppendResult(interp, "Wrong usage of electrokinetics #int print <density|flux> vtk #string\n", (char *)NULL);
          return TCL_ERROR;
        }
        
        if(ARG0_IS_S("density")) {
          if(ek_print_vtk_density(species, argv[2]) == 0) {
            argc -= 3;
            argv += 3;

            if((err = gather_runtime_errors(interp, err)) != TCL_OK)
              return TCL_ERROR;
            else
              return TCL_OK;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error in electrokinetics #int print density vtk #string\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
        else if(ARG0_IS_S("flux")) {
          if(ek_print_vtk_flux(species, argv[2]) == 0) {
            argc -= 3;
            argv += 3;

            if((err = gather_runtime_errors(interp, err)) != TCL_OK)
              return TCL_ERROR;
            else
              return TCL_OK;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error in electrokinetics #int print flux vtk #string\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else {
    	  Tcl_AppendResult(interp, "unknown feature \"", argv[0],"\" of electrokinetics node i j k print\n", (char *)NULL);
    	  return TCL_ERROR ;
      }
    }
  }
  else {
    Tcl_ResetResult(interp);
    
    while(argc > 0) {
      if(ARG0_IS_S("print")) {
        argc--;
        argv++;
        
        if(argc != 3 || !ARG1_IS_S("vtk") || (!ARG0_IS_S("velocity") && !ARG0_IS_S("density") && !ARG0_IS_S("boundary") && !ARG0_IS_S("potential") && !ARG0_IS_S("pressure") && !ARG0_IS_S("lbforce"))) {
          Tcl_AppendResult(interp, "Wrong usage of electrokinetics print <velocity|density|potential|pressure> vtk #string\n", (char *)NULL);
          return TCL_ERROR;
        }
        
        if(ARG0_IS_S("velocity")) {
          if(ek_lb_print_vtk_velocity(argv[2]) == 0) {
            argc -= 3;
            argv += 3;

            if((err = gather_runtime_errors(interp, err)) != TCL_OK)
              return TCL_ERROR;
            else
              return TCL_OK;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error in electrokinetics print velocity vtk #string\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
        else if(ARG0_IS_S("density")) {
          if(ek_lb_print_vtk_density(argv[2]) == 0) {
            argc -= 3;
            argv += 3;

            if((err = gather_runtime_errors(interp, err)) != TCL_OK)
              return TCL_ERROR;
            else
              return TCL_OK;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error in electrokinetics print density vtk #string\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
        else if(ARG0_IS_S("boundary")) {
#ifndef EK_BOUNDARIES
          Tcl_AppendResult(interp, "Feature EK_BOUNDARIES required", (char *) NULL);
          return (TCL_ERROR);
#else
          if(lb_lbfluid_print_vtk_boundary(argv[2]) == 0) {
            argc -= 3;
            argv += 3;

            if((err = gather_runtime_errors(interp, err)) != TCL_OK)
              return TCL_ERROR;
            else
              return TCL_OK;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error in electrokinetics print boundary vtk #string\n", (char *)NULL);
            return TCL_ERROR;
          }
#endif /* EK_BOUNDARIES */
        }
        else if(ARG0_IS_S("potential")) {
          if(ek_print_vtk_potential(argv[2]) == 0) {
            argc -= 3;
            argv += 3;

            if((err = gather_runtime_errors(interp, err)) != TCL_OK)
              return TCL_ERROR;
            else
              return TCL_OK;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error in electrokinetics print lbforce vtk #string\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
        else if(ARG0_IS_S("pressure")) {
#ifndef EK_REACTION
          Tcl_AppendResult(interp, "Feature EK_REACTION required", (char *) NULL);
          return (TCL_ERROR);
#else
          if(ek_print_vtk_pressure(argv[2]) == 0) {
            argc -= 3;
            argv += 3;

            if((err = gather_runtime_errors(interp, err)) != TCL_OK)
              return TCL_ERROR;
            else
              return TCL_OK;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error in electrokinetics print pressure vtk #string\n", (char *)NULL);
            return TCL_ERROR;
          }
#endif /* EK_PRESSURE */
        }
        else if(ARG0_IS_S("lbforce")) {
          if(ek_print_vtk_lbforce(argv[2]) == 0) {
            argc -= 3;
            argv += 3;

            if((err = gather_runtime_errors(interp, err)) != TCL_OK)
              return TCL_ERROR;
            else
              return TCL_OK;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error in electrokinetics print potential vtk #string\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
        else {
            Tcl_AppendResult(interp, "Unknown feature \"", argv[0], "\" in electrokinetics print\n", (char *)NULL);
            return TCL_ERROR;
        }
      }
      else if(ARG0_IS_S("node")) {
        argc--;
        argv++;
        
        if(argc != 5 || !ARG_IS_S(0, "print") || !ARG_IS_I(1, coord[0]) || !ARG_IS_I(2, coord[1]) || !ARG_IS_I(3, coord[2]) || !ARG_IS_S(4, "velocity") ) {
          Tcl_AppendResult(interp, "Wrong usage of electrokinetics node print velocity\n", (char *)NULL);
          return TCL_ERROR;
        }
        
        argc -= 4;
        argv += 4;
        
        if(ARG0_IS_S("velocity")) {
          if(ek_lb_print_velocity(coord[0], coord[1], coord[2], vectarg) == 0) {
            argc --;
            argv ++;
            
            Tcl_PrintDouble(interp, vectarg[0], double_buffer);
            Tcl_AppendResult(interp, double_buffer, " ", (char *) NULL);
            Tcl_PrintDouble(interp, vectarg[1], double_buffer);
            Tcl_AppendResult(interp, double_buffer, " ", (char *) NULL);
            Tcl_PrintDouble(interp, vectarg[2], double_buffer);
            Tcl_AppendResult(interp, double_buffer, " ", (char *) NULL);

            if((err = gather_runtime_errors(interp, err)) != TCL_OK)
              return TCL_ERROR;
            else
              return TCL_OK;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error in electrokinetics #int print density vtk #string\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if(ARG0_IS_S("T")) {
        argc--;
        argv++;
        
        if(argc < 1 || !ARG0_IS_D(floatarg)) {
          Tcl_AppendResult(interp, "electrokinetics T requires one floating point number as argument\n", (char *)NULL);
          return TCL_ERROR;
        }
        else if(floatarg <= 0) {
          Tcl_AppendResult(interp, "electrokinetics T must be positive\n", (char *)NULL);
          return TCL_ERROR;
        }
        else {
          if(ek_set_T(floatarg) == 0) {
            argc --;
            argv ++;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error setting electrokinetics T\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if(ARG0_IS_S("bjerrum_length")) {
        argc--;
        argv++;
      
        if(!ARG0_IS_D(floatarg)) {
          Tcl_AppendResult(interp, "electrokinetics bjerrum_length requires one floating point number as argument\n", (char *)NULL);
          return TCL_ERROR;
        }
        else if(floatarg <= 0) {
          Tcl_AppendResult(interp, "electrokinetics bjerrum_length must be positive\n", (char *)NULL);
          return TCL_ERROR;
        }
        else {
          if(ek_set_bjerrumlength(floatarg) == 0) {
            argc--;
            argv++;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error setting electrokinetics bjerrum_length\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if(ARG0_IS_S("agrid")) {
        argc--;
        argv++;
      
        if(!ARG0_IS_D(floatarg)) {
          Tcl_AppendResult(interp, "electrokinetics agrid requires one floating point number as argument\n", (char *)NULL);
          return TCL_ERROR;
        }
        else if(floatarg <= 0) {
          Tcl_AppendResult(interp, "electrokinetics agrid must be positive\n", (char *)NULL);
          return TCL_ERROR;
        }
        else {
          if(ek_set_agrid(floatarg) == 0) {
            argc--;
            argv++;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error setting electrokinetics agrid\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if(ARG0_IS_S("viscosity")) {
        if(argc < 2 || !ARG1_IS_D(floatarg)) {
          Tcl_AppendResult(interp, "electrokinetics viscosity requires one floating point number as argument\n", (char *)NULL);
          return TCL_ERROR;
        } else if (floatarg <= 0) {
          Tcl_AppendResult(interp, "electrokinetics viscosity must be positive\n", (char *)NULL);
          return TCL_ERROR;
        }
        else {
          if(ek_set_viscosity(floatarg) == 0) {
            argc -= 2;
            argv += 2;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error setting electrokinetics viscosity\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if(ARG0_IS_S("friction")) {
        if(argc < 2 || !ARG1_IS_D(floatarg)) {
          Tcl_AppendResult(interp, "electrokinetics friction requires one floating point number as argument\n", (char *)NULL);
          return TCL_ERROR;
        }
        else if(floatarg <= 0) {
          Tcl_AppendResult(interp, "electrokinetics friction must be positive\n", (char *)NULL);
          return TCL_ERROR;
        }
        else {
          if (ek_set_friction(floatarg) == 0) {
            argc -= 2;
            argv += 2;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error setting electrokinetics friction\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if(ARG0_IS_S("bulk_viscosity")) {
        if(argc < 2 || !ARG1_IS_D(floatarg)) {
          Tcl_AppendResult(interp, "electrokinetics bulk_viscosity requires one floating point number as argument\n", (char *)NULL);
          return TCL_ERROR;
        }
        else if(floatarg <= 0) {
          Tcl_AppendResult(interp, "electrokinetics bulk_viscosity must be positive\n", (char *)NULL);
          return TCL_ERROR;
        }
        else {
          if(ek_set_bulk_viscosity(floatarg) == 0) {
            argc -= 2;
            argv += 2;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error setting electrokinetics bulk_viscosity\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if(ARG0_IS_S("gamma_odd")) {
        if(argc < 2 || !ARG1_IS_D(floatarg)) {
          Tcl_AppendResult(interp, "electrokinetics gamma_odd requires one floating point number as argument\n", (char *)NULL);
          return TCL_ERROR;
        }
        else if(fabs(floatarg) >= 1) {
          Tcl_AppendResult(interp, "electrokinetics gamma_odd must be smaller than 1\n", (char *)NULL);
          return TCL_ERROR;
        }
        else {
          if(ek_set_gamma_odd(floatarg) == 0) {
            argc -= 2;
            argv += 2;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error setting electrokinetics gamma_odd\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if(ARG0_IS_S("gamma_even")) {
        if(argc < 2 || !ARG1_IS_D(floatarg)) {
          Tcl_AppendResult(interp, "electrokinetics gamma_even requires one floating point number as argument\n", (char *)NULL);
          return TCL_ERROR;
        }
        else if(fabs(floatarg) >= 1) {
          Tcl_AppendResult(interp, "electrokinetics gamma_even must be smaller than 1\n", (char *)NULL);
          return TCL_ERROR;
        }
        else {
          if(ek_set_gamma_even(floatarg) == 0) {
            argc -= 2;
            argv += 2;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error setting electrokinetics gamma_even\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if(ARG0_IS_S("accelerated_frame"))
      {
        argc -= 1;
        argv += 1;

        if( ARG_IS_S(0,"off") && argc == 1 )
        {
          argc -= 1;
          argv += 1;

          vectarg[0] = 0.0;
          vectarg[1] = 0.0;
          vectarg[2] = 0.0;

          if( ek_set_accelerated_frame( 0 , -1.0, vectarg ) != 0 )
          {
            Tcl_AppendResult(interp, "Unknown error electrokinetics accelerated_frame off\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
        else if(  ARG_IS_S(0,"on") && ARG_IS_S(1,"boundary_mass_density") && 
                  ARG_IS_D(2,floatarg) && ARG_IS_S(3,"ext_acceleration_force") && 
                  ARG_IS_D(4,vectarg[0]) && ARG_IS_D(5,vectarg[1]) && ARG_IS_D(6,vectarg[2])
                  && argc == 7 )
        {
          argc -= 7;
          argv += 7;

          if ( floatarg > 0.0 ) 
          {
            if( ek_set_accelerated_frame( 1 , floatarg, vectarg ) != 0 )
            {
              Tcl_AppendResult(interp, "Unknown error electrokinetics accelerated_frame on\n", (char *)NULL);
              return TCL_ERROR;
            }
          }
          else
          {
            Tcl_AppendResult(interp, "electrokinetics accelerated_frame boundary_mass_density must be greater than zero\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
        else if( ARG_IS_S(0,"print") && ARG_IS_S(1,"boundary_velocity") && argc == 2 )
        {
          argc -= 2;
          argv += 2;

          if( ek_accelerated_frame_print_boundary_velocity( vectarg ) != 0 )
          {
            Tcl_AppendResult(interp, "Unknown error electrokinetics accelerated_frame print\n", (char *)NULL);
            return TCL_ERROR;
          }

          for (int i = 0; i < 3; i++) 
          {
            Tcl_PrintDouble(interp, vectarg[i], double_buffer);
            Tcl_AppendResult(interp, double_buffer, " ", (char *)NULL);
          }
        }
        else
        {
          Tcl_AppendResult(interp, "electrokinetics accelerated_frame requires the following:\n", (char *)NULL);
          Tcl_AppendResult(interp, "electrokinetics accelerated_frame on\n", (char *)NULL);
          Tcl_AppendResult(interp, "                  boundary_mass_density #double\n", (char *)NULL);
          Tcl_AppendResult(interp, "                  ext_acceleration_force #double #double #double\n", (char *)NULL);
          Tcl_AppendResult(interp, "electrokinetics accelerated_frame off\n", (char *)NULL);
          Tcl_AppendResult(interp, "electrokinetics accelerated_frame print boundary_velocity\n", (char *)NULL);
          return TCL_ERROR;
        }
      }
      else if(ARG0_IS_S("reaction")) {
#ifndef EK_REACTION
  Tcl_AppendResult(interp, "EK_REACTION needs to be compiled in to use electrokinetics reaction\n", (char *)NULL);
  return TCL_ERROR;
#else
        argc--;
        argv++;

        if( argc < 16 ) 
        {
          Tcl_AppendResult(interp, "electrokinetics reaction requires 16 arguments: 8 quantifiers, 3 ints, and 7 floats\n", (char *)NULL);
          return TCL_ERROR;
        }
        else
        {

          int counter = 0;

          while(argc > 0 && counter < 12 ) 
          {
            counter++;

            if ( ARG_IS_S_EXACT(0,"reactant_index") ) 
            {

              if ( !ARG_IS_I(1, reactant) ) 
              {
                Tcl_AppendResult(interp, "electrokinetics reactant_index requires one int as argument\n", (char *)NULL);
                return TCL_ERROR;
              }

              if ( (reactant < 0 || reactant > MAX_NUMBER_OF_SPECIES) ) 
              {
                sprintf(int_buffer, "%d", MAX_NUMBER_OF_SPECIES);
                Tcl_AppendResult(interp, "electrokinetics reactant_index #int requires a number between 0 and", int_buffer, "denoting the species\n", (char *)NULL);
                return TCL_ERROR;
              }

              argc -= 2;
              argv += 2;
            }

            else if( ARG_IS_S_EXACT(0,"product0_index") )
            {

              if ( !ARG_IS_I(1, product0) ) 
              {
                Tcl_AppendResult(interp, "electrokinetics product0_index requires one int as argument\n", (char *)NULL);
                return TCL_ERROR;
              }

              if ( (product0 < 0 || product0 > MAX_NUMBER_OF_SPECIES) ) 
              {
                sprintf(int_buffer, "%d", MAX_NUMBER_OF_SPECIES);
                Tcl_AppendResult(interp, "electrokinetics product0_index #int requires a number between 0 and", int_buffer, "denoting the species\n", (char *)NULL);
                return TCL_ERROR;
              }

              argc -= 2;
              argv += 2;
            }

            else if ( ARG_IS_S_EXACT(0,"product1_index") ) 
            {

              if ( !ARG_IS_I(1, product1) ) 
              {
                Tcl_AppendResult(interp, "electrokinetics product1_index requires one int as argument\n", (char *)NULL);
                return TCL_ERROR;
              }

              if ( (product1 < 0 || product1 > MAX_NUMBER_OF_SPECIES) ) 
              {
                sprintf(int_buffer, "%d", MAX_NUMBER_OF_SPECIES);
                Tcl_AppendResult(interp, "electrokinetics product1_index #int requires a number between 0 and", int_buffer, "denoting the species\n", (char *)NULL);
                return TCL_ERROR;
              }

              argc -= 2;
              argv += 2;
            }

            else if ( ARG_IS_S_EXACT(0,"reactant_resrv_density") ) 
            {

              if ( !ARG_IS_D(1, rho_reactant_reservoir) ) 
              {
                Tcl_AppendResult(interp, "electrokinetics reactant_resrv_density requires one floating point number as argument\n", (char *)NULL);
                return TCL_ERROR;
              }

              if ( rho_reactant_reservoir < 0 ) 
              {
                Tcl_AppendResult(interp, "the reactant reservoir density has to be greater than zero\n", (char *)NULL);
                return TCL_ERROR;
              }

              argc -= 2;
              argv += 2;
            }

            else if (ARG_IS_S_EXACT(0,"product0_resrv_density")) 
            {

              if ( !ARG_IS_D(1, rho_product0_reservoir) ) 
              {
                Tcl_AppendResult(interp, "electrokinetics product0_resrv_density requires one floating point number as argument\n", (char *)NULL);
                return TCL_ERROR;
              }

              if ( rho_product0_reservoir < 0 ) 
              {
                Tcl_AppendResult(interp, "the product0 reservoir density has to be greater than zero\n", (char *)NULL);
                return TCL_ERROR;
              }

              argc -= 2;
              argv += 2;
            }

            else if ( ARG_IS_S_EXACT(0,"product1_resrv_density") ) 
            {

              if (!ARG_IS_D(1, rho_product1_reservoir)) 
              {
                Tcl_AppendResult(interp, "electrokinetics product1_resrv_density requires one floating point number as argument\n", (char *)NULL);
                return TCL_ERROR;
              }

              if ( rho_product1_reservoir < 0 ) 
              {
                Tcl_AppendResult(interp, "the product1 reservoir density has to be greater than zero\n", (char *)NULL);
                return TCL_ERROR;
              }

              argc -= 2;
              argv += 2;
            }

            else if ( ARG_IS_S_EXACT(0,"reaction_rate") ) 
            {

              if ( !ARG_IS_D(1, ct_rate) ) 
              {
                Tcl_AppendResult(interp, "electrokinetics reaction_rate requires one floating point number as argument\n", (char *)NULL);
                return TCL_ERROR;
              }

              if ( ct_rate < 0 ) 
              {
                Tcl_AppendResult(interp, "catalytic reaction rate has to be greater than zero\n", (char *)NULL);
                return TCL_ERROR;
              }

              argc -= 2;
              argv += 2;
            }

            else if (ARG_IS_S_EXACT(0,"reaction_radius")) 
            {

              if ( !ARG_IS_D(1, radius) )
              {
                Tcl_AppendResult(interp, "electrokinetics reaction_radius requires one floating point number as argument\n", (char *)NULL);
                return TCL_ERROR;
              }

              if ( radius < 0 ) 
              {
                Tcl_AppendResult(interp, "the catalytic reaction radius has to be greater than zero\n", (char *)NULL);
                return TCL_ERROR;
              }

              argc -= 2;
              argv += 2;
            }

            else if (ARG_IS_S_EXACT(0,"reaction_fraction_pr_0")) 
            {

              if ( !ARG_IS_D(1, fraction0) )
              {
                Tcl_AppendResult(interp, "electrokinetics reaction_fraction_pr_0 requires one floating point number as argument\n", (char *)NULL);
                return TCL_ERROR;
              }

              if ( fraction0 < 0 ) 
              {
                Tcl_AppendResult(interp, "the fraction has to be greater than zero\n", (char *)NULL);
                return TCL_ERROR;
              }

              argc -= 2;
              argv += 2;
            }

            else if (ARG_IS_S_EXACT(0,"reaction_fraction_pr_1")) 
            {

              if ( !ARG_IS_D(1, fraction1) )
              {
                Tcl_AppendResult(interp, "electrokinetics reaction_fraction_pr_1 requires one floating point number as argument\n", (char *)NULL);
                return TCL_ERROR;
              }

              if ( fraction1 < 0 ) 
              {
                Tcl_AppendResult(interp, "the fraction has to be greater than zero\n", (char *)NULL);
                return TCL_ERROR;
              }

              argc -= 2;
              argv += 2;
            }
          }

          if ( counter == 12 )
          {
            Tcl_AppendResult(interp, "unknown option given, please check for spelling errors\n", (char *)NULL);
            return TCL_ERROR;
          }

          if ( (reactant == product0) ||
               (product0 == product1) ||
               (product1 == reactant) ) 
          {
            Tcl_AppendResult(interp, "the reactant and product species need to be distinct\n", (char *)NULL);
            return TCL_ERROR;
          }          

          if ( ek_set_reaction( reactant, product0, product1, rho_reactant_reservoir, rho_product0_reservoir, rho_product1_reservoir, ct_rate, radius, fraction0, fraction1 ) != 0 ) 
          {
            Tcl_AppendResult(interp, "species are not set before invoking electrokinetics reaction\n", (char *)NULL);
            return TCL_ERROR;
          }  
        }     
#endif 
      }
      else {
    	  Tcl_AppendResult(interp, "unknown feature \"", argv[0],"\" of electrokinetics\n", (char *)NULL);
    	  return TCL_ERROR ;
      }
    }
  }

  if((err = gather_runtime_errors(interp, err)) != TCL_OK)
    return TCL_ERROR;
  
  err = ek_init();
  
  switch(err) {
    case 0:
      return TCL_OK;
      
    case 2:
      Tcl_AppendResult(interp, "electrokinetics agrid not compatible with box size\n", (char *)NULL);
      
    default:
      sprintf(int_buffer, "%d", err);
      Tcl_AppendResult(interp, "Error ", int_buffer, " during initialization of electrokinetics", (char *)NULL);
      break;
  }
  
  return TCL_ERROR;
  
#endif /* defined ELECTROKINETICS */
}
