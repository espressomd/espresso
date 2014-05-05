/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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
/** \file lb_tcl.cpp
 *
 * TCL Interface for the Lattice Boltzmann algorithm for hydrodynamic degrees of freedom.
 *
 */


#include "thermostat.hpp"
#include "lb_tcl.hpp"
#include "lb.hpp"
#include "lbgpu.hpp"
#include "parser.hpp"
#include "electrokinetics.hpp"

#ifdef LB_GPU

#ifdef SHANCHEN
int tclprint_to_result_affinityIA(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  Tcl_PrintDouble(interp, data->affinity[0], buffer);
  Tcl_AppendResult(interp, "affinity ", buffer, " ", (char *) NULL);

  for(int ii=1;ii<LB_COMPONENTS;ii++)
  {
     Tcl_PrintDouble(interp, data->affinity[ii], buffer);
     Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  }

  return TCL_OK;
}

int tclcommand_inter_parse_affinity(Tcl_Interp * interp,int part_type_a,int part_type_b, int argc, char ** argv)
{
  double affinity[LB_COMPONENTS];

  if (argc != LB_COMPONENTS+1 )
  {
    Tcl_AppendResult(interp, "Not enough values for affinity", (char *) NULL);      
    return 0;
  }

  for (int ii=0;ii<LB_COMPONENTS;ii++)
  { 
    if ( ! ARG_IS_D(ii+1, (affinity[ii]) ))
    { 
      Tcl_AppendResult(interp, "list of doubles expected for affinity", (char *) NULL);
      return 0;
    } 
  }

  if(affinity_set_params(part_type_a,part_type_b,affinity) == ES_ERROR)
  { 
    Tcl_AppendResult(interp, "Error setting affinity, values must lie between 0 and 1", (char *) NULL);
    return 0;
  }

  Tcl_AppendResult(interp, "Error setting affinity", (char *) NULL);
  return 3;
}
#endif

static int lbnode_parse_set(Tcl_Interp *interp, int argc, char **argv, int *ind)
{
  double f[3];
  
  while (argc > 0)
  {
    if(ARG0_IS_S_EXACT("force"))
    {
      if ( argc < 4 ||
           !ARG_IS_D(1, f[0]) ||
           !ARG_IS_D(2, f[1]) ||
           !ARG_IS_D(3, f[2])
         ) 
      {
        Tcl_AppendResult(interp, "force expects three doubles as argument", (char *)NULL);
        return TCL_ERROR;
      }

      argc -= 4;
      argv += 4;

      if (argc > 0) 
      {
        Tcl_ResetResult(interp);
        Tcl_AppendResult(interp, "Error in lbnode_extforce force. You can only change one field at the same time.", (char *)NULL);
        return ES_ERROR;
      }
    }
    else 
    {
      Tcl_AppendResult(interp, "unknown parameter \"", argv[0], "\" to set", (char *)NULL);
      return TCL_ERROR;
    }
  }

  if (lb_lbnode_set_extforce_GPU(ind, f) == ES_ERROR)
  {
    Tcl_AppendResult(interp, "position is not in the LB lattice", (char *)NULL);
    return TCL_ERROR;
  }

  return ES_OK;
}

/** Parser for the \ref tclcommand_lbnode_extforce_gpu command. Can be used in future to set more values like rho,u e.g.
*/
int tclcommand_lbnode_extforce_gpu(ClientData data, Tcl_Interp *interp, int argc, char **argv) {

  int err=ES_ERROR;
  int coord[3];

  --argc; ++argv;
  
  if (argc < 3) 
  {
    Tcl_AppendResult(interp, "too few arguments for lbnode_extforce", (char *)NULL);
    return ES_ERROR;
  }

  if (!ARG_IS_I(0,coord[0]) || !ARG_IS_I(1,coord[1]) || !ARG_IS_I(2,coord[2]))
  {
    Tcl_AppendResult(interp, "wrong arguments for lbnode", (char *)NULL);
    return ES_ERROR;
  } 
  argc-=3; argv+=3;

  if (argc == 0 ) 
  { 
    Tcl_AppendResult(interp, "lbnode_extforce syntax: lbnode_extforce X Y Z [ print | set ] [ F(X) | F(Y) | F(Z) ]", (char *)NULL);
    return ES_ERROR;
  }

  if (ARG0_IS_S_EXACT("set")) 
  {
    err = lbnode_parse_set(interp, argc-1, argv+1, coord);
  }
  else 
  {
    Tcl_AppendResult(interp, "unknown feature \"", argv[0], "\" of lbnode_extforce", (char *)NULL);
    return  ES_ERROR;
  }    
 
  return err;
}
#endif /* LB_GPU */

#if defined (LB) || defined (LB_GPU)

/* ********************* TCL Interface part *************************************/
/* ******************************************************************************/

void lbfluid_tcl_print_usage(Tcl_Interp *interp)
{
  Tcl_AppendResult(interp, "Usage of \"lbfluid\":\n", (char *)NULL);
  Tcl_AppendResult(interp, "lbfluid [ agrid #float ] [ dens #float ] [ visc #float ] [ tau #tau ]\n", (char *)NULL);
#ifdef SHANCHEN
  Tcl_AppendResult(interp, "        [ mobility #float ]\n", (char *)NULL);
#endif 
  Tcl_AppendResult(interp, "        [ bulk_visc #float ] [ friction #float ] [ gamma_even #float ] [ gamma_odd #float ]\n", (char *)NULL);
  Tcl_AppendResult(interp, "        [ ext_force #float #float #float ]\n", (char *)NULL);
#ifdef SHANCHEN
  Tcl_AppendResult(interp, "        [ coupling #float ]\n", (char *)NULL);
#endif
}

void lbnode_tcl_print_usage(Tcl_Interp *interp) 
{
  Tcl_AppendResult(interp, "lbnode syntax:\n", (char *)NULL);
  Tcl_AppendResult(interp, "lbnode X Y Z print [ rho | u | pi | pi_neq | boundary | populations ]\n", (char *)NULL);
  Tcl_AppendResult(interp, "     or\n", (char *)NULL);
  Tcl_AppendResult(interp, "lbnode X Y Z set [ rho | u | populations ] #nofloats", (char *)NULL);
}

#endif /* LB || LB_GPU */

#if defined (LB) || defined (LB_GPU)
int tclcommand_lbfluid_print_interpolated_velocity(Tcl_Interp *interp, int argc, char **argv);
#endif

/** TCL Interface: The \ref lbfluid command. */
int tclcommand_lbfluid(ClientData data, Tcl_Interp *interp, int argc, char **argv)
{

#ifdef ELECTROKINETICS
  if ( ek_initialized ) 
  {
    Tcl_AppendResult(interp, "ERROR: Electrokinetics automatically intializes the LB on the GPU and can therefore not be used in conjunction with LB.\n");
    Tcl_AppendResult(interp, "ERROR: Please run either electrokinetics or LB.\n");
    
    return TCL_ERROR;
  }
#endif

#if defined (LB) || defined (LB_GPU)
  argc--; argv++;
  
/**if we have LB the LB cpu is set by default */
#ifdef LB
  if( !(lattice_switch & LATTICE_LB_GPU) ) 
    lattice_switch = lattice_switch | LATTICE_LB;
#else
  lattice_switch = lattice_switch | LATTICE_LB_GPU;
#endif

  int err = TCL_OK;
  double floatarg;
  double vectarg[3+LB_COMPONENTS + (LB_COMPONENTS*(LB_COMPONENTS-1))/2];

  if (argc < 1)
  {
    lbfluid_tcl_print_usage(interp);
    return TCL_ERROR;
  }
  else if (ARG0_IS_S_EXACT("off")) 
  {
    lbfluid_tcl_print_usage(interp);
    return TCL_ERROR;
  }
  else if (ARG0_IS_S_EXACT("init")) 
  {
    lbfluid_tcl_print_usage(interp);
    return TCL_ERROR;
  }
  else
  {
    while (argc > 0) 
    {
      if (ARG0_IS_S_EXACT("gpu") || ARG0_IS_S_EXACT("GPU")) 
      {
#ifdef LB_GPU
        lattice_switch = (lattice_switch &~ LATTICE_LB) | LATTICE_LB_GPU;
        argc--; argv++;
#else
        Tcl_AppendResult(interp, "LB_GPU is not compiled in!", NULL);
        return TCL_ERROR;
#endif
      }
      else if (ARG0_IS_S_EXACT("cpu") || ARG0_IS_S_EXACT("CPU")) 
      {
#ifdef LB
        lattice_switch = (lattice_switch & ~LATTICE_LB_GPU) | LATTICE_LB;
        argc--; argv++;
#else
        Tcl_AppendResult(interp, "LB is not compiled in!", NULL);
        return TCL_ERROR;
#endif
      }
      else if ( ARG0_IS_S_EXACT("grid") || ARG0_IS_S_EXACT("agrid") )
      {
        if ( argc < 2 || !ARG1_IS_D(floatarg) ) 
        {
          Tcl_AppendResult(interp, "agrid requires 1 argument", (char *)NULL);
          return TCL_ERROR;
        }
        else if (floatarg <= 0) 
        {
          Tcl_AppendResult(interp, "agrid must be positive", (char *)NULL);
          return TCL_ERROR;
        }
        else 
        {
          if ( lb_lbfluid_set_agrid(floatarg) == 0 ) 
          {
            argc-=2; argv+=2;
          } 
          else 
          {
            Tcl_AppendResult(interp, "Unknown Error setting agrid", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if (ARG0_IS_S_EXACT("tau")) 
      {
        if ( argc < 2 || !ARG1_IS_D(floatarg) ) 
        {
          Tcl_AppendResult(interp, "tau requires 1 argument", (char *)NULL);
          return TCL_ERROR;
        } 
        else if (floatarg <= 0) 
        {
          Tcl_AppendResult(interp, "tau must be positive", (char *)NULL);
          return TCL_ERROR;
        }
        else if (floatarg < time_step ) 
        {
          Tcl_AppendResult(interp, "tau must larger than the MD time step", (char *)NULL);
          return TCL_ERROR;
        } 
        else 
        {
          if ( lb_lbfluid_set_tau(floatarg) == 0 ) 
          {
            argc-=2; argv+=2;
          } 
          else 
          {
            Tcl_AppendResult(interp, "Unknown Error setting tau", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
#ifdef SHANCHEN
      else if (ARG0_IS_S_EXACT("sc_coupling") ) 
      {
        /* when SHANCHEN==1 we allow liquid/gas phase transitions and require a second parameter (rho_0) besides the coupling */

        int nargs=( (LB_COMPONENTS==1) ? 2 : (LB_COMPONENTS*(LB_COMPONENTS+1))/2);

        if ( argc < 1+nargs ) 
        {
          char str[1024];
          sprintf(str,"sc_coupling requires %d arguments",nargs);
          Tcl_AppendResult(interp, str, (char *)NULL);
          return TCL_ERROR;
        }
        else 
        {
          int i;

          for(i=0;i<nargs;i++)
          {
            if( !ARG_IS_D(i+1,vectarg[i]) ) 
            {
              Tcl_AppendResult(interp, "sc_coupling requires real numbers as arguments", (char *)NULL); 
              return TCL_ERROR;
            }
          }

          if ( lb_lbfluid_set_shanchen_coupling(vectarg) == 0 ) 
          {
            argc-=1+nargs; argv+=1+nargs;
          } 
          else 
          {
            Tcl_AppendResult(interp, "Unknown Error setting coupling", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if (ARG0_IS_S_EXACT("mobility")) 
      {
        argc--; argv++;

        if ( argc < LB_COMPONENTS -1 ) 
        {
          Tcl_AppendResult(interp, "mobility requires argument(s)", (char *)NULL);
          return TCL_ERROR;
        } 
        else 
        { 
          int i;

          for(i=0;i<LB_COMPONENTS-1;i++)
          {
            if(!ARG_IS_D(i,vectarg[i]) ) 
            {
              Tcl_AppendResult(interp, "mobility requires real numbers as arguments", (char *)NULL); // TODO: fix this and similar ones...
              return TCL_ERROR;
            }
            else if (vectarg[i]<=0)
            { 
              Tcl_AppendResult(interp, "mobility must be positive", (char *)NULL);
              return TCL_ERROR;
            }
          } 

          if ( lb_lbfluid_set_mobility(vectarg) == 0 ) 
          {
            argc-=(LB_COMPONENTS-1); argv+=(LB_COMPONENTS-1);
          }
          else
          {
            Tcl_AppendResult(interp, "Unknown Error setting mobility", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if (ARG0_IS_S_EXACT("remove_momentum")) 
      {
        if ( lb_lbfluid_set_remove_momentum() == 0 ) 
        {
          argc-=2; argv+=2;
        } 
        else 
        {
          Tcl_AppendResult(interp, "Unknown Error setting remove_momentum", (char *)NULL);
          return TCL_ERROR;
        }
      }
#endif // SHANCHEN
      else if (ARG0_IS_S_EXACT("density") || ARG0_IS_S_EXACT("dens")) 
      {
        if ( argc < LB_COMPONENTS + 1 ) 
        {
          Tcl_AppendResult(interp, "dens requires \"", LB_COMPONENTS  ,"\"argument(s)", (char *)NULL);
          return TCL_ERROR;
        }
        else 
        {      
          int i;

          for(i=0;i<LB_COMPONENTS;i++)
          {
            if(!ARG_IS_D(i+1,vectarg[i]) ) 
            {
              Tcl_AppendResult(interp, "dens requires real numbers as arguments", (char *)NULL); // TODO: fix this and similar ones...
              return TCL_ERROR;
            }
            else if (vectarg[i]<=0) 
            { 
              Tcl_AppendResult(interp, "dens must be positive", (char *)NULL);
              return TCL_ERROR;
            }
          }

          if ( lb_lbfluid_set_density(vectarg) == 0 ) 
          {
            argc-=(1+LB_COMPONENTS); argv+=(1+LB_COMPONENTS);
          } 
          else 
          {
            Tcl_AppendResult(interp, "Unknown Error setting dens", (char *)NULL);
            return TCL_ERROR;
          }
        } 
      }
      else if (ARG0_IS_S_EXACT("viscosity") || ARG0_IS_S_EXACT("visc")) 
      {
        if ( argc < (1+LB_COMPONENTS)  ) 
        {
          Tcl_AppendResult(interp, "visc requires 1 argument", (char *)NULL);
          return TCL_ERROR;
        } 
        else 
        {  
          int i;

          for(i=0;i<LB_COMPONENTS;i++)
          {
            if(!ARG_IS_D(i+1,vectarg[i]) ) 
            {
              Tcl_AppendResult(interp, "visc requires arguments", (char *)NULL); // TODO: fix this and similar ones...
              return TCL_ERROR;
            } 
            else if (vectarg[i]<=0)
            { 
              Tcl_AppendResult(interp, "visc must be positive", (char *)NULL);
              return TCL_ERROR;
            }
          }

          if ( lb_lbfluid_set_visc(vectarg) == 0 ) 
          {
            argc-=(1+LB_COMPONENTS); argv+=(1+LB_COMPONENTS);
          } 
          else 
          {
            Tcl_AppendResult(interp, "Unknown Error setting viscosity", (char *)NULL);
            return TCL_ERROR;
          }
        }
      } 
      else if (ARG0_IS_S_EXACT("bulk_viscosity")) 
      {
        if ( argc < (LB_COMPONENTS+1) )
        { 
          Tcl_AppendResult(interp, "bulk_viscosity requires arguments", (char *)NULL); // TODO: fix this and similar ones...
          return TCL_ERROR;
        } 
        else  
        {
          int i;

          for(i=0;i<LB_COMPONENTS;i++)
          {
            if(!ARG_IS_D(i+1,vectarg[i]) ) 
            {
              Tcl_AppendResult(interp, "bulk_viscosity requires arguments", (char *)NULL); // TODO: fix this and similar ones...
              return TCL_ERROR;
            } 
            else if (vectarg[i]<=0)
            { 
              Tcl_AppendResult(interp, "bulk_viscosity must be positive", (char *)NULL);
              return TCL_ERROR;
            }
          }

          if ( lb_lbfluid_set_bulk_visc(vectarg) == 0 ) 
          {
            argc-=(1+LB_COMPONENTS); argv+=(1+LB_COMPONENTS) ; 
          } 
          else 
          {
            Tcl_AppendResult(interp, "Unknown Error setting bulk_viscosity", (char *)NULL);
            return TCL_ERROR;
          }
        }     
      } 
      else if (ARG0_IS_S_EXACT("friction") ) 
      {
        if ( argc < (LB_COMPONENTS+1) )
        { 
          Tcl_AppendResult(interp, "friction requires arguments", (char *)NULL); // TODO: fix this and similar ones...
          return TCL_ERROR;
        } 
        else 
        {
          int i;

          for(i=0;i<LB_COMPONENTS;i++)
          {
            if(!ARG_IS_D(i+1,vectarg[i]) ) 
            {
              Tcl_AppendResult(interp, "friction requires arguments", (char *)NULL); // TODO: fix this and similar ones...
              return TCL_ERROR;
            } 
            else if (vectarg[i]<=0)
            { 
              Tcl_AppendResult(interp, "friction must be positive", (char *)NULL);
              return TCL_ERROR;
            }
          }

          if ( lb_lbfluid_set_friction(vectarg) == 0 ) 
          {
            argc-=(1+LB_COMPONENTS); argv+=(1+LB_COMPONENTS); 
          } 
          else 
          {
            Tcl_AppendResult(interp, "Unknown Error setting friction", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if (ARG0_IS_S_EXACT("couple") ) 
      {
        if ( argc < 1 ) 
        { 
          Tcl_AppendResult(interp, "couple requires an argument, either 2pt or 3pt", (char *)NULL);
          return TCL_ERROR;
        }
        else 
        {
          if ( ARG1_IS_S_EXACT("2pt") || ARG1_IS_S_EXACT("2PT") || ARG1_IS_S_EXACT("2Pt") ) 
          {
            lb_lbfluid_set_couple_flag (LB_COUPLE_TWO_POINT);
          }
          else if ( ARG1_IS_S_EXACT("3pt") || ARG1_IS_S_EXACT("3PT") || ARG1_IS_S_EXACT("3Pt") ) 
          {
            lb_lbfluid_set_couple_flag (LB_COUPLE_THREE_POINT);
          }
          else
          {
            Tcl_AppendResult(interp, "Did not understand argument to couple, please send 2pt or 3pt.", (char *)NULL);
            return TCL_ERROR;
          }

          argc-=2; argv+=2;
        }
      }
      else if (ARG0_IS_S_EXACT("gamma_odd") ) 
      {
        if ( argc < (LB_COMPONENTS+1) )
        { 
          Tcl_AppendResult(interp, "gamma_odd requires arguments", (char *)NULL); // TODO: fix this and similar ones...
          return TCL_ERROR;
        }
        else  
        {
          int i;

          for(i=0;i<LB_COMPONENTS;i++)
          {
            if(!ARG_IS_D(i+1,vectarg[i]) ) 
            {
              Tcl_AppendResult(interp, "gamma_odd requires arguments", (char *)NULL); // TODO: fix this and similar ones...
              return TCL_ERROR;
            } 
            else if (vectarg[i] >= 1)
            { 
              Tcl_AppendResult(interp, "gamma_odd must be < 1", (char *)NULL);
              return TCL_ERROR;
            }
          }

          if ( lb_lbfluid_set_gamma_odd(vectarg) == 0 ) 
          {
            argc-=(1+LB_COMPONENTS); argv+=(1+LB_COMPONENTS); 
          }
          else
          {
            Tcl_AppendResult(interp, "Unknown Error setting gamma_odd", (char *)NULL);
            return TCL_ERROR;
          }
        }     
      } 
      else if (ARG0_IS_S_EXACT("gamma_even") ) 
      {
        if ( argc < (LB_COMPONENTS+1) )
        { 
          Tcl_AppendResult(interp, "gamma_even requires arguments", (char *)NULL); // TODO: fix this and similar ones...
          return TCL_ERROR;
        } 
        else  
        {
          int i;

          for(i=0;i<LB_COMPONENTS;i++)
          {
            if(!ARG_IS_D(i+1,vectarg[i]) )
            {
              Tcl_AppendResult(interp, "gamma_even requires arguments", (char *)NULL); // TODO: fix this and similar ones...
              return TCL_ERROR;
            }
            else if (vectarg[i] >= 1)
            { 
              Tcl_AppendResult(interp, "gamma_even must be < 1 ", (char *)NULL);
              return TCL_ERROR;
            }
          }

          if ( lb_lbfluid_set_bulk_visc(vectarg) == 0 )
          {
            argc-=(1+LB_COMPONENTS); argv+=(1+LB_COMPONENTS); 
          }
          else
          {
            Tcl_AppendResult(interp, "Unknown Error setting dens", (char *)NULL);
            return TCL_ERROR;
          }
        }     
      } 
      else if (ARG0_IS_S_EXACT("ext_force")) 
      {
#ifdef EXTERNAL_FORCES
        int i;

        if ( argc <= 3*LB_COMPONENTS ) 
        {
              Tcl_AppendResult(interp, "ext_force requires 3 arguments per component", (char *)NULL);
              return TCL_ERROR;
        } 
        for(i=0;i<LB_COMPONENTS;i++)
        {
            if ( !ARG_IS_D(1, vectarg[0]) || !ARG_IS_D(2, vectarg[1]) ||  !ARG_IS_D(3, vectarg[2]) ) 
            {
              Tcl_AppendResult(interp, "ext_force requires 3 real numbers per component", (char *)NULL);
              return TCL_ERROR;
            } 
            else if (lb_lbfluid_set_ext_force(i,vectarg[0], vectarg[1], vectarg[2]) == 0) 
            {
              argc-=3; argv+=3;
            } 
            else 
            {
              Tcl_AppendResult(interp, "Unknown Error setting ext_force", (char *)NULL);
              return TCL_ERROR;
            }
        }
        argc-=1; argv+=1;
#else
        Tcl_AppendResult(interp, "External Forces not compiled in!", (char *)NULL);
        return TCL_ERROR;
#endif
      }
      else if (ARG0_IS_S_EXACT("print"))
      {
        if ( argc < 3 || (ARG1_IS_S_EXACT("vtk") && argc < 4) )
        {
          Tcl_AppendResult(interp, "lbfluid print requires at least 2 arguments. Usage: lbfluid print [vtk] velocity|boundary filename", (char *)NULL);
          return TCL_ERROR;
        }
        else 
        {
          argc--;
          argv++;

          if (ARG0_IS_S_EXACT("vtk")) 
          {
            if (ARG1_IS_S_EXACT("boundary")) 
            {
              if ( lb_lbfluid_print_vtk_boundary(argv[2]) != 0 ) 
              {
                Tcl_AppendResult(interp, "Unknown Error at lbfluid print vtk boundary", (char *)NULL);
                return TCL_ERROR;
              }

              argc -= 3;
              argv += 3;
            }
            else if (ARG1_IS_S_EXACT("velocity")) 
            {
              if ( lb_lbfluid_print_vtk_velocity(argv[2]) != 0 ) 
              {
                Tcl_AppendResult(interp, "Unknown Error at lbfluid print vtk velocity", (char *)NULL);
                return TCL_ERROR;
              }

              argc -= 3;
              argv += 3;
            }
            else if (ARG1_IS_S_EXACT("density")) 
            {
              argc -= 2;
              argv += 2;

              if ( argc < LB_COMPONENTS ) 
              {
                Tcl_AppendResult(interp, "lbfluid print vtk density requires\"", LB_COMPONENTS,"\" arguments.", (char *)NULL);
                return TCL_ERROR;
              }

              if ( lb_lbfluid_print_vtk_density(argv) != 0 ) 
              {
                Tcl_AppendResult(interp, "Unknown Error at lbfluid print vtk density", (char *)NULL);
                return TCL_ERROR;
              }

              argc -= LB_COMPONENTS;
              argv += LB_COMPONENTS;
            }
            else 
            {
              return TCL_ERROR;
            }
          }
          else
          { 
            // SAW TODO : finish implementing for SHANCHEN

            if (ARG0_IS_S_EXACT("boundary")) 
            {
              if ( lb_lbfluid_print_boundary(argv[1]) != 0 ) 
              {
                Tcl_AppendResult(interp, "Unknown Error at lbfluid print boundary", (char *)NULL);
                return TCL_ERROR;
              }
            }
            else if (ARG0_IS_S_EXACT("velocity")) 
            {
              if ( lb_lbfluid_print_velocity(argv[1]) != 0 ) 
              {
                Tcl_AppendResult(interp, "Unknown Error at lbfluid print velocity", (char *)NULL);
                return TCL_ERROR;
              }
            }
            else 
            {
              return TCL_ERROR;
            }

            argc -= 2;
            argv += 2;
          }
        }
      }
      else if (ARG0_IS_S_EXACT("save_ascii_checkpoint")) 
      { 
        if (argc < 2) 
        {
          Tcl_AppendResult(interp, "usage: lbfluid save_ascii_checkpoint <filename>", (char *)NULL);
          return TCL_ERROR;
        } 
        else 
        {
          return lb_lbfluid_save_checkpoint(argv[1], 0);
        }
      }  
      else if (ARG0_IS_S_EXACT("save_binary_checkpoint")) 
      { 
        if (argc < 2) 
        {
          Tcl_AppendResult(interp, "usage: lbfluid save_binary_checkpoint <filename>", (char *)NULL);
          return TCL_ERROR;
        } 
        else 
        {
          return lb_lbfluid_save_checkpoint(argv[1], 1);
        }
      }  
      else if (ARG0_IS_S_EXACT("load_ascii_checkpoint")) 
      { 
        if (argc < 2) 
        {
          Tcl_AppendResult(interp, "usage: lbfluid load_ascii_checkpoint <filename>", (char *)NULL);
          return TCL_ERROR;
        } 
        else 
        {
          return lb_lbfluid_load_checkpoint(argv[1], 0);
        }
      }  
      else if (ARG0_IS_S_EXACT("load_binary_checkpoint")) 
      { 
        if (argc < 2) 
        {
          Tcl_AppendResult(interp, "usage: lbfluid load_binary_checkpoint <filename>", (char *)NULL);
          return TCL_ERROR;
        } 
        else 
        {
          return lb_lbfluid_load_checkpoint(argv[1], 1);
        }
      }  
#if defined(LB) || defined(LB_GPU)
      else if (ARG0_IS_S_EXACT("print_interpolated_velocity")) 
      { 
        //this has to come after print
        return tclcommand_lbfluid_print_interpolated_velocity(interp, argc-1, argv+1);
      }
#endif
      else
      {
        Tcl_AppendResult(interp, "unknown feature \"", argv[0],"\" of lbfluid", (char *)NULL);
        return TCL_ERROR ;
      }

      if ((err = gather_runtime_errors(interp, err)) != TCL_OK)
      return TCL_ERROR;
    } /* END OF WHILE LOOP */
  } /* END OF ELSE */

  mpi_bcast_parameter(FIELD_LATTICE_SWITCH);

  /* thermo_switch is retained for backwards compatibility */
  thermo_switch = (thermo_switch | THERMO_LB);
  mpi_bcast_parameter(FIELD_THERMO_SWITCH);

  return TCL_OK;
#else /* !defined LB */
  Tcl_AppendResult(interp, "LB is not compiled in!", NULL);
  return TCL_ERROR;
#endif
}

/** Parser for the \ref tclcommand_lbnode command. */
int tclcommand_lbnode(ClientData data, Tcl_Interp *interp, int argc, char **argv)
{
#if defined (LB) || defined (LB_GPU)
  int coord[3];
  int counter;
  double double_return[19*LB_COMPONENTS];
  char double_buffer[TCL_DOUBLE_SPACE];

  for (counter = 0; counter < 19; counter++) 
  double_return[counter]=0;

  --argc; ++argv;

  if (lattice_switch & LATTICE_LB_GPU) 
  {
    // Do nothing
  } 
  else 
  {
#ifdef LB
    if (lbfluid[0][0]==0) 
    {
      Tcl_AppendResult(interp, "lbnode: lbfluid not correctly initialized", (char *)NULL);
      return TCL_ERROR;
    }
#endif
  }

  if (argc < 3) 
  {
    lbnode_tcl_print_usage(interp);
    return TCL_ERROR;
  }

  if (!ARG_IS_I(0,coord[0]) || !ARG_IS_I(1,coord[1]) || !ARG_IS_I(2,coord[2])) 
  {
    Tcl_AppendResult(interp, "Coordinates are not integer.", (char *)NULL);
    return TCL_ERROR;
  } 

  if (lattice_switch & LATTICE_LB_GPU) 
  {
#ifdef LB_GPU
    if (coord[0]<0 || coord[0]>(box_l[0])/lbpar_gpu.agrid-1 || coord[1]<0 || coord[1]>(box_l[1])/lbpar_gpu.agrid-1 || coord[2]<0 || coord[2]>(box_l[2])/lbpar_gpu.agrid-1)
    {
      Tcl_AppendResult(interp, "Coordinates do not correspond to a valid LB node index", (char *)NULL);
      return TCL_ERROR;
    }
#endif
  } 
  else 
  {
#ifdef LB
  if (coord[0]<0 || coord[0]>(box_l[0])/lbpar.agrid-1 || coord[1]<0 || coord[1]>(box_l[1])/lbpar.agrid-1 || coord[2]<0 || coord[2]>(box_l[2])/lbpar.agrid-1) 
  {
    Tcl_AppendResult(interp, "Coordinates do not correspond to a valid LB node index", (char *)NULL);
    return TCL_ERROR;
  } 
#endif
  }

  argc-=3;
  argv+=3;

  if (ARG0_IS_S_EXACT("print")) 
  {
    argc--;
    argv++;

    while (argc > 0) 
    {
      if (ARG0_IS_S_EXACT("rho") || ARG0_IS_S_EXACT("density")) 
      {
        lb_lbnode_get_rho(coord, double_return);

        for (counter = 0; counter < LB_COMPONENTS; counter++) 
        {
          Tcl_PrintDouble(interp, double_return[counter], double_buffer);
          Tcl_AppendResult(interp, double_buffer, " ", (char *)NULL);
        }

        argc--; argv++;
      }
      else if (ARG0_IS_S_EXACT("u") || ARG0_IS_S_EXACT("v") || ARG0_IS_S_EXACT("velocity")) 
      {
        lb_lbnode_get_u(coord, double_return);

        for (counter = 0; counter < 3; counter++) 
        {
          Tcl_PrintDouble(interp, double_return[counter], double_buffer);
          Tcl_AppendResult(interp, double_buffer, " ", (char *)NULL);
        }

        argc--; argv++;
      }
      else if (ARG0_IS_S_EXACT("pi") || ARG0_IS_S_EXACT("pressure")) 
      {
        lb_lbnode_get_pi(coord, double_return);

        for (counter = 0; counter < 6; counter++) 
        {
          Tcl_PrintDouble(interp, double_return[counter], double_buffer);
          Tcl_AppendResult(interp, double_buffer, " ", (char *)NULL);
        }

        argc--; argv++;
      }
      else if (ARG0_IS_S_EXACT("pi_neq")) 
      { 
        /* this has to come after pi */

        lb_lbnode_get_pi_neq(coord, double_return);

        for (counter = 0; counter < 6; counter++) 
        {
          Tcl_PrintDouble(interp, double_return[counter], double_buffer);
          Tcl_AppendResult(interp, double_buffer, " ", (char *)NULL);
        }

        argc--; argv++;
      }
#ifndef SHANCHEN
      else if (ARG0_IS_S_EXACT("boundary")) 
      {
        char integer_buffer[TCL_INTEGER_SPACE];
        int integer_return = 0;

        lb_lbnode_get_boundary(coord, &integer_return);
        sprintf(integer_buffer, "%d", integer_return);
        Tcl_AppendResult(interp, integer_buffer, " ", (char *)NULL);

        argc--; argv++;
      }
#endif // SHANCHEN
      else if (ARG0_IS_S_EXACT("populations") || ARG0_IS_S_EXACT("pop")) 
      { 
        lb_lbnode_get_pop(coord, double_return);

        for (counter = 0; counter < 19; counter++) 
        {
          Tcl_PrintDouble(interp, double_return[counter], double_buffer);
          Tcl_AppendResult(interp, double_buffer, " ", (char *)NULL);
        }

        argc--; argv++;
      }
      else 
      {
        Tcl_ResetResult(interp);
        Tcl_AppendResult(interp, "unknown fluid data \"", argv[0], "\" requested", (char *)NULL);

        return TCL_ERROR;
      }
    }
  }
  else if (ARG0_IS_S_EXACT("set")) 
  {
    argc--;
    argv++;

    if (ARG0_IS_S_EXACT("rho") || ARG0_IS_S_EXACT("density")) 
    {
      argc--; argv++;

      if (argc!=LB_COMPONENTS) 
      {
        char integer_buffer[TCL_INTEGER_SPACE];

        sprintf(integer_buffer, "%d", LB_COMPONENTS);
        Tcl_AppendResult(interp, "LB node set density requires ", integer_buffer, " double(s).", (char *)NULL);

        return TCL_ERROR;
      }

      for (counter = 0; counter < LB_COMPONENTS; counter++) 
      {
        if (!ARG0_IS_D(double_return[counter])) 
        {
          Tcl_AppendResult(interp, "received not a double but \"", argv[0], "\" requested", (char *)NULL);
          return TCL_ERROR;
        }

        argc--; argv++;
      }

      //SAW  TODO: clean code/naming, there is ambiguity with lbfluid_set_dens
      if (lb_lbnode_set_rho(coord, double_return) != 0) 
      {
        Tcl_AppendResult(interp, "General Error on lbnode set rho.", (char *)NULL);
        return TCL_ERROR;
      }
    }
    else if (ARG0_IS_S_EXACT("u") || ARG0_IS_S_EXACT("v") || ARG0_IS_S_EXACT("velocity")) 
    {
      argc--; argv++;

      if(argc!=3) 
      {
        Tcl_AppendResult(interp, "lbnode set velocity|u|v needs three arguments (vx vy vz)", (char *)NULL);
        return TCL_ERROR;
      }

      for (counter = 0; counter < 3; counter++) 
      {
        if (!ARG0_IS_D(double_return[counter])) 
        {
          Tcl_AppendResult(interp, "received not a double but \"", argv[0], "\" requested", (char *)NULL);
          return TCL_ERROR;
        }

        argc--; argv++;
      }

      if (lb_lbnode_set_u(coord, double_return) != 0) 
      {
        Tcl_AppendResult(interp, "General Error on lbnode set u.", (char *)NULL);
        return TCL_ERROR;
      }
    }
    else if (ARG0_IS_S_EXACT("pop") || ARG0_IS_S_EXACT("populations") ) 
    {
      argc--; argv++;

      if (argc!=19) 
        {
        Tcl_AppendResult(interp, "LB node set populations requires 19 doubles.", (char *)NULL);
        return TCL_ERROR;
      }

      for (counter = 0; counter < 19; counter++) 
      {
        if (!ARG0_IS_D(double_return[counter])) 
        {
          Tcl_AppendResult(interp, "received not a double but \"", argv[0], "\" requested", (char *)NULL);
          return TCL_ERROR;
        }

        argc--; argv++;
      }

      if (lb_lbnode_set_pop(coord, double_return) != 0) 
      {
        Tcl_AppendResult(interp, "General Error on lbnode set pop.", (char *)NULL);
        return TCL_ERROR;
      }
    }
    else 
    {
      Tcl_AppendResult(interp, "unknown feature \"", argv[0], "\" of lbnode x y z set", (char *)NULL);
      return  TCL_ERROR;
    }
  } 
  else 
  {
    Tcl_AppendResult(interp, "unknown feature \"", argv[0], "\" of lbnode", (char *)NULL);
    return  TCL_ERROR;
  }

  return TCL_OK;

#else /* !defined LB */
  Tcl_AppendResult(interp, "LB is not compiled in!", NULL);
  return TCL_ERROR;
#endif
}

int tclcommand_lbfluid_print_interpolated_velocity(Tcl_Interp *interp, int argc, char **argv) 
{
#if defined (LB) || defined (LB_GPU)
  double p[3], v[3];
  char buffer[3*TCL_DOUBLE_SPACE+3];

  if (argc!=3) 
  {
    printf("usage: print_interpolated_velocity $x $y $z");
    return TCL_ERROR;
  }

  for (int i = 0; i < 3; i++) 
  {
    if (!ARG_IS_D(i, p[i]))
      printf("usage: print_interpolated_velocity $x $y $z");
  }

  lb_lbfluid_get_interpolated_velocity_global(p, v);

  for (int i = 0; i < 3; i++) 
  {
    Tcl_PrintDouble(interp, v[i], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  }

  return TCL_OK;

#else
  Tcl_AppendResult(interp, "LB is not compiled in!", NULL);
  return TCL_ERROR;
#endif
}
