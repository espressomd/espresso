/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
/** \file interaction_data.cpp
    Implementation of interaction_data.hpp
 */
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <algorithm>
#include "interaction_data_tcl.hpp"
#include "interaction_data.hpp"
#include "communication.hpp"
#include "global.hpp"

#include "comforce_tcl.hpp"
#include "comfixed_tcl.hpp"
#include "rattle_tcl.hpp"
#include "mol_cut_tcl.hpp"
#include "actor/DipolarDirectSum_tcl.hpp" 


// for the force caps
#include "lj.hpp"
#include "ljangle.hpp"
#include "morse.hpp"
#include "tab.hpp"
#include "buckingham.hpp"

//for surface charge output
#include "mmm2d.hpp"
#include "mmm-common.hpp"
#include "elc.hpp"

// Nonbonded
#include "bmhtf-nacl_tcl.hpp"
#include "buckingham_tcl.hpp"
#include "gb_tcl.hpp"
#include "gaussian_tcl.hpp"
#include "hat_tcl.hpp"
#include "lb_tcl.hpp"
#include "lj_tcl.hpp"
#include "ljangle_tcl.hpp"
#include "ljcos_tcl.hpp"
#include "ljcos2_tcl.hpp"
#include "cos2_tcl.hpp"
#include "ljgen_tcl.hpp"
#include "hertzian_tcl.hpp"
#include "morse_tcl.hpp"
#include "dpd_tcl.hpp"
#include "soft_sphere_tcl.hpp"
#include "object-in-fluid/affinity_tcl.hpp"
#include "object-in-fluid/membrane_collision_tcl.hpp"
#include "steppot_tcl.hpp"
#include "tab_tcl.hpp"
#include "tunable_slip_tcl.hpp"

// Coulomb
#include "debye_hueckel_tcl.hpp"
#include "elc_tcl.hpp"
#include "maggs_tcl.hpp"
#include "mmm1d_tcl.hpp"
#include "mmm2d_tcl.hpp"
#include "p3m_tcl.hpp"
#include "reaction_field_tcl.hpp"
#include "actor/Mmm1dgpu_tcl.hpp"
#include "actor/Ewaldgpu_tcl.hpp"
#ifdef SCAFACOS
#include "scafacos.hpp"
#endif

// Magnetostatics
#include "mdlc_correction_tcl.hpp"
#include "p3m-dipolar_tcl.hpp"
#include "magnetic_non_p3m_methods_tcl.hpp"

// Bonded
#include "angle_tcl.hpp"
#include "angle_harmonic_tcl.hpp"
#include "angle_cosine_tcl.hpp"
#include "angle_cossquare_tcl.hpp"
#include "angledist_tcl.hpp"
#include "dihedral_tcl.hpp"
#include "endangledist_tcl.hpp"
#include "fene_tcl.hpp"
#include "overlap_tcl.hpp"
#include "harmonic_dumbbell_tcl.hpp"
#include "harmonic_tcl.hpp"
#include "quartic_tcl.hpp"
#include "bonded_coulomb_tcl.hpp"
#include "subt_lj_tcl.hpp"
#include "umbrella_tcl.hpp"
#include "tcl/object-in-fluid/oif_global_forces_tcl.hpp"
#include "tcl/object-in-fluid/oif_local_forces_tcl.hpp"
#include "tcl/object-in-fluid/out_direction_tcl.hpp"
#ifdef TWIST_STACK
#include "twist_stack_tcl.hpp"
#endif
#ifdef HYDROGEN_BOND
#include "hydrogen_bond_tcl.hpp"
#endif
#include "immersed_boundary/ibm_triel_tcl.hpp"
#include "immersed_boundary/ibm_volume_conservation_tcl.hpp"
#include "immersed_boundary/ibm_tribend_tcl.hpp"

#ifdef DIPOLES
int tclprint_to_result_DipolarIA(Tcl_Interp *interp);
#endif

#ifdef ELECTROSTATICS
int tclprint_to_result_CoulombIA(Tcl_Interp *interp);

/********************************************************************************/
/*                                 electrostatics                               */
/********************************************************************************/

#ifdef SCAFACOS
int tclcommand_scafacos_methods(ClientData data, Tcl_Interp * interp, int argc, char ** argv) {
  std::list<std::string> methods = Electrostatics::Scafacos::available_methods();
  
  Tcl_AppendResult(interp, "{ ", 0);
  for(auto &m: methods) {
    Tcl_AppendResult(interp, m.c_str(), " ", 0);
  }
  Tcl_AppendResult(interp, "}", 0);

  return TCL_OK;
}

int tclcommand_inter_coulomb_parse_scafacos(Tcl_Interp *interp, int argc, char ** argv) {
  if(argc < 1)
    return TCL_ERROR;

  const std::string method(argv[0]);
  
  std::stringstream params;

  if(argc > 1) {
    for(int i = 1; i < argc; i++) {
      params << std::string(argv[i]);
      if(i != (argc-1))
        params << ",";
    }
  }

  coulomb.method  = COULOMB_SCAFACOS;

  mpi_bcast_coulomb_params();
  
  Electrostatics::Scafacos::set_parameters(method, params.str());
  
  return TCL_OK;
}

int tclprint_to_result_scafacos(Tcl_Interp *interp) {
  std::string p = Electrostatics::Scafacos::get_parameters();
  Tcl_AppendResult(interp, p.c_str(), 0);
  return TCL_OK;
}

#endif

int tclcommand_inter_parse_coulomb(Tcl_Interp * interp, int argc, char ** argv)
{
  double d1;

  Tcl_ResetResult(interp);

  if(argc == 0) {
    tclprint_to_result_CoulombIA(interp);
    return TCL_OK;
  }
  
  if (! ARG0_IS_D(d1)) {
#ifdef P3M
    Tcl_ResetResult(interp);
    if (ARG0_IS_S("elc") && ((coulomb.method == COULOMB_P3M) || (coulomb.method == COULOMB_ELC_P3M)))
      return tclcommand_inter_coulomb_parse_elc_params(interp, argc - 1, argv + 1);
    if (ARG0_IS_S("elc") && coulomb.method == COULOMB_P3M_GPU) {
      Tcl_AppendResult(interp, "elc can not be used in conjunction with the gpu p3m",
		       (char *) NULL);
      return TCL_ERROR;
    }
    if (coulomb.method == COULOMB_P3M || coulomb.method == COULOMB_P3M_GPU)
      return tclcommand_inter_coulomb_parse_p3m_opt_params(interp, argc, argv);
    else {
      Tcl_AppendResult(interp, "expect: inter coulomb <bjerrum>",
		       (char *) NULL);
      return TCL_ERROR;
    }
#else
    return TCL_ERROR;
#endif
  }

  if (coulomb_set_bjerrum(d1) == TCL_ERROR) {
    Tcl_AppendResult(interp, argv[0], "bjerrum length must be positive",
		     (char *) NULL);
    return TCL_ERROR;
  }
    
  argc -= 1;
  argv += 1;

  if (d1 == 0.0 && argc == 0) {
    mpi_bcast_coulomb_params();
    return TCL_OK;
  }

  if(argc < 1) {
    Tcl_AppendResult(interp, "wrong # args for inter coulomb.",
		     (char *) NULL);
    mpi_bcast_coulomb_params();
    return TCL_ERROR;
  }

  /* check method */

#define REGISTER_COULOMB(name, parser)			\
  if(ARG0_IS_S(name))					\
    return parser(interp, argc-1, argv+1);

#ifdef P3M
  REGISTER_COULOMB("p3m", tclcommand_inter_coulomb_parse_p3m);
#endif

  REGISTER_COULOMB("dh", tclcommand_inter_coulomb_parse_dh);    

  if(ARG0_IS_S("rf")) return tclcommand_inter_coulomb_parse_rf(interp, argc-1, argv+1,COULOMB_RF);

  if(ARG0_IS_S("inter_rf")) return tclcommand_inter_coulomb_parse_rf(interp, argc-1, argv+1,COULOMB_INTER_RF);

  REGISTER_COULOMB("mmm1d", tclcommand_inter_coulomb_parse_mmm1d);

  REGISTER_COULOMB("mmm2d", tclcommand_inter_coulomb_parse_mmm2d);

  REGISTER_COULOMB("maggs", tclcommand_inter_coulomb_parse_maggs);

  REGISTER_COULOMB("memd", tclcommand_inter_coulomb_parse_maggs);

  #ifdef MMM1D_GPU
  REGISTER_COULOMB("mmm1dgpu", tclcommand_inter_coulomb_parse_mmm1dgpu);
  #endif

  #ifdef EWALD_GPU
  REGISTER_COULOMB("ewaldgpu", tclcommand_inter_coulomb_parse_ewaldgpu);
  #endif

  #ifdef SCAFACOS
  REGISTER_COULOMB("scafacos", tclcommand_inter_coulomb_parse_scafacos);
  #endif
  
  /* fallback */
  coulomb.method  = COULOMB_NONE;
  coulomb.bjerrum = 0.0;

  mpi_bcast_coulomb_params();

  Tcl_AppendResult(interp, "do not know coulomb method \"",argv[0],
		   "\": coulomb switched off", (char *) NULL);
  
  return TCL_ERROR;
}

/* =========================================================
   ========================================================= */
#endif /*ifdef ELECTROSTATICS */


#ifdef DIPOLES
int tclcommand_inter_parse_magnetic(Tcl_Interp * interp, int argc, char ** argv)
{
  double d1;

  Tcl_ResetResult(interp);

  if(argc == 0) {
    tclprint_to_result_DipolarIA(interp);
    return TCL_OK;
  }
  
  if (! ARG0_IS_D(d1)) {
    Tcl_ResetResult(interp);
    
    if (ARG0_IS_S("mdlc") && ((coulomb.Dmethod == DIPOLAR_DS) || (coulomb.Dmethod == DIPOLAR_MDLC_DS)))
      return tclcommand_inter_magnetic_parse_mdlc_params(interp, argc - 1, argv + 1);

#ifdef DP3M
    if (ARG0_IS_S("mdlc") && ((coulomb.Dmethod == DIPOLAR_P3M) || (coulomb.Dmethod == DIPOLAR_MDLC_P3M)))
      return tclcommand_inter_magnetic_parse_mdlc_params(interp, argc - 1, argv + 1);
    
    if (coulomb.Dmethod == DIPOLAR_P3M)
      return tclcommand_inter_magnetic_parse_dp3m_opt_params(interp, argc, argv);
    else {
      Tcl_AppendResult(interp, "expect: inter magnetic <Dbjerrum>",
		       (char *) NULL);
      return TCL_ERROR;
    }
#else
    return TCL_ERROR;
#endif
  }


  if (dipolar_set_Dbjerrum(d1) == TCL_ERROR) {
    Tcl_AppendResult(interp, argv[0], "Dbjerrum length must be positive",
		     (char *) NULL);
    return TCL_ERROR;
  }
    
  argc -= 1;
  argv += 1;

  if (d1 == 0.0 && argc == 0) {
    mpi_bcast_coulomb_params();
    return TCL_OK;
  }

  if(argc < 1) {
    Tcl_AppendResult(interp, "wrong # args for inter magnetic.",
		     (char *) NULL);
    mpi_bcast_coulomb_params();
    return TCL_ERROR;
  }

  /* check method */

#define REGISTER_DIPOLAR(name, parser)			\
  if(ARG0_IS_S(name))					\
    return parser(interp, argc-1, argv+1);

#ifdef DP3M
  REGISTER_DIPOLAR("p3m", tclcommand_inter_magnetic_parse_dp3m);
#endif

  REGISTER_DIPOLAR("dawaanr", tclcommand_inter_magnetic_parse_dawaanr);

  REGISTER_DIPOLAR("mdds", tclcommand_inter_magnetic_parse_mdds);
  
#ifdef DIPOLAR_DIRECT_SUM
  REGISTER_DIPOLAR("dds-gpu", tclcommand_inter_magnetic_parse_dds_gpu);
#endif

  /* fallback */
  set_dipolar_method_local(DIPOLAR_NONE);
  coulomb.Dbjerrum = 0.0;

  mpi_bcast_coulomb_params();

  Tcl_AppendResult(interp, "do not know magnetic method \"",argv[0],
		   "\": magnetic switched off", (char *) NULL);
  
  return TCL_ERROR;
}
#endif   /* ifdef  DIPOLES */


/********************************************************************************/
/*                                       printing                               */
/********************************************************************************/

int tclprint_to_result_BondedIA(Tcl_Interp *interp, int i)
{
  Bonded_ia_parameters *params = &bonded_ia_params[i];
  char buffer[TCL_INTEGER_SPACE];

  sprintf(buffer, "%d ", i);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  
  switch (params->type) {
  case BONDED_IA_FENE:
    return tclprint_to_result_feneIA(interp, params);
  case BONDED_IA_OIF_LOCAL_FORCES:
	return tclprint_to_result_oiflocalforcesIA(interp, params);
#ifdef MEMBRANE_COLLISION
  case BONDED_IA_OIF_OUT_DIRECTION:
    return tclprint_to_result_oifoutdirectionIA(interp, params);
#endif
#ifdef OIF_GLOBAL_FORCES
  case BONDED_IA_OIF_GLOBAL_FORCES:						
	return tclprint_to_result_oifglobalforcesIA(interp, params);
#endif
      
#ifdef IMMERSED_BOUNDARY
    case BONDED_IA_IBM_TRIEL:
      return tclprint_to_result_ibm_triel(interp, params);
    case BONDED_IA_IBM_VOLUME_CONSERVATION:
      return tclprint_to_result_ibm_volume_conservation(interp, params);
    case BONDED_IA_IBM_TRIBEND:
      return tclprint_to_result_ibm_tribend(interp, params);
      
#endif
      
#ifdef ROTATION
  case BONDED_IA_HARMONIC_DUMBBELL:
    return tclprint_to_result_harmonic_dumbbellIA(interp, params);
#endif
  case BONDED_IA_HARMONIC:
    return tclprint_to_result_harmonicIA(interp, params);
  case BONDED_IA_QUARTIC:
    return tclprint_to_result_quarticIA(interp, params);
#ifdef ELECTROSTATICS
  case BONDED_IA_BONDED_COULOMB:
    return tclprint_to_result_bonded_coulombIA(interp, params);
#endif
#ifdef BOND_ANGLE_OLD
  case BONDED_IA_ANGLE_OLD:
    return tclprint_to_result_angleIA(interp, params);
#endif
#ifdef BOND_ANGLE
  case BONDED_IA_ANGLE_HARMONIC:
    return tclprint_to_result_angle_harmonicIA(interp, params);
  case BONDED_IA_ANGLE_COSINE:
    return tclprint_to_result_angle_cosineIA(interp, params);
  case BONDED_IA_ANGLE_COSSQUARE:
    return tclprint_to_result_angle_cossquareIA(interp, params);
#endif
#ifdef BOND_ANGLEDIST
  case BONDED_IA_ANGLEDIST:
    return tclprint_to_result_angledistIA(interp, params);
#endif
  case BONDED_IA_DIHEDRAL:
    return tclprint_to_result_dihedralIA(interp, params);
#ifdef BOND_ENDANGLEDIST
  case BONDED_IA_ENDANGLEDIST:
    return tclprint_to_result_endangledistIA(interp, params);
#endif
#ifdef TABULATED
  case BONDED_IA_TABULATED:
    return tclprint_to_result_tabulated_bondedIA(interp, params);
#endif
#ifdef OVERLAPPED
  case BONDED_IA_OVERLAPPED:
    return tclprint_to_result_overlapIA(interp, params);
#endif
#ifdef UMBRELLA
  case BONDED_IA_UMBRELLA:
    return tclprint_to_result_umbrellaIA(interp, params);
#endif
#ifdef BOND_CONSTRAINT
  case BONDED_IA_RIGID_BOND:
    return tclprint_to_result_rigid_bond(interp, params);
#endif
#ifdef LENNARD_JONES
  case BONDED_IA_SUBT_LJ:
    return tclprint_to_result_subt_ljIA(interp, params);
#endif
#ifdef BOND_VIRTUAL
  case BONDED_IA_VIRTUAL_BOND:
    Tcl_AppendResult(interp, "VIRTUAL_BOND ", (char *) NULL);
    return (TCL_OK);
#endif
  case BONDED_IA_NONE:
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "unknown bonded interaction number ",buffer,
		     (char *) NULL);
    return (TCL_ERROR);
  default: // keep the compiler happy if interactions are compiled out
    break;
  }
  /* if none of the above */
  Tcl_ResetResult(interp);
  Tcl_AppendResult(interp, "unknown bonded interaction type",(char *) NULL);
  return (TCL_ERROR);
}

int tclprint_to_result_NonbondedIA(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE + 2*TCL_INTEGER_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  if (!data) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "interaction does not exist",
		     (char *) NULL);
    return (TCL_ERROR);
  }

  sprintf(buffer, "%d %d ", i, j);
  Tcl_AppendResult(interp, buffer, (char *) NULL);

#ifdef LENNARD_JONES
  if (data->LJ_cut > 0.0) tclprint_to_result_ljIA(interp,i,j);
#endif

#ifdef LENNARD_JONES_GENERIC
  if (data->LJGEN_cut > 0.0) tclprint_to_result_ljgenIA(interp,i,j);
#endif

#ifdef LJ_ANGLE
  if (data->LJANGLE_cut > 0.0) tclprint_to_result_ljangleIA(interp,i,j);
#endif

#ifdef SMOOTH_STEP
  if (data->SmSt_cut > 0.0) tclprint_to_result_SmStIA(interp,i,j);
#endif

#ifdef HERTZIAN
  if (data->Hertzian_sig > 0.0) tclprint_to_result_HertzianIA(interp,i,j);
#endif

#ifdef GAUSSIAN
  if (data->Gaussian_cut > 0.0) tclprint_to_result_GaussianIA(interp,i,j);
#endif

#ifdef BMHTF_NACL
  if (data->BMHTF_cut > 0.0) tclprint_to_result_BMHTFIA(interp,i,j);
#endif

#ifdef MORSE
  if (data->MORSE_cut > 0.0) tclprint_to_result_morseIA(interp,i,j);
#endif

#ifdef LJCOS
  if (data->LJCOS_cut > 0.0) tclprint_to_result_ljcosIA(interp,i,j);
#endif

#ifdef BUCKINGHAM
  if (data->BUCK_cut > 0.0) tclprint_to_result_buckIA(interp,i,j);
#endif

#ifdef SOFT_SPHERE
  if (data->soft_cut > 0.0) tclprint_to_result_softIA(interp,i,j);
#endif

#ifdef AFFINITY
  if (data->affinity_cut > 0.0) tclprint_to_result_affinityIA(interp,i,j);
#endif
    
#ifdef MEMBRANE_COLLISION
    if (data->membrane_cut > 0.0) tclprint_to_result_membraneIA(interp,i,j);
#endif

#ifdef HAT
  if (data->HAT_r > 0.0) tclprint_to_result_hatIA(interp,i,j);
#endif

#ifdef LJCOS2
  if (data->LJCOS2_cut > 0.0) tclprint_to_result_ljcos2IA(interp,i,j);
#endif

#ifdef COS2
  if (data->COS2_cut > 0.0) tclprint_to_result_cos2IA(interp,i,j);
#endif

#ifdef GAY_BERNE
  if (data->GB_cut > 0.0) tclprint_to_result_gbIA(interp,i,j);
#endif

#ifdef TABULATED
  if (data->TAB_maxval > 0.0)
    Tcl_AppendResult(interp, "tabulated \"", data->TAB_filename,"\"", (char *) NULL);
#endif

#ifdef COMFORCE
  if (data->COMFORCE_flag > 0.0) tclprint_to_result_comforceIA(interp,i,j);
#endif

#ifdef COMFIXED
  if (data->COMFIXED_flag > 0.0) tclprint_to_result_comfixedIA(interp,i,j);
#endif

#ifdef INTER_DPD
  if ((data->dpd_r_cut > 0.0)||(data->dpd_tr_cut > 0.0)) tclprint_to_result_inter_dpdIA(interp,i,j);
#endif

#ifdef INTER_RF
  if (data->rf_on == 1) tclprint_to_result_interrfIA(interp,i,j);
#endif
  
#ifdef MOL_CUT
  if (data->mol_cut_type > 0.0) tclprint_to_result_molcutIA(interp,i,j);
#endif

#ifdef TUNABLE_SLIP
  if (data->TUNABLE_SLIP_r_cut > 0.0) tclprint_to_result_tunable_slipIA(interp,i,j);
#endif
#ifdef SHANCHEN
  if (data->affinity_on == 1 ) tclprint_to_result_affinityIA(interp,i,j);
#endif

  return (TCL_OK);
}

int tclprint_to_result_CoulombIA(Tcl_Interp *interp) 
{
#ifdef ELECTROSTATICS
  char buffer[TCL_DOUBLE_SPACE + 2*TCL_INTEGER_SPACE];
  if (coulomb.method == COULOMB_NONE) {
    Tcl_AppendResult(interp, "coulomb 0.0", (char *) NULL);
    return (TCL_OK);
  }
  Tcl_PrintDouble(interp, coulomb.bjerrum, buffer);
  Tcl_AppendResult(interp, "{coulomb ", buffer, " ", (char *) NULL);
  switch (coulomb.method) {
#ifdef P3M
  case COULOMB_ELC_P3M:
    tclprint_to_result_p3m(interp);
    tclprint_to_result_ELC(interp);
    break;
  case COULOMB_P3M_GPU:
  case COULOMB_P3M: tclprint_to_result_p3m(interp); break;
#endif
  case COULOMB_DH: tclprint_to_result_dh(interp); break;
  case COULOMB_RF: tclprint_to_result_rf(interp,"rf"); break;
  case COULOMB_INTER_RF: tclprint_to_result_rf(interp,"inter_rf"); break;
  case COULOMB_MMM1D: tclprint_to_result_MMM1D(interp); break;
#ifdef MMM1D_GPU
  case COULOMB_MMM1D_GPU: tclprint_to_result_MMM1DGPU(interp); break;
#endif
  case COULOMB_MMM2D: tclprint_to_result_MMM2D(interp); break;
  case COULOMB_MAGGS: tclprint_to_result_Maggs(interp); break;
#ifdef EWALD_GPU
  case COULOMB_EWALD_GPU: tclprint_to_result_ewaldgpu(interp); break;
#endif
#ifdef SCAFACOS
    case COULOMB_SCAFACOS: tclprint_to_result_scafacos(interp); break;
#endif
  default: break;
  }
  Tcl_AppendResult(interp, "}",(char *) NULL);

#else
  Tcl_AppendResult(interp, "ELECTROSTATICS not compiled (see config.hpp)",(char *) NULL);
#endif
  return (TCL_OK);
}

#ifdef DIPOLES
int tclprint_to_result_DipolarIA(Tcl_Interp *interp) 
{
  char buffer[TCL_DOUBLE_SPACE + 2*TCL_INTEGER_SPACE];
  if (coulomb.Dmethod == DIPOLAR_NONE) {
	    Tcl_AppendResult(interp, "magnetic 0.0", (char *) NULL);
    return (TCL_OK);
  }
 
  Tcl_PrintDouble(interp, coulomb.Dbjerrum, buffer);
  Tcl_AppendResult(interp, "{magnetic ", buffer, " ", (char *) NULL);
  switch (coulomb.Dmethod) {
#ifdef DP3M
  case DIPOLAR_MDLC_P3M:
    tclprint_to_result_dp3m(interp);   
    tclprint_to_result_MDLC(interp);
    break;
  case DIPOLAR_P3M: tclprint_to_result_dp3m(interp); break;
#endif
  case DIPOLAR_MDLC_DS:
    tclprint_to_result_Magnetic_dipolar_direct_sum_(interp);
    tclprint_to_result_MDLC(interp);
    break;
  case DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA: tclprint_to_result_DAWAANR(interp); break;
  case DIPOLAR_DS: tclprint_to_result_Magnetic_dipolar_direct_sum_(interp); break;
#ifdef DIPOLAR_DIRECT_SUM
  case DIPOLAR_DS_GPU: tclprint_to_result_dds_gpu(interp); break;
#endif
  default: break;
  }
  Tcl_AppendResult(interp, "}",(char *) NULL);

  return (TCL_OK);
}
#endif

int tclcommand_inter_print_all(Tcl_Interp *interp)
{
  int i, j, start = 1;

  for (i = 0; i < n_bonded_ia; i++) {
    if (bonded_ia_params[i].type != BONDED_IA_NONE) {
      if (start) {
	Tcl_AppendResult(interp, "{", (char *)NULL);
	start = 0;
      }
      else
	Tcl_AppendResult(interp, " {", (char *)NULL);

      if (tclprint_to_result_BondedIA(interp, i) == TCL_ERROR) {
        return TCL_ERROR;
      }
      Tcl_AppendResult(interp, "}", (char *)NULL);
    }
  }
  for (i = 0; i < n_particle_types; i++)
    for (j = i; j < n_particle_types; j++) {
      if (checkIfParticlesInteract(i, j)) {
	if (start) {
	  Tcl_AppendResult(interp, "{", (char *)NULL);
	  start = 0;
	}
	else
	  Tcl_AppendResult(interp, " {", (char *)NULL);
	if (tclprint_to_result_NonbondedIA(interp, i, j) == TCL_ERROR) {
          return TCL_ERROR;
        }
	Tcl_AppendResult(interp, "}", (char *)NULL);
      }
    }
#ifdef ELECTROSTATICS
  if(coulomb.method != COULOMB_NONE) {
    if (start) 
      start = 0;
    else
      Tcl_AppendResult(interp, " ", (char *)NULL);
    /* here the curled braces will be set inside \ref tclprint_to_result_CoulombIA
       because electrostatics might be using several lists */
    tclprint_to_result_CoulombIA(interp);
  }
#endif

#ifdef DIPOLES
  if(coulomb.Dmethod != DIPOLAR_NONE) {
    if (start) 
      start = 0;
    else
      Tcl_AppendResult(interp, " ", (char *)NULL);
    /* here the curled braces will be set inside \ref tclprint_to_result_DipolarIA
       because magnetostatics might be using several lists */
    tclprint_to_result_DipolarIA(interp);
  }
#endif

  if(force_cap != 0.0) {
    char buffer[TCL_DOUBLE_SPACE];
    
    if (start) {
      Tcl_AppendResult(interp, "{", (char *)NULL);
      start = 0;
    }
    else
      Tcl_AppendResult(interp, " {", (char *)NULL);
    if (force_cap == -1.0)
      Tcl_AppendResult(interp, "forcecap individual", (char *)NULL);
    else {
      Tcl_PrintDouble(interp, force_cap, buffer);
      Tcl_AppendResult(interp, "forcecap ", buffer, (char *) NULL);
    }
    Tcl_AppendResult(interp, "}", (char *)NULL);
  }

  return (TCL_OK);
}

int tclcommand_inter_print_bonded(Tcl_Interp *interp, int i)
{
  char buffer[TCL_INTEGER_SPACE];

  Tcl_ResetResult(interp);

  if(i < 0) {
    Tcl_AppendResult(interp, "interaction type must be nonnegative",
		     (char *) NULL);
    return (TCL_ERROR);
  }
  
  /* print specific interaction information */
  if(i<n_bonded_ia) {
    tclprint_to_result_BondedIA(interp, i);
    return TCL_OK;
  }

  sprintf(buffer, "%d", i);
  Tcl_AppendResult(interp, "unknown bonded interaction number ", buffer,
		   (char *) NULL);
  return TCL_ERROR;
}

int tclcommand_inter_print_non_bonded(Tcl_Interp * interp,
			   int part_type_a, int part_type_b)
{
  IA_parameters *data, *data_sym;

  Tcl_ResetResult(interp);
  make_particle_type_exist(part_type_a);
  make_particle_type_exist(part_type_b);
    
  data     = get_ia_param(part_type_a, part_type_b);
  data_sym = get_ia_param(part_type_b, part_type_a);

  if (!data || !data_sym) {
    Tcl_AppendResult(interp, "particle types must be nonnegative",
		     (char *) NULL);
    return TCL_ERROR;
  }

  return tclprint_to_result_NonbondedIA(interp, part_type_a, part_type_b);
}

int tclcommand_inter_parse_non_bonded(Tcl_Interp * interp,
			   int part_type_a, int part_type_b,
			   int argc, char ** argv)
{
  int change;
  
  Tcl_ResetResult(interp);

  if (argc <= 0) {
    Tcl_AppendResult(interp, "wrong # args:  should be \"",
		     "inter <type 1> <type 2> ?interaction? ?values?\"",
		     (char *) NULL);
    return TCL_ERROR;
  }

  /* get interaction parameters */

  while (argc > 0) {
    /* The various parsers return the number of parsed parameters.
       If an error occured, 0 should be returned, since none of the parameters were
       understood */

    /* that's just for the else below... */
    if (0);

#define REGISTER_NONBONDED(name, parser)				\
    else if (ARG0_IS_S(name))						\
      change = parser(interp, part_type_a, part_type_b, argc, argv)

#ifdef LENNARD_JONES
    REGISTER_NONBONDED("lennard-jones", tclcommand_inter_parse_lj);
#endif

#ifdef LENNARD_JONES_GENERIC
    REGISTER_NONBONDED("lj-gen", tclcommand_inter_parse_ljgen);
#endif

#ifdef LJ_ANGLE
    REGISTER_NONBONDED("lj-angle", tclcommand_inter_parse_ljangle);
#endif

#ifdef SMOOTH_STEP
    REGISTER_NONBONDED("smooth-step", tclcommand_inter_parse_SmSt);
#endif

#ifdef HERTZIAN
    REGISTER_NONBONDED("hertzian", tclcommand_inter_parse_hertzian);
#endif

#ifdef GAUSSIAN
    REGISTER_NONBONDED("gaussian", tclcommand_inter_parse_gaussian);
#endif

#ifdef BMHTF_NACL
    REGISTER_NONBONDED("bmhtf-nacl", tclcommand_inter_parse_BMHTF);
#endif

#ifdef MORSE
    REGISTER_NONBONDED("morse", tclcommand_inter_parse_morse);
#endif

#ifdef LJCOS
    REGISTER_NONBONDED("lj-cos", tclcommand_inter_parse_ljcos);
#endif

#ifdef BUCKINGHAM
    REGISTER_NONBONDED("buckingham", tclcommand_inter_parse_buckingham);
#endif

#ifdef SOFT_SPHERE
    REGISTER_NONBONDED("soft-sphere", tclcommand_inter_parse_soft);
#endif

#ifdef AFFINITY
    REGISTER_NONBONDED("affinity", tclcommand_inter_parse_affinity);
#endif
      
#ifdef MEMBRANE_COLLISION
    REGISTER_NONBONDED("membrane", tclcommand_inter_parse_membrane);
#endif

#ifdef HAT
    REGISTER_NONBONDED("hat", tclcommand_inter_parse_hat);
#endif

#ifdef COMFORCE
    REGISTER_NONBONDED("comforce", tclcommand_inter_parse_comforce);
#endif

#ifdef LJCOS2
    REGISTER_NONBONDED("lj-cos2", tclcommand_inter_parse_ljcos2);
#endif

#ifdef COS2
    REGISTER_NONBONDED("cos2", tclcommand_inter_parse_cos2);
#endif


#ifdef COMFIXED
    REGISTER_NONBONDED("comfixed", tclcommand_inter_parse_comfixed);
#endif

#ifdef GAY_BERNE
    REGISTER_NONBONDED("gay-berne", tclcommand_inter_parse_gb);
#endif

#ifdef TABULATED
    REGISTER_NONBONDED("tabulated", tclcommand_inter_parse_tab);
#endif
#ifdef INTER_DPD
    REGISTER_NONBONDED("inter_dpd", tclcommand_inter_parse_inter_dpd);
#endif
#ifdef INTER_RF
    REGISTER_NONBONDED("inter_rf", tclcommand_inter_parse_interrf);
#endif
#ifdef TUNABLE_SLIP
    REGISTER_NONBONDED("tunable_slip", tclcommand_inter_parse_tunable_slip);
#endif
#ifdef MOL_CUT
    REGISTER_NONBONDED("molcut", tclcommand_inter_parse_molcut);
#endif
  
#ifdef SHANCHEN
    REGISTER_NONBONDED("affinity",tclcommand_inter_parse_affinity);
#endif
 
    else {
      Tcl_AppendResult(interp, "excessive parameter/unknown interaction type \"", argv[0],
		       "\" in parsing non bonded interaction",
		       (char *) NULL);
      return TCL_ERROR;
    }

    if (change <= 0)
      return TCL_ERROR;

    argc -= change;
    argv += change;
  }
  return TCL_OK;
}

int tclcommand_inter_print_partner_num(Tcl_Interp *interp, int bond_type)
{
  Bonded_ia_parameters * params;
  char buffer[TCL_INTEGER_SPACE];

  if(bond_type < 0) {
    Tcl_AppendResult(interp, "interaction type must be nonnegative",
		     (char *) NULL);
    return TCL_ERROR;
  }

  if(bond_type < n_bonded_ia) {
    params = &bonded_ia_params[bond_type];
    sprintf(buffer, "%d", params->num);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    return TCL_OK;
  }
 
  sprintf(buffer, "%d", bond_type);
  Tcl_AppendResult(interp, "unknown bonded interaction number ", buffer,
		   (char *) NULL);
  return TCL_ERROR;
}

#ifdef ELECTROSTATICS
/* print applied/induced efields for capacitor feature in ic-mmm2d and ic-elc */
int tclcommand_print_efield_capacitors(ClientData data, Tcl_Interp * interp, int argc, char ** argv)
{
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE + 2];
  double value=0;

  if (!(
        (coulomb.method == COULOMB_MMM2D   && mmm2d_params.const_pot_on)
#ifdef P3M
     || (coulomb.method == COULOMB_ELC_P3M && elc_params.const_pot_on)
#endif
     )) {
    Tcl_AppendResult(interp, "Electric field output only available for mmm2d or p3m+elc with capacitor feature", (char *) NULL);
    return TCL_ERROR;
  }
  else if (argc > 2) {
    Tcl_AppendResult(interp, "wrong # arguments: efield_caps <{total} | {induced} | {applied}>", (char *) NULL);
    return TCL_ERROR;
  }
  else if (argc == 1 || ARG1_IS_S("total")) {
	  value = field_induced + field_applied;
  }
  else if (ARG1_IS_S("induced")) {
	  value = field_induced;
  }
  else if (ARG1_IS_S("applied")) {
	  value = field_applied;
  }

  Tcl_PrintDouble(interp, value, buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);

  return TCL_OK;
}
#endif

/********************************************************************************/
/*                                       parsing                                */
/********************************************************************************/

#ifdef BOND_VIRTUAL
int tclcommand_inter_parse_virtual_bonds(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  CHECK_VALUE(virtual_set_params(bond_type), "bond type must be nonnegative");
}
#endif

int tclcommand_inter_parse_bonded(Tcl_Interp *interp,
		       int bond_type,
		       int argc, char ** argv)
{
  if (ARG0_IS_S("num")) {
    if (argc == 1)
      return tclcommand_inter_print_partner_num(interp, bond_type);
    else {
	Tcl_AppendResult(interp, "too many parameters",
			 (char *) NULL);
	return TCL_ERROR;
    }
  }

#define REGISTER_BONDED(name,parser)			\
  if (ARG0_IS_S(name)) return parser(interp, bond_type, argc, argv);
  
  REGISTER_BONDED("fene", tclcommand_inter_parse_fene);
  REGISTER_BONDED("oif_local_forces", tclcommand_inter_parse_oif_local_forces);
#ifdef OIF_GLOBAL_FORCES
  REGISTER_BONDED("oif_global_forces", tclcommand_inter_parse_oif_global_forces);
#endif
#ifdef MEMBRANE_COLLISION
  REGISTER_BONDED("oif_out_direction", tclcommand_inter_parse_oif_out_direction);
#endif
  // IMMERSED_BOUNDARY
#ifdef IMMERSED_BOUNDARY
  REGISTER_BONDED("ibm_triel", tclcommand_inter_parse_ibm_triel);
  REGISTER_BONDED("ibm_volcons", tclcommand_inter_parse_ibm_volume_conservation);
  REGISTER_BONDED("ibm_tribend", tclcommand_inter_parse_ibm_tribend);
#endif
  
  REGISTER_BONDED("harmonic", tclcommand_inter_parse_harmonic);
#ifdef ROTATION
  REGISTER_BONDED("harmonic_dumbbell", tclcommand_inter_parse_harmonic_dumbbell);
#endif
  REGISTER_BONDED("quartic", tclcommand_inter_parse_quartic);
#ifdef ELECTROSTATICS
  REGISTER_BONDED("bonded_coulomb", tclcommand_inter_parse_bonded_coulomb);  
#endif
#ifdef LENNARD_JONES  
  REGISTER_BONDED("subt_lj", tclcommand_inter_parse_subt_lj);
#endif
#ifdef BOND_ANGLE_OLD
  REGISTER_BONDED("angle", tclcommand_inter_parse_angle);
#endif
#ifdef BOND_ANGLE
  REGISTER_BONDED("angle_harmonic", tclcommand_inter_parse_angle_harmonic);
  REGISTER_BONDED("angle_cosine", tclcommand_inter_parse_angle_cosine);
  REGISTER_BONDED("angle_cossquare", tclcommand_inter_parse_angle_cossquare);
#endif
#ifdef BOND_ANGLEDIST
  REGISTER_BONDED("angledist", tclcommand_inter_parse_angledist);
#endif
  REGISTER_BONDED("dihedral", tclcommand_inter_parse_dihedral);
#ifdef BOND_ENDANGLEDIST
  REGISTER_BONDED("endangledist", tclcommand_inter_parse_endangledist);
#endif
#ifdef TABULATED
  REGISTER_BONDED("tabulated", tclcommand_inter_parse_tabulated_bonded);
#endif
#ifdef OVERLAPPED
  REGISTER_BONDED("overlapped", tclcommand_inter_parse_overlapped_bonded);
#endif
#ifdef UMBRELLA
  REGISTER_BONDED("umbrella", tclcommand_inter_parse_umbrella);
#endif
#ifdef BOND_CONSTRAINT
  REGISTER_BONDED("rigid_bond", tclcommand_inter_parse_rigid_bond);
#endif
#ifdef BOND_VIRTUAL
  REGISTER_BONDED("virtual_bond", tclcommand_inter_parse_virtual_bonds);
#endif
#ifdef HYDROGEN_BOND
  REGISTER_BONDED("hydrogen_bond", tclcommand_inter_parse_hydrogen_bond);
#endif
#ifdef TWIST_STACK
  REGISTER_BONDED("twist_stack", tclcommand_inter_parse_twist_stack);  
#endif
  Tcl_AppendResult(interp, "unknown bonded interaction type \"", argv[0],
		   "\"", (char *) NULL);
  return TCL_ERROR;
}

int tclcommand_inter_parse_rest(Tcl_Interp * interp, int argc, char ** argv)
{
  if(ARG0_IS_S("forcecap"))
    return tclcommand_inter_parse_forcecap(interp,argc-1, argv+1);

#if defined(LENNARD_JONES) || defined(LENNARD_JONES_GENERIC)
  if(ARG0_IS_S("ljforcecap"))
    return tclcommand_inter_parse_ljforcecap(interp, argc-1, argv+1);
#endif

#ifdef LJ_ANGLE
  if(ARG0_IS_S("ljangleforcecap"))
    return tclcommand_inter_parse_ljangleforcecap(interp, argc-1, argv+1);
#endif

  
#ifdef MORSE
  if(ARG0_IS_S("morseforcecap"))
    return tclcommand_inter_parse_morseforcecap(interp, argc-1, argv+1);
#endif

#ifdef BUCKINGHAM
  if(ARG0_IS_S("buckforcecap"))
    return tclcommand_inter_parse_buckforcecap(interp, argc-1, argv+1);
#endif

#ifdef TABULATED
  if(ARG0_IS_S("tabforcecap"))
    return tclcommand_inter_parse_tabforcecap(interp, argc-1, argv+1);
#endif

  if(ARG0_IS_S("coulomb")) {
    #ifdef ELECTROSTATICS
      return tclcommand_inter_parse_coulomb(interp, argc-1, argv+1);
   #else
       Tcl_AppendResult(interp, "ELECTROSTATICS not compiled (see config.hpp)", (char *) NULL);
    #endif
  }
  
  if(ARG0_IS_S("magnetic")) {
   #ifdef DIPOLES
      return tclcommand_inter_parse_magnetic(interp, argc-1, argv+1);
    #else
      Tcl_AppendResult(interp, "DIPOLES not compiled (see config.hpp)", (char *) NULL);
    #endif
  }
  
  Tcl_AppendResult(interp, "unknown interaction type \"", argv[0],
		   "\"", (char *) NULL);

  return TCL_ERROR;
}

int tclcommand_inter(ClientData _data, Tcl_Interp *interp,
	  int argc, char **argv)
{
  int i, j, err_code = TCL_OK, is_i1, is_i2;

  Tcl_ResetResult(interp);

  /* first we handle the special cases

     1. no parameters
     2. one parameter
     3. two parameters

     then the bonded interactions
     then the non bonded interactions
     then the rest
   */

  if (argc == 1) {
    /* no argument -> print all interaction informations. */
    err_code = tclcommand_inter_print_all(interp);
  }
  else if (argc == 2) {
    /* There is only 1 parameter, bonded ia printing or force caps */

    if (ARG1_IS_I(i))
      err_code = tclcommand_inter_print_bonded(interp, i);
    else {
      Tcl_ResetResult(interp);
      err_code = tclcommand_inter_parse_rest(interp, argc-1, argv+1);
    }
  }
  else if (argc == 3) {
    /* There are only 2 parameters, non_bonded printing */
    
    is_i1 = ARG_IS_I(1, i);
    is_i2 = ARG_IS_I(2, j);

    Tcl_ResetResult(interp);

    if (is_i1 && is_i2)
      err_code = tclcommand_inter_print_non_bonded(interp, i, j);
    else if (is_i1)
      err_code = tclcommand_inter_parse_bonded(interp, i, argc-2, argv+2);
    else
      err_code = tclcommand_inter_parse_rest(interp, argc-1, argv+1);
  }
  else {
    /****************************************************
     * Here we have more than 2 parameters
     ****************************************************/

    is_i1 = ARG_IS_I(1, i);
    is_i2 = ARG_IS_I(2, j);

    Tcl_ResetResult(interp);

    // non bonded interactions
    if (is_i1 && is_i2)
      err_code = tclcommand_inter_parse_non_bonded(interp, i, j, argc-3, argv+3);
    else if (is_i1)
      // bonded interactions
      err_code = tclcommand_inter_parse_bonded(interp, i, argc-2, argv+2);
    else
      // named interactions
      err_code = tclcommand_inter_parse_rest(interp, argc-1, argv+1);
  }
  /* check for background errors which have not been handled so far */
  return gather_runtime_errors(interp, err_code);
}

int tclcallback_min_global_cut(Tcl_Interp *interp, void *_data)
{
  min_global_cut = *((double *)_data);
  mpi_bcast_parameter(FIELD_MIN_GLOBAL_CUT);
  return TCL_OK;
}
