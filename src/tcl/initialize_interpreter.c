
#include <tcl.h>

#include "adresso_tcl.h"
#include "bin_tcl.h"
#include "binary_file_tcl.h"
#include "blockfile_tcl.h"
#include "cells_tcl.h"
#include "constraint_tcl.h"
#include "domain_decomposition_tcl.h"
#include "dpd_tcl.h"
#include "global_tcl.h"
#include "grid_tcl.h"
#include "iccp3m_tcl.h"
#include "imd_tcl.h"
#include "integrate_tcl.h"
#include "interaction_data_tcl.h"
#include "lb-boundaries_tcl.h"
#include "lb_tcl.h"
#include "lj_tcl.h"
#include "maggs_tcl.h"
#include "metadynamics_tcl.h"
#include "nemd_tcl.h"
#include "p3m-dipolar_tcl.h"
#include "p3m_tcl.h"
#include "particle_data_tcl.h"
#include "polymer_tcl.h"
#include "pressure_tcl.h"
#include "random_tcl.h"
#include "statistics_chain_tcl.h"
#include "statistics_cluster_tcl.h"
#include "statistics_fluid_tcl.h"
#include "statistics_tcl.h"
#include "thermostat_tcl.h"
#include "topology_tcl.h"
#include "uwerr_tcl.h"
#include "virtual_sites_com_tcl.h"
#include "cuda_init_tcl.h"
#include "tcl/cuda_init_tcl.h"


#define REGISTER_COMMAND(name, routine)					\
  Tcl_CreateCommand(interp, name, (Tcl_CmdProc *)routine, 0, NULL);

void register_tcl_commands(Tcl_Interp* interp) {
  /* in cells.c */
  REGISTER_COMMAND("cellsystem", tclcommand_cellsystem);
  /* in integrate.c */
  REGISTER_COMMAND("invalidate_system", tclcommand_invalidate_system);
  REGISTER_COMMAND("integrate", tclcommand_integrate);
  /* in global.c */
  REGISTER_COMMAND("setmd", tclcommand_setmd);
  /* in grid.c */
  REGISTER_COMMAND("change_volume", tclcommand_change_volume);
  /* in global.c */
  REGISTER_COMMAND("code_info", tclcommand_code_info);
  /* in interaction_data.c */
  REGISTER_COMMAND("inter",tclcommand_inter);
  /* in particle_data.c */
  REGISTER_COMMAND("part",tclcommand_part);
  /* in file binaryfile.c */
  REGISTER_COMMAND("writemd", tclcommand_writemd);
  REGISTER_COMMAND("readmd", tclcommand_readmd);
  /* in file statistics.c */
  REGISTER_COMMAND("analyze", tclcommand_analyze);
  /* in file polymer.c */
  REGISTER_COMMAND("polymer", tclcommand_polymer);
  REGISTER_COMMAND("counterions", tclcommand_counterions);
  REGISTER_COMMAND("salt", tclcommand_salt);
  REGISTER_COMMAND("velocities", tclcommand_velocities);
  REGISTER_COMMAND("maxwell_velocities", tclcommand_maxwell_velocities);
  REGISTER_COMMAND("crosslink", tclcommand_crosslink);
  REGISTER_COMMAND("diamond", tclcommand_diamond);
  REGISTER_COMMAND("icosaeder", tclcommand_icosaeder);
  /* in file imd.c */
  REGISTER_COMMAND("imd", tclcommand_imd);
  /* in file random.c */
  REGISTER_COMMAND("t_random", tclcommand_t_random);
  REGISTER_COMMAND("bit_random", tclcommand_bit_random);
  /* in file blockfile_tcl.c */
  REGISTER_COMMAND("blockfile", tclcommand_blockfile);
  /* in constraint.c */
  REGISTER_COMMAND("constraint", tclcommand_constraint);
  /* in uwerr.c */
  REGISTER_COMMAND("uwerr", tclcommand_uwerr);
  /* in nemd.c */
  REGISTER_COMMAND("nemd", tclcommand_nemd);
  /* in thermostat.c */
  REGISTER_COMMAND("thermostat", tclcommand_thermostat);
  /* in bin.c */
  REGISTER_COMMAND("bin", tclcommand_bin);
  /* in lb.c */

  REGISTER_COMMAND("lbfluid", tclcommand_lbfluid);
  REGISTER_COMMAND("lbnode", tclcommand_lbnode);
  REGISTER_COMMAND("lbboundary", tclcommand_lbboundary);
  /* in utils.h */
  REGISTER_COMMAND("replacestdchannel", tclcommand_replacestdchannel);
  /* in iccp3m.h */
#ifdef ELECTROSTATICS
#ifdef P3M
  REGISTER_COMMAND("iccp3m", tclcommand_iccp3m);
#endif 
#endif 
  /* in adresso.h */
  REGISTER_COMMAND("adress", tclcommand_adress);
#ifdef ADRESS
  /* #ifdef THERMODYNAMIC_FORCE */
  REGISTER_COMMAND("thermodynamic_force", tclcommand_thermodynamic_force);
  /* #endif */
  REGISTER_COMMAND("update_adress_weights", tclcommand_update_adress_weights);
#endif
#ifdef METADYNAMICS
  /* in metadynamics.c */
  REGISTER_COMMAND("metadynamics", tclcommand_metadynamics);
#endif
#ifdef LB_GPU
  /* in lbgpu_cfile.c */
  REGISTER_COMMAND("lbnode_extforce", tclcommand_lbnode_extforce_gpu);

  //REGISTER_COMMAND("lbprint", tclcommand_lbprint_gpu);
#endif
#ifdef CUDA
  REGISTER_COMMAND("cuda", tclcommand_cuda);
#endif

}
