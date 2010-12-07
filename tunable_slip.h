#include "utils.h"
#include "parser.h"

/** \file tunable_slip.h
 *  Routines to generate tunable-slip boundary conditions.
 *  J.Smiatek, M.P. Allen, F. Schmid:
 *  "Tunable-slip boundaries for coarse-grained simulations of fluid flow", Europ. Phys. J. E 26, 115 (2008) 
*/


#ifdef TUNABLE_SLIP

MDINLINE int tunable_slip_set_params(int part_type_a, int part_type_b,
		      double temp, double gamma, double r_cut,
		      double time, double vx, double vy, double vz)
{
  IA_parameters *data, *data_sym;

  make_particle_type_exist(part_type_a);
  make_particle_type_exist(part_type_b);
    
  data     = get_ia_param(part_type_a, part_type_b);
  data_sym = get_ia_param(part_type_b, part_type_a);
  
  if (!data || !data_sym) {
      return TCL_ERROR;
  }
  
  /* TUNABLE SLIP should be symmetrical! */
  data_sym->TUNABLE_SLIP_temp  = data->TUNABLE_SLIP_temp     = temp;
  data_sym->TUNABLE_SLIP_gamma = data->TUNABLE_SLIP_gamma    = gamma;
  data_sym->TUNABLE_SLIP_r_cut = data->TUNABLE_SLIP_r_cut    = r_cut;
  data_sym->TUNABLE_SLIP_time  = data->TUNABLE_SLIP_time     = time;
  data_sym->TUNABLE_SLIP_vx  = data->TUNABLE_SLIP_vx     = vx;
  data_sym->TUNABLE_SLIP_vy  = data->TUNABLE_SLIP_vy     = vy;
  data_sym->TUNABLE_SLIP_vz  = data->TUNABLE_SLIP_vz     = vz;
 

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);
  mpi_bcast_ia_params(part_type_b, part_type_a);
  
  return TCL_OK;
}

MDINLINE int tclprint_to_result_tunable_slipIA(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  Tcl_PrintDouble(interp, data->TUNABLE_SLIP_temp, buffer);
  Tcl_AppendResult(interp, "tunable_slip ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->TUNABLE_SLIP_gamma, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->TUNABLE_SLIP_r_cut, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->TUNABLE_SLIP_time, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->TUNABLE_SLIP_vx, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->TUNABLE_SLIP_vy, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->TUNABLE_SLIP_vz, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  return TCL_OK;
}

MDINLINE int tclcommand_inter_parse_tunable_slip(Tcl_Interp * interp,
			    int part_type_a, int part_type_b,
			    int argc, char ** argv)
{
  double temp, gamma, r_cut, time, vx, vy, vz;
  int change;

  if (argc < 8) {
    Tcl_AppendResult(interp, "tunable_slip needs 7 parameters: "
		     "<tunable_slip_temp> <tunable_slip_gamma> <tunable_slip_r_cut> <tunable_slip_time> <vx> <vy> <vz>",
		     (char *) NULL);
    return 0;
  }

  /* copy lj-cos parameters */
  if ((! ARG_IS_D(1, temp))    ||
      (! ARG_IS_D(2, gamma))   ||
      (! ARG_IS_D(3, r_cut))   ||
      (! ARG_IS_D(4, time))   ||
      (! ARG_IS_D(5, vx))   ||
      (! ARG_IS_D(6, vy))   ||
      (! ARG_IS_D(7, vz)    )) {
    Tcl_AppendResult(interp, "tunable_slip needs 7 DOUBLE parameters: "
		     "<tunable_slip_temp> <tunable_slip_gamma> <tunable_slip_r_cut> <tunable_slip_time> <vx> <vy> <vz>",
		     (char *) NULL);
    return 0;
  }
  change = 8;

  if (tunable_slip_set_params(part_type_a, part_type_b, temp, gamma, r_cut, time, vx, vy, vz) == TCL_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  }

  return change;
}

MDINLINE void add_tunable_slip_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params, double d[3], double dist, double force[3])
{
  double gamma_t=0.0;
  double gamma_t_sqrt=0.0;
  double pre_diss=0.0;
  double pre_rand=0.0;
  double scale_vx=0.0;
  double scale_vy=0.0;
  double scale_vz=0.0; 

  if(dist < ia_params->TUNABLE_SLIP_r_cut) {
    gamma_t += 1-(dist/(ia_params->TUNABLE_SLIP_r_cut));
   }else{
    gamma_t += 0.0;
   }
  gamma_t_sqrt += sqrt(gamma_t);
  pre_diss += ((ia_params->TUNABLE_SLIP_gamma)/ia_params->TUNABLE_SLIP_time);
  pre_rand += sqrt((24.0*ia_params->TUNABLE_SLIP_temp*ia_params->TUNABLE_SLIP_gamma)/ia_params->TUNABLE_SLIP_time);

  /* Rescaling of reference velocity */ 
  scale_vx += ia_params->TUNABLE_SLIP_vx*ia_params->TUNABLE_SLIP_time;
  scale_vy += ia_params->TUNABLE_SLIP_vy*ia_params->TUNABLE_SLIP_time;
  scale_vz += ia_params->TUNABLE_SLIP_vz*ia_params->TUNABLE_SLIP_time;

  force[0] += -(pre_diss*gamma_t*(p1->m.v[0]-scale_vx))+(pre_rand*gamma_t_sqrt*(d_random()-0.5));
  force[1] += -(pre_diss*gamma_t*(p1->m.v[1]-scale_vy))+(pre_rand*gamma_t_sqrt*(d_random()-0.5));
  force[2] += -(pre_diss*gamma_t*(p1->m.v[2]-scale_vz))+(pre_rand*gamma_t_sqrt*(d_random()-0.5));
    
  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: TUNABLE_SLIP   f = (%.3e,%.3e,%.3e) with part id=%d\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: TUNABLE_SLIP   f = (%.3e,%.3e,%.3e) with part id=%d\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity));
 
}

#endif

