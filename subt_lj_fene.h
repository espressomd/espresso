#ifndef SUBT_LJ_FENE_H
#define SUBT_LJ_FENE_H
/** \file subt_lj_fene.h
 *  Routines to subtract the LENNARD-JONES Energy and/or the LENNARD-JONES force 
 *  from the FENE Energy and/or the force for a particle pair.
 *  \ref forces.c
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:sathish@mpip-mainz.mpg.de">sathish</a>
*/

#ifdef LENNARD_JONES

/************************************************************/

/// parameters for the subtract LJ from a fene bond potential. Deprecated.
MDINLINE int subt_lj_fene_set_params(int bond_type, double k, double r)
{
  if(bond_type < 0)
    return TCL_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].p.subt_lj_fene.k = k;
  bonded_ia_params[bond_type].p.subt_lj_fene.r = r;
  bonded_ia_params[bond_type].type = BONDED_IA_SUBT_LJ_FENE;
  bonded_ia_params[bond_type].p.subt_lj_fene.r2 = SQR(bonded_ia_params[bond_type].p.subt_lj_fene.r);
  bonded_ia_params[bond_type].num  = 1;

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(bond_type, -1); 

  return TCL_OK;
}

/** Computes the difference between the FENE and the LENNARD-JONES pair forces 
    and adds this force to the particle forces (see \ref #inter). 
    @param p1        Pointer to first particle.
    @param p2        Pointer to second/middle particle.
    @param type_num  bond type number of this interaction (see \ref #inter).
*/
MDINLINE int calc_subt_lj_fene_pair_force(Particle *p1, Particle *p2, Bonded_ia_parameters *iaparams, double dx[3], double force[3])
{
  int i;
  double fac_fene=0.0, fac_lj=0.0, fac;
  IA_parameters *ia_params;
  double r_off, frac2, frac6;
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);

  if(dist >= iaparams->p.subt_lj_fene.r)
    return 1;

  fac_fene = iaparams->p.subt_lj_fene.k / (1.0 - dist2/iaparams->p.subt_lj_fene.r2);
  
  FENE_TRACE(if(fac_fene > 50) fprintf(stderr,"WARNING: FENE force factor between Pair (%d,%d) large: %f at distance %f\n", p1->p.identity,p2->p.identity,fac_fene,dist) );

  ia_params = get_ia_param(p1->p.type,p2->p.type);
  if(dist < ia_params->LJ_cut+ia_params->LJ_offset) { 
    r_off = dist - ia_params->LJ_offset;

    if(r_off > ia_params->LJ_capradius) {
      /* normal case: resulting force/energy smaller than capping. */
      frac2 = SQR(ia_params->LJ_sig/r_off);
      frac6 = frac2*frac2*frac2;
      fac_lj   = 48.0 * ia_params->LJ_eps * frac6*(frac6 - 0.5) / (r_off * dist);
    }
    else if(dist > 0.0) {
      /* capped part of lj potential. */
      frac2 = SQR(ia_params->LJ_sig/ia_params->LJ_capradius);
      frac6 = frac2*frac2*frac2;
      fac_lj   = 48.0 * ia_params->LJ_eps * frac6*(frac6 - 0.5) / (ia_params->LJ_capradius * dist);
    }
  }
  
  fac = -(fac_fene + fac_lj);
          
  for(i=0;i<3;i++)
    force[i] = fac*dx[i];
  
  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: SUBT_LJ_FENE f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,fac));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: SUBT_LJ_FENE f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,fac));

  return 0;
}

MDINLINE int subt_lj_fene_pair_energy(Particle *p1, Particle *p2, Bonded_ia_parameters *iaparams, double dx[3], double *_energy)
{
  double energy_fene = 0.0, energy_lj = 0.0;
  IA_parameters *ia_params;
  double r_off, frac2, frac6;
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);

  if(dist >= iaparams->p.subt_lj_fene.r)
    return 1;
  
  energy_fene = -0.5*iaparams->p.subt_lj_fene.k*iaparams->p.subt_lj_fene.r2;
  energy_fene *= log((1.0 - dist2/iaparams->p.subt_lj_fene.r2));
  
  ia_params = get_ia_param(p1->p.type,p2->p.type);
  if(dist < ia_params->LJ_cut+ia_params->LJ_offset) {
    r_off = dist - ia_params->LJ_offset;
    if(r_off > ia_params->LJ_capradius) {
      /* normal case: resulting force/energy smaller than capping. */
      frac2 = SQR(ia_params->LJ_sig/r_off);
      frac6 = frac2*frac2*frac2;
      energy_lj = 4.0*ia_params->LJ_eps*(SQR(frac6)-frac6+ia_params->LJ_shift);
    }   
    else if(dist > 0.0) {
      /* capped part of lj potential. */
      frac2 = SQR(ia_params->LJ_sig/ia_params->LJ_capradius);
      frac6 = frac2*frac2*frac2;
      energy_lj = 4.0*ia_params->LJ_eps*(SQR(frac6)-frac6+ia_params->LJ_shift);
    }
  }
  *_energy = energy_fene-energy_lj;
  return 0;
}

#endif
#endif
