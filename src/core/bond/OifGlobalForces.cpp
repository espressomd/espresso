#include"OifGlobalForces.hpp"
#include "grid.hpp" // get_mi_vector //unfold_position
#include "interaction_data.hpp"

//---Static Variables---
//initialize the static variables for bonds
double Bond::OifGlobalForces::m_area = 0.0;
double Bond::OifGlobalForces::m_VOL_volume = 0.0;

//---Functions---
int Bond::OifGlobalForces::calc_oif_global(Particle* p1, int bl_id, double* partArea, 
					    double* VOL_partVol)

{
  // z volume
  double VOL_A, VOL_norm[3], VOL_dn, VOL_hz;
  Particle *p2, *p3;
  double p11[3], p22[3], p33[3];
  int img[3];
  double AA[3], BB[3];

  if(auto bond_partners =  get_n_bond_partners<2>(p1, bl_id)){

    p2 = (*bond_partners)[0];
    p3 = (*bond_partners)[1];

    // getting unfolded positions of all particles
#ifdef GHOST_FLAG
    // first find out which particle out of p1, p2 (possibly p3, p4) is not
    // a ghost particle. In almost all cases it is p1, however, it might be
    // other one. we call this particle reference particle.
    if (p1->l.ghost != 1) {
      // unfold non-ghost particle using image, because for physical
      // particles, the structure p->l.i is correctly set
      memmove(p11, p1->r.p, 3 * sizeof(double));
      memmove(img, p1->l.i, 3 * sizeof(int));
      unfold_position(p11, img);
      // other coordinates are obtained from its relative positions to the
      // reference particle
      get_mi_vector(AA, p2->r.p, p11);
      get_mi_vector(BB, p3->r.p, p11);
      for (int i = 0; i < 3; i++) {
	p22[i] = p11[i] + AA[i];
	p33[i] = p11[i] + BB[i];
      }
    } else {
      // in case the first particle is a ghost particle
      if (p2->l.ghost != 1) {
	memmove(p22, p2->r.p, 3 * sizeof(double));
	memmove(img, p2->l.i, 3 * sizeof(int));
	unfold_position(p22, img);
	get_mi_vector(AA, p1->r.p, p22);
	get_mi_vector(BB, p3->r.p, p22);
	for (int i = 0; i < 3; i++) {
	  p11[i] = p22[i] + AA[i];
	  p33[i] = p22[i] + BB[i];
	}
      } else {
	// in case the first and the second particle are ghost particles
	if (p3->l.ghost != 1) {

	  memmove(p33, p3->r.p, 3 * sizeof(double));
	  memmove(img, p3->l.i, 3 * sizeof(int));
	  unfold_position(p33, img);
	  get_mi_vector(AA, p1->r.p, p33);
	  get_mi_vector(BB, p2->r.p, p33);
	  for (int i = 0; i < 3; i++) {
	    p11[i] = p33[i] + AA[i];
	    p22[i] = p33[i] + BB[i];
	  }
	} else {
	  printf("Something wrong in oif_global_forces.hpp: All particles "
		 "in a bond are ghost particles, impossible to unfold the "
		 "positions...");
	  return 2;
	}
      }
    }
#endif
#ifndef GHOST_FLAG
    // if ghost flag was not defined we have no other option than to assume
    // the first particle is a physical one.
    memmove(p11, p1->r.p, 3 * sizeof(double));
    memmove(img, p1->l.i, 3 * sizeof(int));
    unfold_position(p11, img);
    // other coordinates are obtained from its relative positions to the
    // reference particle
    get_mi_vector(AA, p2->r.p, p11);
    get_mi_vector(BB, p3->r.p, p11);
    for (int i = 0; i < 3; i++) {
      p22[i] = p11[i] + AA[i];
      p33[i] = p11[i] + BB[i];
    }
#endif
    // unfolded positions correct
    VOL_A = area_triangle(p11, p22, p33);
    *partArea += VOL_A;

    get_n_triangle(p11, p22, p33, VOL_norm);
    VOL_dn = normr(VOL_norm);
    VOL_hz = 1.0 / 3.0 * (p11[2] + p22[2] + p33[2]);
    *VOL_partVol += VOL_A * -1 * VOL_norm[2] / VOL_dn * VOL_hz;

  }
  else{
    return 2;
  };

  return 1;
}

int Bond::OifGlobalForces::add_bonded_energy(Particle *p1, int bl_id)
{

  return 0;

}


int Bond::OifGlobalForces::add_bonded_force(Particle *p1, int bl_id)
{

  double VOL_force[3];
  double VOL_A, VOL_norm[3], VOL_dn, VOL_vv;

  double aa, force1[3], force2[3], force3[3], rh[3], hn, h[3];
  int k;

  Particle *p2, *p3;
  double p11[3], p22[3], p33[3];
  double AA[3], BB[3];
  int img[3];



  if(auto bond_partners =  get_n_bond_partners<2>(p1, bl_id)){
    p2 = (*bond_partners)[0];
    p3 = (*bond_partners)[1];


#ifdef GHOST_FLAG
    // first find out which particle out of p1, p2 (possibly p3, p4) is not
    // a ghost particle. In almost all cases it is p1, however, it might be
    // other one. we call this particle reference particle.
    if (p1->l.ghost != 1) {
      // unfold non-ghost particle using image, because for physical
      // particles, the structure p->l.i is correctly set
      memmove(p11, p1->r.p, 3 * sizeof(double));
      memmove(img, p1->l.i, 3 * sizeof(int));
      unfold_position(p11, img);
      // other coordinates are obtained from its relative positions to the
      // reference particle
      get_mi_vector(AA, p2->r.p, p11);
      get_mi_vector(BB, p3->r.p, p11);
      for (int i = 0; i < 3; i++) {
	p22[i] = p11[i] + AA[i];
	p33[i] = p11[i] + BB[i];
      }
    } else {
      // in case the first particle is a ghost particle
      if (p2->l.ghost != 1) {
	memmove(p22, p2->r.p, 3 * sizeof(double));
	memmove(img, p2->l.i, 3 * sizeof(int));
	unfold_position(p22, img);
	get_mi_vector(AA, p1->r.p, p22);
	get_mi_vector(BB, p3->r.p, p22);
	for (int i = 0; i < 3; i++) {
	  p11[i] = p22[i] + AA[i];
	  p33[i] = p22[i] + BB[i];
	}
      } else {
	// in case the first and the second particle are ghost particles
	if (p3->l.ghost != 1) {
	  memmove(p33, p3->r.p, 3 * sizeof(double));
	  memmove(img, p3->l.i, 3 * sizeof(int));
	  unfold_position(p33, img);
	  get_mi_vector(AA, p1->r.p, p33);
	  get_mi_vector(BB, p2->r.p, p33);
	  for (int i = 0; i < 3; i++) {
	    p11[i] = p33[i] + AA[i];
	    p22[i] = p33[i] + BB[i];
	  }
	} else {
	  printf("Something wrong in oif_global_forces.hpp: All particles "
		 "in a bond are ghost particles, impossible to unfold the "
		 "positions...");
	  return 2;
	}
      }
    }
#endif
#ifndef GHOST_FLAG
    // if ghost flag was not defined we have no other option than to assume
    // the first particle is a physical one.
    memmove(p11, p1->r.p, 3 * sizeof(double));
    memmove(img, p1->l.i, 3 * sizeof(int));
    unfold_position(p11, img);
    // other coordinates are obtained from its relative positions to the
    // reference particle
    get_mi_vector(AA, p2->r.p, p11);
    get_mi_vector(BB, p3->r.p, p11);
    for (int i = 0; i < 3; i++) {
      p22[i] = p11[i] + AA[i];
      p33[i] = p11[i] + BB[i];
    }
#endif
    // unfolded positions correct
    /// starting code from volume force
    get_n_triangle(p11, p22, p33, VOL_norm);
    VOL_dn = normr(VOL_norm);
    VOL_A = area_triangle(p11, p22, p33);
    VOL_vv = (m_VOL_volume - m_V0) /
      m_V0;
    for (k = 0; k < 3; k++) {
      VOL_force[k] = m_kv * VOL_vv * VOL_A *
	VOL_norm[k] / VOL_dn * 1.0 / 3.0;
      // printf("%e ",force[k]);
      p1->f.f[k] += VOL_force[k];
      p2->f.f[k] += VOL_force[k];
      p3->f.f[k] += VOL_force[k];
    }
    ///  ending code from volume force

    for (k = 0; k < 3; k++) {
      h[k] = 1.0 / 3.0 * (p11[k] + p22[k] + p33[k]);
    }

    aa = (m_area - m_A0_g) /
      m_A0_g;

    // aminusb(3,h,p11,rh);				// area_forces for each
    // triangle node
    vecsub(h, p11, rh); // area_forces for each triangle node
    hn = normr(rh);
    for (k = 0; k < 3; k++) {
      force1[k] = m_ka_g * aa * rh[k] / hn;
      //(&part1)->f.f[k]+=force[k];
    }
    // aminusb(3,h,p22,rh);				// area_forces for each
    // triangle node
    vecsub(h, p22, rh); // area_forces for each triangle node
    hn = normr(rh);
    for (k = 0; k < 3; k++) {
      force2[k] = m_ka_g * aa * rh[k] / hn;
      //(&part2)->f.f[k]+=force[k];
    }
    // aminusb(3,h,p33,rh);				// area_forces for each
    // triangle node
    vecsub(h, p33, rh); // area_forces for each triangle node
    hn = normr(rh);
    for (k = 0; k < 3; k++) {
      force3[k] = m_ka_g * aa * rh[k] / hn;
      //(&part3)->f.f[k]+=force[k];
    }

    for (k = 0; k < 3; k++) {
      p1->f.f[k] += force1[k];
      p2->f.f[k] += force2[k];
      p3->f.f[k] += force3[k];
    };

  }
  else{
    return 2;
  };

  return 1;
}

boost::any Bond::OifGlobalForces::get_bond_parameters_from_bond() const
{

  Oif_global_forces_bond_parameters params = {m_A0_g, m_ka_g, m_V0, m_kv};
  return boost::any(params);
  
}
