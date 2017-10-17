#include "ThreeParticleBond.hpp" 
#include "statistics.hpp" //obsstat_nonbonded function
#include "pressure.hpp" // virials, p_tensor


int Bond::ThreeParticleBond::add_bonded_force(Particle *p1, int bl_id) {
  // variables
  double force[3] = {0., 0., 0.};
  double force2[3] = {0., 0., 0.};
  Particle *p2 = NULL, *p3 = NULL;
  int bond_broken;

  // get bond partners
  if(get_n_bond_partners(p1, bl_id) == 1){
    return 2;
  };
  p2 = m_bond_partners[0];
  p3 = m_bond_partners[1];

  bond_broken = add_bonded_three_particle_force(p1, p2, p3, force, force2);

  if (bond_broken) {
    runtimeErrorMsg() << "bond broken between particles " << p1->p.identity
		      << ", " << p2->p.identity << " and "
		      << p3->p.identity;
    return bond_broken;
  }

  write_force_to_particle(p1, p2, p3, force, force2);
  return bond_broken;

}

int Bond::ThreeParticleBond::add_bonded_energy(Particle *p1, int bl_id, double* _energy) {

  Particle *p2 = NULL, *p3 = NULL;

  // get bond partners
  if(get_n_bond_partners(p1, bl_id) == 1){
    return 2;
  };
  p2 = m_bond_partners[0];
  p3 = m_bond_partners[1];
  
  return add_bonded_three_particle_energy(p1, p2, p3, _energy);

}

// default function
void Bond::ThreeParticleBond::write_force_to_particle(Particle *p1, Particle *p2, Particle *p3, double force[3], double force2[3]) const {
  for(int j=0;j<3;j++){
    p1->f.f[j] += force[j];
    p2->f.f[j] += force2[j];
    p3->f.f[j] -= (force[j] + force2[j]);
  }
}

// as default: no force => no contribution to pressure
// exception: angle bonds and tabulated angle
int Bond::ThreeParticleBond::calc_3body_forces(Particle *p_mid, Particle *p_left,
					       Particle *p_right, double force1[3],
					       double force2[3], double force3[3]) const {
  for(int i = 0; i<3; i++){
    force1[i] = 0.0;
    force2[i] = 0.0;
    force3[i] = 0.0;    
  };
  return 0;

}

int Bond::ThreeParticleBond::add_three_body_pressure(Particle *p1, int bl_id)
{

  double dx12[3]; // espresso notation
  double dx21[3];
  double dx31[3];
  double force1[3];
  double force2[3];
  double force3[3];
  Particle *p2;
  Particle *p3;
  int bond_broken;
  int bond_map_id = p1->bl.e[bl_id];

  // get bond partners
  if(get_n_bond_partners(p1, bl_id) == 1){
    return 2;
  };
  p2 = m_bond_partners[0];
  p3 = m_bond_partners[1];

  // get distances
  get_mi_vector(dx12, p1->r.p, p2->r.p);
  for(int j = 0; j < 3; j++)
    dx21[j] = -dx12[j];

  get_mi_vector(dx31, p3->r.p, p1->r.p);

  for(int j = 0; j < 3; j++) {
    force1[j] = 0.0;
    force2[j] = 0.0;
    force3[j] = 0.0;
  }

  bond_broken = calc_3body_forces(p1, p2, p3, force1, force2, force3);

  /* three-body bonded interactions contribute to the stress but not the scalar pressure */
  for(int k = 0; k < 3; k++) {
    for(int l = 0; l < 3; l++) {
      obsstat_bonded(&p_tensor, 
		     bond_map_id)[3 * k + l] += force2[k] * dx21[l] + force3[k] * dx31[l];
    }
  }
  
  return bond_broken;

}
