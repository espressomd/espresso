#include "FourParticleBond.hpp"


int Bond::FourParticleBond::add_bonded_force(Particle *p1, int bl_id){
  // variables
  double force[3] = {0., 0., 0.};
  double force2[3] = {0., 0., 0.};
  double force3[3] = {0., 0., 0.};
  double force4[3] = {0., 0., 0.};
  Particle *p2 = NULL, *p3 = NULL, *p4 = NULL;
  int bond_broken;

  // get bond partners
  if(get_n_bond_partners(p1, bl_id) == 1){
    return 2;
  };
  p2 = m_bond_partners[0];
  p3 = m_bond_partners[1];
  p4 = m_bond_partners[2];

  bond_broken = add_bonded_four_particle_force(p1, p2, p3, p4, force, force2, force3, force4);

  if (bond_broken) {
    runtimeErrorMsg() << "bond broken between particles " << p1->p.identity
		      << ", " << p2->p.identity << ", " << p3->p.identity
		      << " and " << p4->p.identity;
    return bond_broken;
  };

  write_force_to_particle(p1, p2, p3, p4, force, force2, force3, force4);
  return bond_broken;

}

int Bond::FourParticleBond::add_bonded_energy(Particle *p1, int bl_id, double *_energy){

  Particle *p2, *p3, *p4 = NULL;

  // get bond partners
  if(get_n_bond_partners(p1, bl_id) == 1){
    return 2;
  };
  p2 = m_bond_partners[0];
  p3 = m_bond_partners[1];
  p4 = m_bond_partners[2];

  return add_bonded_four_particle_energy(p1, p2, p3, p4, _energy);

}

// default for CG_DNA
void Bond::FourParticleBond::write_force_to_particle(Particle *p1, Particle *p2, Particle *p3, 
						     Particle *p4, double force[3], double force2[3],
						     double force3[3], double force4[3]) const
{
        for (int j = 0; j < 3; j++) {
          p1->f.f[j] += force[j];
          p2->f.f[j] += force2[j];
          p3->f.f[j] += force3[j];
          p4->f.f[j] += force4[j];
        }
}

