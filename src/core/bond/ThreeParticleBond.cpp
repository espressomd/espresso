#include "ThreeParticleBond.hpp" 


int Bond::ThreeParticleBond::add_bonded_force(Particle *p1, int bl_id) const {
  // variables
  double force[3] = {0., 0., 0.};
  double force2[3] = {0., 0., 0.};
  Particle *p2 = NULL, *p3 = NULL;
  int bond_broken;

  /* fetch particle 2, which is always needed */
  p2 = local_particles[bl_id+1];
  if (!p2) {
    runtimeErrorMsg() << "bond broken between particles " << p1->p.identity
		      << " and " << p1->bl.e[bl_id+1]
		      << " (particles are not stored on the same node)";
    return 2;
  }

  /* fetch particle 3 eventually */
  p3 = local_particles[bl_id+2];
  if (!p3) {
    runtimeErrorMsg() << "bond broken between particles " << p1->p.identity
		      << ", " << p1->bl.e[bl_id+1] << " and "
		      << p1->bl.e[bl_id+2]
		      << " (particles are not stored on the same node)";
    return 2;
  }

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

int Bond::ThreeParticleBond::add_bonded_energy(Particle *p1, int bl_id, double* _energy) const {

  Particle *p2 = NULL, *p3 = NULL;
  /* fetch particle 2, which is always needed */
  p2 = local_particles[bl_id+1];
  if (!p2) {
    runtimeErrorMsg() << "bond broken between particles " << p1->p.identity
		      << " and " << p1->bl.e[bl_id+1]
		      << " (particles are not stored on the same node)";
    return 2;
  }

  /* fetch particle 3 eventually */

  p3 = local_particles[bl_id+2];
  if (!p3) {
    runtimeErrorMsg() << "bond broken between particles " << p1->p.identity
		      << ", " << p1->bl.e[bl_id+1] << " and "
		      << p1->bl.e[bl_id+2]
		      << " (particles are not stored on the same node)";
    return 2;
  }
  
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

int Bond::ThreeParticleBond::get_number_of_bond_partners() const {
  return 2;
}
