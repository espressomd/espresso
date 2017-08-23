#include "EightParticleBond.hpp"

void Bond::EightParticleBond::write_force_to_particle(Particle *p1, Particle *p2, Particle *p3, 
						      Particle *p4, Particle *p5, Particle *p6, 
						      Particle *p7, Particle *p8, double force[3], 
						      double force2[3], double force3[3], 
						      double force4[3], double force5[3], 
						      double force6[3], double force7[3],
						      double force8[3]) const
{
  for (int j = 0; j < 3; j++) {
    p1->f.f[j] += force[j];
    p2->f.f[j] += force2[j];
    p3->f.f[j] += force3[j];
    p4->f.f[j] += force4[j];
    p5->f.f[j] += force5[j];
    p6->f.f[j] += force6[j];
    p7->f.f[j] += force7[j];
    p8->f.f[j] += force8[j];
  };
}

int Bond::EightParticleBond::get_number_of_bond_partners() const
{
  return 7;
}

int Bond::EightParticleBond::add_bonded_force(Particle *p1, int bl_id) const
{

  double force[3] = {0., 0., 0.};
  double force2[3] = {0., 0., 0.};
  double force3[3] = {0., 0., 0.};
  double force4[3] = {0., 0., 0.};
  double force5[3] = {0., 0., 0.};
  double force6[3] = {0., 0., 0.};
  double force7[3] = {0., 0., 0.};
  double force8[3] = {0., 0., 0.};
  Particle *p2, *p3, *p4, *p5, *p6, *p7, *p8 = NULL;
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

  p4 = local_particles[bl_id+3];
  if (!p4) {
    runtimeErrorMsg() << "bond broken between particles " << p1->p.identity
		      << ", " << p1->bl.e[bl_id+1] << ", " << p1->bl.e[bl_id+2]
		      << " and " << p1->bl.e[bl_id+3]
		      << " (particles not stored on the same node)";
    return 2;
  }

  p5 = local_particles[bl_id+4];
  p6 = local_particles[bl_id+5];
  p7 = local_particles[bl_id+6];
  p8 = local_particles[bl_id+7];

  if (!p4 || !p5 || !p6 || !p7 || !p8) {
    runtimeErrorMsg() << "bond broken between particles" << p1->p.identity
		      << ", " << p1->bl.e[bl_id+1] << ", " << p1->bl.e[bl_id+2]
		      << ", " << p1->bl.e[bl_id+3] << ", " << p1->bl.e[bl_id+4]
		      << ", " << p1->bl.e[bl_id+5] << ", " << p1->bl.e[bl_id+6]
		      << ", " << p1->bl.e[bl_id+7]
		      << " (particles not stored on the same node)";
    return 2; 

  };

    bond_broken = add_bonded_eight_particle_force(p1, p2, p3, p4, p5, p6, p7, p8,
						  force, force2, force3, force4,
						  force5, force6, force7, force8);

    if (bond_broken) {
      runtimeErrorMsg() << "bond broken between particles " << p1->p.identity
			<< ", " << p2->p.identity << ", " << p3->p.identity
			<< " and " << p4->p.identity;
      return bond_broken;
    };

    write_force_to_particle(p1, p2, p3, p4, p5, p6, p7, p8, force, force2, force3,
			    force4, force5, force6, force7, force8);
    return bond_broken;
}
int Bond::EightParticleBond::add_bonded_energy(Particle *p1, int bl_id, double* _energy) const
{

  Particle *p2, *p3, *p4, *p5, *p6, *p7, *p8 = NULL;

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

  p4 = local_particles[bl_id+3];
  if (!p4) {
    runtimeErrorMsg() << "bond broken between particles " << p1->p.identity
		      << ", " << p1->bl.e[bl_id+1] << ", " << p1->bl.e[bl_id+2]
		      << " and " << p1->bl.e[bl_id+3]
		      << " (particles not stored on the same node)";
    return 2;
  }

  p5 = local_particles[bl_id+4];
  p6 = local_particles[bl_id+5];
  p7 = local_particles[bl_id+6];
  p8 = local_particles[bl_id+7];

  if (!p4 || !p5 || !p6 || !p7 || !p8) {
    runtimeErrorMsg() << "bond broken between particles" << p1->p.identity
		      << ", " << p1->bl.e[bl_id+1] << ", " << p1->bl.e[bl_id+2]
		      << ", " << p1->bl.e[bl_id+3] << ", " << p1->bl.e[bl_id+4]
		      << ", " << p1->bl.e[bl_id+5] << ", " << p1->bl.e[bl_id+6]
		      << ", " << p1->bl.e[bl_id+7]
		      << " (particles not stored on the same node)";
    return 2; 

  };

  return add_bonded_eight_particle_energy(p1, p2, p3, p4, p5, p6, p7, p8, _energy);

}
