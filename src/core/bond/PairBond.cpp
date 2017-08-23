#include "PairBond.hpp"
#include "grid.hpp" //get_mi_vector

// force calculation
int Bond::PairBond::add_bonded_force(Particle *p1, int bl_id) const {
    // variables
    double dx[3] = {0., 0., 0.};
    double force[3] = {0., 0., 0.};
    int bond_broken;
    Particle* p2 = NULL;
#ifdef ROTATION
    double torque1[3] = {0., 0., 0.};
    double torque2[3] = {0., 0., 0.};
#endif

    /* fetch particle 2, which is always needed */
    p2 = local_particles[p1->bl.e[bl_id+1]];
    if (!p2) {
      runtimeErrorMsg() << "bond broken between particles " << p1->p.identity
                        << " and " << p1->bl.e[bl_id+1]
                        << " (particles are not stored on the same node)";
      return 2;
    } 

    /* because of the NPT pressure calculation for pair forces, we need the
         1->2 distance vector here. For many body interactions this vector is
         not needed,
         and the pressure calculation not yet clear. */
    get_mi_vector(dx, p1->r.p, p2->r.p);

    bond_broken = add_bonded_pair_force(p1, p2, dx, force);

    if (bond_broken) {
      runtimeErrorMsg() << "bond broken between particles " << p1->p.identity
			<< " and " << p2->p.identity
			<< ". Distance vector: " << dx[0] << " " << dx[1]
			<< " " << dx[2];
      return bond_broken;
    }
    // add forces to particles
    write_force_to_particle(p1, p2, force);

    // add torques to particles
#ifdef ROTATION
	  for (int j=0;j++;j<3){
	    p1->f.torque[j] += torque1[j];
	    p2->f.torque[j] += torque2[j];
	  }
#endif

#ifdef NPT
    if (integ_switch == INTEG_METHOD_NPT_ISO)
          nptiso.p_vir[j] += force[j] * dx[j];
#endif

    return bond_broken;

}

// energy calculation
int Bond::PairBond::add_bonded_energy(Particle *p1, int bl_id, double *_energy) const {

    // variables
    double dx[3] = {0., 0., 0.};

    /* fetch particle 2, which is always needed */
    Particle *p2 = local_particles[p1->bl.e[bl_id+1]];
    if (!p2) {
      runtimeErrorMsg() << "bond broken between particles " << p1->p.identity
                        << " and " << p1->bl.e[bl_id+1]
                        << " (particles are not stored on the same node)";
      return 2;
    } 

    /* because of the NPT pressure calculation for pair forces, we need the
         1->2 distance vector here. For many body interactions this vector is
         not needed,
         and the pressure calculation not yet clear. */
    get_mi_vector(dx, p1->r.p, p2->r.p);

    return add_bonded_pair_energy(p1, p2, dx, _energy);
}

// default function for writing forces to particles
void Bond::PairBond::write_force_to_particle(Particle *p1, Particle *p2, double force[3]) const {
      for (int j = 0; j < 3; j++) {
          p1->f.f[j] += force[j];
          p2->f.f[j] -= force[j];
      }

}


// Pair Bond -> 2 partners
int Bond::PairBond::get_number_of_bond_partners() const {
  return 1;
}

