#include "ThreeParticleBond.hpp" 
#include "energy.hpp"//for energy observable

int Bond::ThreeParticleBond::add_bonded_force(Particle *p1, int bl_id) {
  // variables
  double force[3] = {0., 0., 0.};
  double force2[3] = {0., 0., 0.};
  Particle *p2 = NULL, *p3 = NULL;
  int bond_broken;

  // get bond partners
  if(auto bond_partners =  get_n_bond_partners<2>(p1, bl_id)){

    p2 = (*bond_partners)[0];
    p3 = (*bond_partners)[1];

    bond_broken = calc_bonded_three_particle_force(p1, p2, p3, force, force2);

    if (bond_broken) {
      runtimeErrorMsg() << "bond broken between particles " << p1->p.identity
			<< ", " << p2->p.identity << " and "
			<< p3->p.identity;
      return bond_broken;
    }

    write_force_to_particle(p1, p2, p3, force, force2);
    return bond_broken;
  }
  else{
    return 2;
  };

}

int Bond::ThreeParticleBond::add_bonded_energy(Particle *p1, int bl_id) {

  Particle *p2 = NULL, *p3 = NULL;
  double _energy;
  // get bond partners
  if(auto bond_partners =  get_n_bond_partners<2>(p1, bl_id)){

    p2 = (*bond_partners)[0];
    p3 = (*bond_partners)[1];
  
    // calc energy
    int bond_broken =  calc_bonded_three_particle_energy(p1, p2, p3, &_energy);
    // add energy
    *obsstat_bonded(&energy, p1->bl.e[bl_id]) += _energy;
    // return
    return bond_broken;
  }
  else{
    return 2;
  };

}

// default function
void Bond::ThreeParticleBond::write_force_to_particle(Particle *p1, Particle *p2, Particle *p3, double force[3], double force2[3]) const {
  for(int j=0;j<3;j++){
    p1->f.f[j] += force[j];
    p2->f.f[j] += force2[j];
    p3->f.f[j] -= (force[j] + force2[j]);
  }
}

