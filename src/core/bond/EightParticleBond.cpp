#include "EightParticleBond.hpp"
#include "energy.hpp"//for energy observable

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

int Bond::EightParticleBond::add_bonded_force(Particle *p1, int bl_id)
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

  // get bond partners
  if(auto bond_partners =  get_n_bond_partners<7>(p1, bl_id)){
    p2 = (*bond_partners)[0];
    p3 = (*bond_partners)[1];
    p4 = (*bond_partners)[2];
    p5 = (*bond_partners)[3];
    p6 = (*bond_partners)[4];
    p7 = (*bond_partners)[5];
    p8 = (*bond_partners)[6];
  
    bond_broken = calc_bonded_eight_particle_force(p1, p2, p3, p4, p5, p6, p7, p8,
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
  else{
    return 2;
    };
}

int Bond::EightParticleBond::add_bonded_energy(Particle *p1, int bl_id)
{

  Particle *p2, *p3, *p4, *p5, *p6, *p7, *p8 = NULL;
  double _energy;
  // get bond partners
  if(auto bond_partners =  get_n_bond_partners<7>(p1, bl_id)){
    p2 = (*bond_partners)[0];
    p3 = (*bond_partners)[1];
    p4 = (*bond_partners)[2];
    p5 = (*bond_partners)[3];
    p6 = (*bond_partners)[4];
    p7 = (*bond_partners)[5];
    p8 = (*bond_partners)[6];

    // calc energy
    int bond_broken = calc_bonded_eight_particle_energy(p1, p2, p3, p4, p5, p6, p7, p8, &_energy);
    // add energy
    *obsstat_bonded(&energy, p1->bl.e[bl_id]) += _energy;
    // return
    return bond_broken;
  }
  else{
    return 2;
  };
}
