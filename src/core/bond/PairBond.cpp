#include "PairBond.hpp"
#include "grid.hpp" //get_mi_vector
#include "integrate.hpp" //integ_switch
#include "npt.hpp" //nptiso struct
#include "statistics.hpp" //obsstat_nonbonded function
#include "pressure.hpp" // virials, p_tensor

// force calculation
int Bond::PairBond::add_bonded_force(Particle *p1, int bl_id){
    // variables
    double dx[3] = {0., 0., 0.};
    double force[3] = {0., 0., 0.};
    int bond_broken;
    Particle* p2 = NULL;
#ifdef ROTATION
    double torque1[3] = {0., 0., 0.};
    double torque2[3] = {0., 0., 0.};
#endif
    
    // get the bond partners
    if(get_n_bond_partners(p1, bl_id) == 1){
      return 2;
    };
    p2 = m_bond_partners[0];

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
    };
#endif

#ifdef NPT
    if (integ_switch == INTEG_METHOD_NPT_ISO){
      for (int j=0;j++;j<3){
	nptiso.p_vir[j] += force[j] * dx[j];
      };
    };
#endif
    return bond_broken;

}

// energy calculation
int Bond::PairBond::add_bonded_energy(Particle *p1, int bl_id, double *_energy){

    // variables
    double dx[3] = {0., 0., 0.};
    Particle *p2 = NULL;

    // get bond partners
    if(get_n_bond_partners(p1, bl_id) == 1){
      return 2;
    };
    p2 = m_bond_partners[0];

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

int Bond::PairBond::add_virial(Particle *p1, int bl_id)
{

  double dx[3], force[3] = {0,0,0};
  Particle *p2;
  int bond_broken;
  int bond_map_id = p1->bl.e[bl_id];

  //get bond partners
  if(get_n_bond_partners(p1, bl_id) == 1){
    return 2;
  };
  p2 = m_bond_partners[0];  

  /* because of the NPT pressure calculation for pair forces, we need the
     1->2 distance vector here. For many body interactions this vector is
     not needed,
     and the pressure calculation not yet clear. */
  get_mi_vector(dx, p1->r.p, p2->r.p);

  bond_broken = add_bonded_pair_force(p1, p2, dx, force);

  *obsstat_bonded(&virials, bond_map_id) += dx[0]*force[0] + dx[1]*force[1] + dx[2]*force[2];
  
  /* stress tensor part */
  for(int k=0;k<3;k++)
    for(int l=0;l<3;l++)
      obsstat_bonded(&p_tensor, bond_map_id)[k*3 + l] += force[k]*dx[l];
  
  return bond_broken;
}

int Bond::PairBond::calc_pair_force(Particle *p1, Particle *p2, int bl_id, double force[3])
{

  double dx[3];
  int bond_broken;

  for(int i=0; i<3;i++){
    force[i] = 0.0;
  };

  //get partner
  if(get_n_bond_partners(p1, bl_id) == 1)
    {
      p2 = NULL;
      return 2;
    };
  // assign value
  p2 = m_bond_partners[0];
  // get dx
  get_mi_vector(dx, p1->r.p, p2->r.p);

  // calculate force
  bond_broken = add_bonded_pair_force(p1, p2, dx, force);

  // return
  if(bond_broken){return bond_broken;};
  return 0;

}
