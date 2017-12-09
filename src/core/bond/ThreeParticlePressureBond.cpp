#include "ThreeParticlePressureBond.hpp"
#include "statistics.hpp" //obsstat_nonbonded function
#include "pressure.hpp" // virials, p_tensor

int Bond::ThreeParticlePressureBond::add_three_body_pressure(Particle *p1, int bl_id)
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
  if(auto bond_partners =  get_n_bond_partners<2>(p1, bl_id)){
    
    p2 = (*bond_partners)[0];
    p3 = (*bond_partners)[1];

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
  else{
    return 2;
  };

}
