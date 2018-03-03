#include"ThermalizedBond.hpp"
#include "random.hpp"
#include "debug.hpp"
#include "utils.hpp"
#include "interaction_data.hpp"
#include "integrate.hpp"

int Bond::ThermalizedBond::calc_bonded_pair_force(Particle *p1, Particle *p2, double dx[3],
						  double force[3]) const {

  //Bond broke?
  double dist = Utils::veclen(dx);
  if (m_r_cut > 0.0 && dist > m_r_cut) {
    return 1;
  }

  double force1[3], force2[3];
  double force_lv_com, force_lv_dist, com_vel, dist_vel;
  double mass_tot = p1->p.mass + p2->p.mass;
  double mass_tot_inv = 1.0 / mass_tot;
  double sqrt_mass_tot = sqrt(mass_tot);
  double sqrt_mass_red = sqrt(p1->p.mass * p2->p.mass / mass_tot);

  //printf("forcecalc pref1 %f pref2 %f\n", m_pref1_com, m_pref2_com);
  for (int i=0; i<3; i++) {
    //Langevin thermostat for center of mass
    com_vel = mass_tot_inv * (p1->p.mass * p1->m.v[i] + p2->p.mass * p2->m.v[i]);
    force_lv_com = -m_pref1_com * com_vel + sqrt_mass_tot * m_pref2_com * (d_random()-0.5);

    //Langevin thermostat for distance p1->p2
    dist_vel = p2->m.v[i] - p1->m.v[i];
    force_lv_dist =  -m_pref1_dist * dist_vel + sqrt_mass_red * m_pref2_dist * (d_random()-0.5);

    //Add forces
    force1[i] = p1->p.mass * mass_tot_inv * force_lv_com - force_lv_dist; 
    force2[i] = p2->p.mass * mass_tot_inv * force_lv_com + force_lv_dist; 
  }

  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: THERMALIZED BOND f = (%.3e,%.3e,%.3e)\n",this_node,p1->f.f[0]+force1[0],p1->f.f[1]+force1[1],p1->f.f[2]+force1[2]));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: THERMALIZED BOND f = (%.3e,%.3e,%.3e)\n",this_node,p2->f.f[0]+force2[0],p2->f.f[1]+force2[1],p2->f.f[2]+force2[2]));

  for (int j = 0; j < 3; j++) {
    p1->f.f[j] += force1[j];
    p2->f.f[j] += force2[j];
  };

      
  return 0;
  
}

int Bond::ThermalizedBond::calc_bonded_pair_energy(Particle *p1, Particle *p2, double dx[3],
						   double *_energy) const{

  return 0;
    
}

void Bond::ThermalizedBond::thermalized_bond_heat_up() {
  double pref_scale = sqrt(3);
  thermalized_bond_update_params(pref_scale);
}

void Bond::ThermalizedBond::thermalized_bond_cool_down() {
  double pref_scale = 1.0 / sqrt(3);
  thermalized_bond_update_params(pref_scale);
}

void Bond::ThermalizedBond::thermalized_bond_init()
{

  m_pref1_com = m_gamma_com / time_step;
  m_pref2_com = sqrt(24.0 * m_gamma_com / time_step * m_temp_com);
  m_pref1_dist = m_gamma_distance / time_step;
  m_pref2_dist = sqrt(24.0 * m_gamma_distance / time_step * m_temp_distance);
}


void Bond::ThermalizedBond::thermalized_bond_update_params(double pref_scale){
  
  m_pref2_com *= pref_scale;
  m_pref2_dist *= pref_scale;
	
}

void Bond::ThermalizedBond::write_force_to_particle(Particle *p1, Particle *p2, double force[3])
  const
{

  return;
  
}
