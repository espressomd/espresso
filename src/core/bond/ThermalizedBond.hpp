#ifndef THERMALIZED_DIST_BOND_H
#define THERMALIZED_DIST_BOND_H
#include "PairBond.hpp"
#include "CutoffBond.hpp"

namespace Bond {

  class ThermalizedBond : public PairBond, public CutoffBond {

  public:
    ThermalizedBond(double temp_com, double gamma_com, double temp_distance, double gamma_distance,
		    double r_cut, double pref1_com, double pref2_com, double pref1_dist,
		    double pref2_dist) : PairBond(), CutoffBond(r_cut), m_temp_com{temp_com},
      m_gamma_com{gamma_com}, m_temp_distance{temp_distance}, m_r_cut{r_cut}, m_pref1_com{pref1_com},
      m_pref2_com{pref2_com}, m_pref1_dist{pref1_dist}, m_pref2_dist{pref2_dist}
    {m_bondtype=BondType::BONDED_IA_THERMALIZED_DIST;}

    // Member function
    int calc_bonded_pair_force(Particle *p1, Particle *p2, double dx[3],
			       double force[3]) const override;
    int calc_bonded_pair_energy(Particle *p1, Particle *p2, double dx[3],
				double *_energy) const override;
    
    //force is written to particles in calc_bonded_pair_force
    void write_force_to_particle(Particle *p1, Particle *p2, double force[3]) const override;
    
    void thermalized_bond_heat_up();
    void thermalized_bond_cool_down();
    void thermalized_bond_update_params(double pref_scale);
    void thermalized_bond_init();

    boost::any get_bond_parameters_from_bond() const override;
    
  private:
    double m_temp_com;
    double m_gamma_com;
    double m_temp_distance;
    double m_gamma_distance;
    double m_r_cut;
    double m_pref1_com;
    double m_pref2_com;
    double m_pref1_dist;
    double m_pref2_dist;
  };
  
}

#endif
