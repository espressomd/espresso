#ifndef OIF_GLOBAL_FORCES_BOND_CLASS_H
#define OIF_GLOBAL_FORCES_BOND_CLASS_H
#include "ThreeParticleBond.hpp"

namespace Bond {

  class OifGlobalForces : public Bond {
    
  public:
    OifGlobalForces(double A0_g, double ka_g, double V0, double kv) : 
      Bond(2), m_A0_g(A0_g), m_ka_g(ka_g), m_V0(V0), m_kv(kv) 
    {m_bondtype = BondType::BONDED_IA_OIF_GLOBAL_FORCES;}

    /*first the area volume has to be calculated
      for this we have to loop over all oifglobalforce bonds*/
    int calc_oif_global(Particle* p1, int bl_id, double* partArea, double* VOL_partVol);
    //then the force can be calculated
    int add_bonded_force(Particle *p1, int bl_id) override;
    //There is no energy contribution
    int add_bonded_energy(Particle *p1, int bl_id) override;
    //set the static variables for calculating the force acting on the particles
    static void set_area_VOL(double area, double Vol){m_area=area;m_VOL_volume=Vol;};

  private:
    /*for every force m_area and Vol has to be calculated for all bonds
     every bond needs these two variables for calculating the force*/
    static double m_area;
    static double m_VOL_volume;
    //parameters for the bond
    double m_A0_g;
    double m_ka_g;
    double m_V0;
    double m_kv;

  };
}

#endif

