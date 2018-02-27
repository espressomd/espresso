#ifndef IBM_VOLUME_CONSERVATION_BOND_H
#define IBM_VOLUME_CONSERVATION_BOND_H

#include "Bond.hpp"

namespace Bond{

  class IbmVolumeConservation : public Bond {

  public:
    IbmVolumeConservation(int softID, double volRef, double kappaV) : 
      Bond(1), m_softID{softID}, m_volRef{volRef}, m_kappaV{kappaV} 
    {m_bondtype = BondType::BONDED_IA_IBM_VOLUME_CONSERVATION;}

    //virtual functions of bond
    int add_bonded_force(Particle *p1, int bl_id) override;
    int add_bonded_energy(Particle *p1, int bl_id) override;
    
    //specific functions of this bond
    //get soft ID of the bond
    //bl_id not necesary but we need it because we use this function in BondContainer
    int get_soft_ID(Particle *p1, int bl_id, int *softID, int *bond_map_id);
    //to be called directly without BondContainer
    int calc_volumes(Particle *p1, int bl_id, double *tempVol);//bl_id is here now the id for a IBM_Triel_Bend_Bond!
    void ResetParams(double VolRef);

    int m_softID; // ID of the large soft particle to which this node belongs
    // Reference volume
    double m_volRef;
    // Spring constant for volume force
    double m_kappaV;
    // Whether to write out center-of-mass at each time step
    // Actually this is more of an analysis function and does not strictly belong
    // to volume conservation
    //  bool writeCOM;

  };

}

#endif
