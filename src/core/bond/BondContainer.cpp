#include "BondContainer.hpp"
#include "energy.hpp"//for energy observable
#include "pressure.hpp"//distribute_tensors->local_Stress_tensor arguments
#include "immersed_boundary/ibm_volume_conservation.hpp" //VolumesCurrent[]>

//help function
//cast base class into a derived class
template<class BaseClass, class DerivedClass>
DerivedClass* cast_base_class(BaseClass* base_class){

  return dynamic_cast<DerivedClass*>(base_class);

}

//insert bond
void Bond::BondContainer::set_bond_by_type(int type, std::unique_ptr<Bond> && bond)
{

  m_all_bonds.insert(std::pair<int, std::unique_ptr<Bond>>(type, std::move(bond)));
  sort_bond_into_lists(type);

}

//delete bond 
void Bond::BondContainer::delete_bond(int bond_map_id)
{

  m_force_loop_bonds.erase(bond_map_id);
  m_energy_loop_bonds.erase(bond_map_id);
  m_virial_loop_bonds.erase(bond_map_id);
  m_three_body_pressure_bonds.erase(bond_map_id);
  m_oif_global_forces_bonds.erase(bond_map_id);
  m_ibm_vol_con_bonds.erase(bond_map_id);
  m_rigid_bonds.erase(bond_map_id);
}

// private function which casts base classes into concrete classes
void Bond::BondContainer::sort_bond_into_lists(int bond_map_id)
{

  //try to find bond and cast it into derived classes
  if(m_all_bonds.count(bond_map_id) != 0){

    //---add Bonds with specific interface---

    //pairbond
    PairBond* derived_pair_bond = cast_base_class<Bond, PairBond>(m_all_bonds[bond_map_id].get());
    //3_particle_pressure bond
    ThreeParticlePressureBond* derived_3_pressure_bond = cast_base_class<Bond, 
									 ThreeParticlePressureBond>
      (m_all_bonds[bond_map_id].get());
    //oifglobal forces
    OifGlobalForces* derived_oif_global_bond = cast_base_class<Bond,
							       OifGlobalForces>
      (m_all_bonds[bond_map_id].get());
    //IbmVolumeConservation Bond
    IbmVolumeConservation* derived_ibm_vol_con_bond = cast_base_class<Bond,
								      IbmVolumeConservation>
      (m_all_bonds[bond_map_id].get());

    //rigid bond
    RigidBond* derived_rigid_bond = cast_base_class<Bond,
						    RigidBond>(m_all_bonds[bond_map_id].get());

    //add if cast was successful
    if(derived_pair_bond){
      m_virial_loop_bonds.insert(std::pair<int, PairBond*>(bond_map_id, derived_pair_bond));
    };
    if(derived_3_pressure_bond){
      m_three_body_pressure_bonds.insert(std::pair<int, ThreeParticlePressureBond*>
					 (bond_map_id, derived_3_pressure_bond));
    };
    if(derived_oif_global_bond){
      m_oif_global_forces_bonds.insert(std::pair<int, OifGlobalForces*>
				       (bond_map_id, derived_oif_global_bond));
    };
    if(derived_ibm_vol_con_bond){
      m_ibm_vol_con_bonds.insert(std::pair<int, IbmVolumeConservation*>
				 (bond_map_id, derived_ibm_vol_con_bond));
    };
    if(derived_rigid_bond){
      m_rigid_bonds.insert(std::pair<int, RigidBond*>(bond_map_id, derived_rigid_bond));
    };

    //---add bonds which add energies and forces---
    Bond* insert_to = m_all_bonds[bond_map_id].get();
    BondType bond_type = insert_to->get_Bond_Type();
    //some special bonds don't contribute to forces and/or energies
    switch(bond_type){
    case BondType::BONDED_IA_IBM_TRIBEND:
      //only force contribution
      m_force_loop_bonds.insert(std::pair<int,Bond*>(bond_map_id, insert_to));
      break;
    case BondType::BONDED_IA_IBM_TRIEL:
      //only force contribution
      m_force_loop_bonds.insert(std::pair<int,Bond*>(bond_map_id, insert_to));
      break;
    case BondType::BONDED_IA_OIF_OUT_DIRECTION:
      //only force contribution
      m_force_loop_bonds.insert(std::pair<int,Bond*>(bond_map_id, insert_to));
      break;
    case BondType::BONDED_IA_OIF_LOCAL_FORCES:
      //only force contribution
      m_force_loop_bonds.insert(std::pair<int,Bond*>(bond_map_id, insert_to));
      break;
    case BondType::BONDED_IA_VIRTUAL_BOND:
      //no force and energy contribution
      break;
    case BondType::BONDED_IA_OIF_GLOBAL_FORCES:
      //no energy contribution
      //force is calculated separately
      break;
    case BondType::BONDED_IA_IBM_VOLUME_CONSERVATION:
      //no energy contribution
      //force is calculated separately
      break;
    case BondType::BONDED_IA_RIGID_BOND:
      //neither direct energy nor force
      break;
    default:
      m_force_loop_bonds.insert(std::pair<int,Bond*>(bond_map_id, insert_to));
      m_energy_loop_bonds.insert(std::pair<int,Bond*>(bond_map_id, insert_to));
      break;
    };

  }
  else{
    runtimeErrorMsg() << "BondContainer.sort_bond_into_lists(): bond map id '" << bond_map_id << "' is unknown!";
    return;
  };
  
}

//force and energy calculation
int Bond::BondContainer::force_loop(Particle *p1)
{

  return loop_over_bond_partners(m_force_loop_bonds, &Bond::add_bonded_force, p1);
}

int Bond::BondContainer::energy_loop(Particle *p1)
{

  return loop_over_bond_partners(m_energy_loop_bonds, &Bond::add_bonded_energy, p1);

}

//pressure and stress tensor calculation
int Bond::BondContainer::virial_loop(Particle* p1)
{

  return loop_over_bond_partners(m_virial_loop_bonds, &PairBond::add_virial, p1);

}

int Bond::BondContainer::three_body_stress_loop(Particle *p1)
{

  return loop_over_bond_partners(m_three_body_pressure_bonds,
				 &ThreeParticlePressureBond::add_three_body_pressure, p1);

}

int Bond::BondContainer::local_stress_tensor_loop(Particle *p1, DoubleList *TensorInBin, int bins[3],
						  double range_start[3], double range[3])
{

  return loop_over_bond_partners(m_virial_loop_bonds, &PairBond::add_local_stress_tensor, p1,
				 TensorInBin, bins, range_start, range);

}

//oif global forces
int Bond::BondContainer::oif_global_loop(Particle *p1, double* partArea, double* VOL_partVol)
{

  return loop_over_bond_partners(m_oif_global_forces_bonds, &OifGlobalForces::calc_oif_global,
				 p1, partArea, VOL_partVol);

}

int Bond::BondContainer::oif_global_force_loop(Particle *p1)
{

  return loop_over_bond_partners(m_oif_global_forces_bonds, &OifGlobalForces::add_bonded_force,
				 p1);

}

//ibm volume conservation
int Bond::BondContainer::ibm_vol_con_softID_loop(Particle *p1, int *softID, int *bond_map_id)
{

  return loop_over_bond_partners(m_ibm_vol_con_bonds, &IbmVolumeConservation::get_soft_ID, p1,
				 softID, bond_map_id);

}

Bond::Bond* Bond::BondContainer::get_Bond(int bond_map_id)
{

  if(m_all_bonds.count(bond_map_id) == 0){
    return NULL;
  }
  else{
    return m_all_bonds[bond_map_id].get();
  };

}

Bond::IbmVolumeConservation* Bond::BondContainer::get_IBM_Vol_Con_Bond(int bond_map_id)
{

  if(m_ibm_vol_con_bonds.count(bond_map_id) == 0){
    return NULL;
  }
  else{
    return m_ibm_vol_con_bonds[bond_map_id];
  };
}

void Bond::BondContainer::init_Vol_Con()
{

  size_t size = m_ibm_vol_con_bonds.size();
  if(size == 0){return;};
  // Loop through all bonded interactions and check if we need to set the reference volume
  for(int i=0;i<size;i++){
    // This check is important because InitVolumeConservation may be called accidentally
    // during the integration. Then we must not reset the reference
    if(m_ibm_vol_con_bonds[i]->m_volRef == 0){
      int softID = m_ibm_vol_con_bonds[i]->m_softID;
      m_ibm_vol_con_bonds[i]->m_volRef = VolumesCurrent[softID];
      mpi_bcast_ia_params(i, -1);
    };
  }; 
}

//rigid bond
int Bond::BondContainer::RB_pos_corr(Particle *p1, int* repeat, int &cnt)
{
  return loop_over_bond_partners(m_rigid_bonds, &RigidBond::pos_corr, p1, repeat, cnt);
}

int Bond::BondContainer::RB_vel_corr(Particle *p1, int* repeat)
{
  return loop_over_bond_partners(m_rigid_bonds, &RigidBond::vel_corr, p1, repeat);
}
