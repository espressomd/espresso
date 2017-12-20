#include "BondContainer.hpp"
#include "energy.hpp"//for energy observable
#include "pressure.hpp"//distribute_tensors->local_Stress_tensor arguments

//cast base class into a derived class
template<class BaseClass, class DerivedClass>
DerivedClass* cast_base_class(BaseClass* base_class){

  return dynamic_cast<DerivedClass*>(base_class);

}

void Bond::BondContainer::set_bond_by_type(int type, std::unique_ptr<Bond> && bond)
{

  m_all_bonds.insert(std::pair<int, std::unique_ptr<Bond>>(type, std::move(bond)));
  sort_bond_into_lists(type);
}

// private function which casts base classes into concrete classes
void Bond::BondContainer::sort_bond_into_lists(int type)
{

  //try to find bond and cast it into derived classes
  if(m_all_bonds.count(type) != 0){

    //---add Bonds with specific interface---

    //pairbond
    PairBond* derived_pair_bond = cast_base_class<Bond, PairBond>(m_all_bonds[type].get());
    //3_particle_pressure bond
    ThreeParticlePressureBond* derived_3_pressure_bond = cast_base_class<Bond, 
									 ThreeParticlePressureBond>
      (m_all_bonds[type].get());
    //oifglobal forces
    OifGlobalForces* derived_oif_global_bond = cast_base_class<Bond,
							       OifGlobalForces>
      (m_all_bonds[type].get());

    //add if cast was successful
    if(derived_pair_bond){
      m_virial_loop_bonds.insert(std::pair<int, PairBond*>(type, derived_pair_bond));
    };
    if(derived_3_pressure_bond){
      m_three_body_pressure_bonds.insert(std::pair<int, ThreeParticlePressureBond*>
					 (type, derived_3_pressure_bond));
    };
    if(derived_oif_global_bond){
      m_oif_global_forces_bonds.insert(std::pair<int, OifGlobalForces*>
				       (type, derived_oif_global_bond));
    };

    //---add bonds which add energies and forces---
    Bond* insert_to = m_all_bonds[type].get();
    BondType bond_type = insert_to->get_Bond_Type();
    //some special bonds don't contribute to forces and/or energies
    switch(bond_type){
    case BondType::BONDED_IA_IBM_TRIBEND:
      //only force contribution
      m_force_loop_bonds.insert(std::pair<int,Bond*>(type, insert_to));
      break;
    case BondType::BONDED_IA_IBM_TRIEL:
      //only force contribution
      m_force_loop_bonds.insert(std::pair<int,Bond*>(type, insert_to));
      break;
    case BondType::BONDED_IA_OIF_OUT_DIRECTION:
      //only force contribution
      m_force_loop_bonds.insert(std::pair<int,Bond*>(type, insert_to));
      break;
    case BondType::BONDED_IA_OIF_LOCAL_FORCES:
      //only force contribution
      m_force_loop_bonds.insert(std::pair<int,Bond*>(type, insert_to));
      break;
    case BondType::BONDED_IA_VIRTUAL_BOND:
      //no force and energy contribution
      break;
    case BondType::BONDED_IA_OIF_GLOBAL_FORCES:
      //no energy contribution
      //force is calculated separately
      break;
    default:
      m_force_loop_bonds.insert(std::pair<int,Bond*>(type, insert_to));
      m_energy_loop_bonds.insert(std::pair<int,Bond*>(type, insert_to));
      break;
    };

  }
  else{
    runtimeErrorMsg() << "sort_bond_into_lists(): bond type '" << type << "' is unknown!";
    return;
  };
  
}

int Bond::BondContainer::force_loop(Particle *p1)
{

  return loop_over_bond_partners(m_force_loop_bonds, &Bond::add_bonded_force, p1);
}

int Bond::BondContainer::energy_loop(Particle *p1)
{

  return loop_over_bond_partners(m_energy_loop_bonds, &Bond::add_bonded_energy, p1);

}

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

int Bond::BondContainer::oif_global_loop(Particle *p1, double* partArea, double* VOL_partVol)
{

  return loop_over_bond_partners(m_oif_global_forces_bonds, &OifGlobalForces::calc_oif_global,
				 p1, partArea, VOL_partVol);

};

int Bond::BondContainer::oif_global_force_loop(Particle *p1)
{

  return loop_over_bond_partners(m_oif_global_forces_bonds, &OifGlobalForces::add_bonded_force,
				 p1);

};

