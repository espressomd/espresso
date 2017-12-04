#include "BondContainer.hpp"
#include "energy.hpp"//for energy observable
#include "pressure.hpp"//distribute_tensors->local_Stress_tensor arguments

/*
  zum casten
  derived_ptr = dynamic_cast<Derived *>(base_ptr);

  if(derived_ptr) {

  }

  template<typename Base, typename Derived>
  smart_ptr<Derived> std::dynamic_pointer_cast<Derived>(smart_ptr<Base> base) {
  auto bare_ptr = base.get();

  return smart_ptr<Derived>(bare_ptr);
  }
*/

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
    //add if cast was successful
    if(derived_pair_bond){
      m_virial_loop_bonds.insert(std::pair<int, PairBond*>(type, derived_pair_bond));
    };
    if(derived_3_pressure_bond){
      m_three_body_pressure_bonds.insert(std::pair<int, ThreeParticlePressureBond*>
					 (type, derived_3_pressure_bond));
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
    case BondType::BONDED_IA_VIRTUAL_BOND:
      //no force and energy contribution
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

  /*
  // variables
  // bond_list_id: id of bond_map_id in p1->bl.e[]
  // bond_map_id: id of bond in bond_map
  // n_partners: number of partners, e.g. for a pair bond -> 1
  int i, bond_list_id, bond_map_id, n_partners, bond_broken;
  i = 0;
  //We are going through the bond list of particle p1
  while(i < p1->bl.n){

    //--first assign the variables for this step--
    //get the bond id in bond list
    bond_list_id = i;
    //get the bond map id
    bond_map_id = p1->bl.e[bond_list_id];

    //now the number of partners can be determined
    //try to find bond
    if(m_all_bonds.count(bond_map_id)==0){
      runtimeErrorMsg() << "add_bonded_force: bond type of atom "
                        << p1->p.identity << " unknown\n";
      return 1;
    }
    
    // if there is a type, get the number of bond partners
    n_partners = m_all_bonds[bond_map_id]->get_number_of_bond_partners();

    //try to find bond which has no vanishing force exerted on particles
    if(m_force_loop_bonds.count(bond_map_id)==0){
      i += n_partners + 1;
      continue;
    }
    
    //--Then calculate forces--
    //in this function the forces are written directly into the particles
    //because it determines the bond partners by itself
    bond_broken = m_force_loop_bonds[bond_map_id]->add_bonded_force(p1, bond_list_id);
    // if there are no bond partners get out of this function
    // it is easier to return 1 for the future!
    if(bond_broken == 2)
      return bond_broken;

    //--Now we are finished and have to go to the next bond--
    //in bond list: bond itself and number of partners 
    //next bond in -> m_partners + 1
    i += n_partners + 1;
    
  };

  return 0;
  */
}

int Bond::BondContainer::energy_loop(Particle *p1)
{

  return loop_over_bond_partners(m_energy_loop_bonds, &Bond::add_bonded_energy, p1);

  /*

  // variables
  // bond_list_id: id of bond_map_id in p1->bl.e[]
  // bond_map_id: id of bond in bond_map
  // n_partners: number of partners, e.g. for a pair bond -> 1
  int i, bond_list_id, bond_map_id, n_partners, bond_broken;
  //energy of that bond
  double ret;
  i = 0;
  //We are going through the bond list of particle p1
  while(i < p1->bl.n){

    //--first assign the variables for this step--
    //get the bond id in bond list
    bond_list_id = i;
    //get the bond map id
    bond_map_id = p1->bl.e[bond_list_id];

    //now the number of partners can be determined
    //try to find bond
    if(m_all_bonds.count(bond_map_id)==0){
      runtimeErrorMsg() << "add_bonded_energy: bond type of atom "
                        << p1->p.identity << " unknown\n";
      return 1;
    }
    
    // if there is a type, get the number of bond partners
    n_partners = m_all_bonds[bond_map_id]->get_number_of_bond_partners();

    //try to find bond which has no vanishing energy exerted on particles
    if(m_energy_loop_bonds.count(bond_map_id)==0){
      i += n_partners + 1;
      continue;
    }
    
    //--Then calculate energys--
    bond_broken = m_energy_loop_bonds[bond_map_id]->add_bonded_energy(p1, bond_list_id, &ret);
    //--Finally add that value to global variable
    *obsstat_bonded(&energy, bond_map_id) += ret;
    // if there are no bond partners get out of this function
    // it is easier to return 1 for the future!
    if(bond_broken == 2)
      return bond_broken;

    //--Now we are finished and have to go to the next bond--
    //in bond list: bond itself and number of partners 
    //next bond in -> m_partners + 1
    i += n_partners + 1;
    
  };

  return 0;
  */
}

int Bond::BondContainer::virial_loop(Particle* p1)
{

  return loop_over_bond_partners(m_virial_loop_bonds, &PairBond::add_virial, p1);
  /*
  int i, bond_list_id, bond_broken, bond_map_id, n_partners;
  i=0;
  while(i<p1->bl.n){

    //--first assign the variables for this step--
    //get the bond id in bond list
    bond_list_id = i;
    //get the bond map id
    bond_map_id = p1->bl.e[bond_list_id];

    //try to find bond
    if(m_all_bonds.count(bond_map_id)==0){
      runtimeErrorMsg() << "add_bonded_virials(): bond type of atom "
                        << p1->p.identity << " unknown\n";
      return 1;
    }

    // if there is a type, get the number of bond partners
    n_partners = m_all_bonds[bond_map_id]->get_number_of_bond_partners();

    //try to find bond which has no vanishing energy exerted on particles
    if(m_virial_loop_bonds.count(bond_map_id)==0){
      i += n_partners + 1;
      continue;
    }


    //--Calculate pressure--
    bond_broken = m_virial_loop_bonds[bond_map_id]->add_virial(p1, bond_list_id);

    // if there are no bond partners get out of this function
    // it is easier to return 1 for the future!
    if(bond_broken == 2)
       return 2;

    //--Now we are finished and have to go to the next bond--
    //in bond list: bond itself and number of partners 
    //next bond in -> n_partners + 1
    i += n_partners + 1;
    
  };

  return 0;
*/
}

int Bond::BondContainer::three_body_stress_loop(Particle *p1)
{

  return loop_over_bond_partners(m_three_body_pressure_bonds,
				 &ThreeParticlePressureBond::add_three_body_pressure, p1);
  /*
  int i, bond_list_id, bond_broken, bond_map_id, n_partners;
  i=0;
  while(i<p1->bl.n){

    //--first assign the variables for this step--
    //get the bond id in bond list
    bond_list_id = i;
    //get the bond map id
    bond_map_id = p1->bl.e[bond_list_id];
    if(m_all_bonds.count(bond_map_id) == 0){
      runtimeErrorMsg() << "add_three_body_bonded_stress(): bond type of atom "
                        << p1->p.identity << " unknown\n";
      return 1;
    }
    //now the number of partners can be determined
    n_partners = m_all_bonds[bond_map_id]->get_number_of_bond_partners();

    //try to find bond which has no vanishing energy exerted on particles
    if(m_three_body_pressure_bonds.count(bond_map_id)==0){
      i += n_partners + 1;
      continue;
    }
    
    //--Calculate pressure--
    bond_broken = m_three_body_pressure_bonds[bond_map_id]->add_three_body_pressure
      (p1, bond_list_id);

    // if there are no bond partners get out of this function
    // it is easier to return 1 for the future!
    if(bond_broken == 2)
       return 2;

    //--Now we are finished and have to go to the next bond--
    //in bond list: bond itself and number of partners 
    //next bond in -> n_partners + 1
    i += n_partners + 1;
    
  };

  return 0;*/
}

int Bond::BondContainer::local_stress_tensor_loop(Particle *p1, DoubleList *TensorInBin, int bins[3],
						  double range_start[3], double range[3])
{

  return loop_over_bond_partners(m_virial_loop_bonds, &PairBond::add_local_stress_tensor, p1,
				 TensorInBin, bins, range_start, range);

  /*
  int j, bond_list_id, bond_map_id, n_partners, bond_broken;
  double force[3];
  Particle *p2 = NULL;
  j = 0;
  while(j < p1->bl.n){

    bond_list_id = j;
    bond_map_id = p1->bl.e[bond_list_id];
    
    if(m_all_bonds.count(bond_map_id) == 0){
      runtimeErrorMsg() << "local_stress_tensor_calc(): bond type of atom "
                        << p1->p.identity << " unknown\n";
      return 1;
    };

    n_partners = m_all_bonds[bond_map_id]->get_number_of_bond_partners();

    if(m_virial_loop_bonds.count(bond_map_id)==0){
      j += n_partners + 1;
      continue;
    }

    bond_broken = m_virial_loop_bonds[bond_map_id]->calc_pair_force(p1, p2, bond_list_id, force);

    if(p2 != NULL){
      //write into stress tensor
      PTENSOR_TRACE(fprintf(stderr,"%d: Bonded to particle %d with force %f %f %f\n",
			    this_node,p2->p.identity,force[0],force[1],force[2]));
      if ((pow(force[0],2)+pow(force[1],2)+pow(force[2],2)) > 0) {

	if (distribute_tensors(TensorInBin,force,bins,
			       range_start,range,p1->r.p, p2->r.p) != 1) return 0;

      };

    };

    j+= n_partners + 1;
  };

  return 1;
  */
}
