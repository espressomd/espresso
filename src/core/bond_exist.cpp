#include"bond_exist.hpp"
#include"interaction_data.hpp"
#include"bond.hpp"
bool make_bond_exist_bond_class(int bond_type){
  //total size of bond vector
  int bond_size = bonds_ia.size();
  //possible new size if bond_number is well chosen
  int ns = bond_type + 1;
  
  if(ns <= bond_size){
    // ifdef cases need to be here!!!
    return false;
  }
  else{
    // actually ns-1 is bond_type but this is for understanding algorithm
    for(int i = bond_size; i < ns-1; i++){
      NO_BOND* placeholder = new NO_BOND;
      bonds_ia.push_back(placeholder);
    };
    return true;
  };
}
