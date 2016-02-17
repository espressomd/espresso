#include "MpiCallbacks.hpp"
#include "ParallelObject.hpp"

std::unordered_map<int, ParallelObject *> ParallelObject::address_table;

ParallelObject::ParallelObject() {
 {
    using namespace std::placeholders;    
    MpiCallbacks::Function f;
    /** Bind member function to this instance */
    f = std::bind(&ParallelObject::callback, this, _1);

    m_callback_id = MpiCallbacks::add(f);

    address_table[m_callback_id] =this;        
  }
}

void ParallelObject::call_slaves(int par) {
    MpiCallbacks::call(m_callback_id, par);
  }
 
int ParallelObject::get_id() const { return m_callback_id; }

ParallelObject *ParallelObject::get_local_address(int id) {
  return address_table.at(id);
}
