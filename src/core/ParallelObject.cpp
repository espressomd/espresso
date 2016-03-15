#include "MpiCallbacks.hpp"
#include "ParallelObject.hpp"

std::unordered_map<int, ParallelObject *> ParallelObject::address_table;

/** Ctor and Dtor run on all nodes in parallel */
ParallelObject::ParallelObject() {
 {
    using namespace std::placeholders;    
    MpiCallbacks::function_type f;
    /** Bind member function to this instance */
    f = std::bind(&ParallelObject::callback, this, _1, _2);

    m_callback_id = mpiCallbacks().add(f);

    address_table[m_callback_id] =this;        
  }
}

ParallelObject::~ParallelObject() {
  /** Remove the callback when deleting the object */
  mpiCallbacks().remove(m_callback_id);
}

void ParallelObject::call_slaves(int par1, int par2) {
  mpiCallbacks().call(m_callback_id, par1, par2);
}
 
int ParallelObject::get_id() const { return m_callback_id; }

ParallelObject *ParallelObject::get_local_address(int id) {
  return address_table.at(id);
}
