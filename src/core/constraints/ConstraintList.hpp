#ifndef __CONSTRAINT_CONTAINER_HPP
#define __CONSTRAINT_CONTAINER_HPP

#include "Constraint.hpp"

#include <unordered_set>
#include <unordered_map>
#include <functional>
#include <memory>

#include "ParallelObject.hpp"
#include "initialize.hpp"

namespace Constraints {

class ConstraintList : public ParallelObject {  
 public:
  typedef typename std::unordered_set<Constraint *>::iterator iterator;
  ConstraintList();
  
  /** Add constraint (master only) */
  void add_constraint(std::shared_ptr<Constraint> c) {
    assert(c != nullptr);
    
  /** On the master we keep a copy of the shared_ptr to keep the
      object around. The existance of the objects on the slaves
      is coupled to the master. */
    int id = c->get_id();
    m_constraints_master.insert(std::make_pair(id, c));

    call_slaves(ADD, id);
    Constraint *c_p = get_pointer_from_id(id);
    m_constraints.insert(c_p);

    on_constraint_change();
  }

  /** Remove constraint (master only) */
  void remove_constraint(std::shared_ptr<Constraint> c) {
    int id = c->get_id();

    call_slaves(DELETE, id);    
    Constraint *c_p = get_pointer_from_id(id);
    m_constraints.erase(c_p);

    on_constraint_change();
  }

  iterator begin() {
    return m_constraints.begin();
  }
  iterator end() {
    return m_constraints.end();
  }
  
  
  void add_forces(Particle *p);
  void add_energies(Particle *p);
  void init_forces();
  double min_dist(double pos[3]);
  
 private:
  Constraint *get_pointer_from_id(const int &id) {
    /** Get pointer for local object.
        We know that this is actually a Constraint, because
        it was typechecked at the master */
    ParallelObject *p = ParallelObject::get_local_address(id);
    assert(p);
    Constraint *c = dynamic_cast<Constraint *>(p);
    assert(c);

    return c;
  }
  
  void callback(int action, int id) {
    switch(action) {
      case ADD:
        {
          /** Get id from master */
          Constraint *c = get_pointer_from_id(id);
          
          m_constraints.insert(c);

          on_constraint_change();
        }
        break;
      case DELETE:
        {
          /** Get id from master */
          int id;
          Constraint *c = get_pointer_from_id(id);

          m_constraints.erase(c);

          on_constraint_change();
        }
        break;
      default:
        break;
    }    
  }
  
  enum Action { ADD, DELETE };

  std::unordered_map<int, std::shared_ptr<Constraint> > m_constraints_master;
  std::unordered_set<Constraint *> m_constraints;
};

/**
 * Workaround to ensure construction order.
 * (master only). Will construct instances on
 * slaves and set the list_p to the local instance.
 */
ConstraintList &list();
}
  
#endif
