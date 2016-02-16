#ifndef __PARALLEL_SCRIPTOBJECTP_HPP
#define __PARALLEL_SCRIPTOBJECTP_HPP

#include "ScriptObject.hpp"

class ParallelScriptObject : public ScriptObject {
 public:
  ParallelScriptObject() {
    using namespace std::placeholders;    
    MpiCallbacks::Function f;
    /** Bind member function to this instance */
    f = std::bind(&ScriptObject::callback, this, _1);

    m_callback_id = MpiCallbacks::add(f);
  }

 virtual void set_parameters(Parameters &parameters) {
    std::cout << this_node << ": ScriptObject::set_parameters(), m_callback_id = "
              << m_callback_id << ", size = " << parameters.size() << "\n";

    call_slaves(SET_PARAMETERS);

    boost::mpi::communicator comm;
    boost::mpi::broadcast(comm, parameters, 0);
    
    set_parameters_slave(parameters);
  }
  
 private:
  virtual void set_parameters_slave(Parameters &parameters) {
    std::cout << this_node << ": ScriptObject::set_parameters_slave(), m_callback_id = "
              << m_callback_id << ", size = " << parameters.size() << "\n";
    for(Parameters::iterator it = parameters.begin(); it != parameters.end(); ++it) {
      if(it->second.set != 0)
        set_parameter(it->first, it->second.value);
    }
  }

  virtual void set_parameter_slave(const std::string &name, const Variant &value) = 0;
  
  virtual void callback(int par) {
    std::cout << this_node << ": ScriptObject::callback(" << par << "), m_callback_id = " << m_callback_id << "\n";
    switch(par) {
      case SET_PARAMETERS:
        {
          Parameters param;
          boost::mpi::communicator comm;
          boost::mpi::broadcast(comm, param, 0);
          set_parameters_slave(param);          
        }
        break;
      case DELETE:
        delete this;
        return;
      default:
        break;
    }
  }
  virtual void call_slaves(int par) {
    std::cout << this_node << ": call_slaves(" << par << "), m_callback_id = " << m_callback_id << "\n";
    MpiCallbacks::call(m_callback_id, par);
  }

  ~ParallelScriptObject() {
    call_slaves(DELETE);

    delete this;
  }
  
  int m_callback_id;
  enum CallbackActions { DELETE, SET_PARAMETERS, SET_PARAMETER };
}

#endif
