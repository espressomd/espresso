#ifndef __SCRIPTOBJECT_HPP
#define __SCRIPTOBJECT_HPP

#include <iostream>
#include <functional>

#include "Parameters.hpp"
#include "MpiCallbacks.hpp"

#include <boost/mpi.hpp>

#define SET_PARAMETER_HELPER(PARAMETER_NAME, MEMBER_NAME) if(name == PARAMETER_NAME) MEMBER_NAME = value

class ScriptObject {
public:  
  ScriptObject() {
    using namespace std::placeholders;    
    MpiCallbacks::Function f;
    /** Bind member function to this instance */
    f = std::bind(&ScriptObject::callback, this, _1);

    m_callback_id = MpiCallbacks::add(f);
  }
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
      default:
        break;
    }
  }
  void call_slaves(int par) {
    std::cout << this_node << ": call_slaves(" << par << "), m_callback_id = " << m_callback_id << "\n";
    MpiCallbacks::call(m_callback_id, par);
  }
  virtual const std::string name() const = 0;
  virtual void set_parameters(Parameters &parameters) {
    std::cout << this_node << ": ScriptObject::set_parameters(), m_callback_id = "
              << m_callback_id << ", size = " << parameters.size() << "\n";

    call_slaves(SET_PARAMETERS);

    boost::mpi::communicator comm;
    boost::mpi::broadcast(comm, parameters, 0);
    
    set_parameters_slave(parameters);
  }
  virtual void set_parameters_slave(Parameters &parameters) {
    std::cout << this_node << ": ScriptObject::set_parameters_slave(), m_callback_id = "
              << m_callback_id << ", size = " << parameters.size() << "\n";
    for(Parameters::iterator it = parameters.begin(); it != parameters.end(); ++it) {
      if(it->second.set != 0)
        set_parameter(it->first, it->second.value);
    }
  }
  
  virtual Parameters get_parameters() = 0;
  virtual Parameters all_parameters() const = 0;
  virtual const Variant get_parameter(const std::string &name) { return get_parameters()[name].value; };
  virtual void set_parameter(const std::string &name, const Variant &value) = 0;
 protected:
  enum CallbackActions { SET_PARAMETERS };
  ~ScriptObject() {
  }
 private:
  int m_callback_id;
};

#endif
