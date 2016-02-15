#ifndef __SCRIPTOBJECT_HPP
#define __SCRIPTOBJECT_HPP

#include "Parameters.hpp"
#include "MpiCallbacks.hpp"

#include <iostream>
#include <functional>

#define SET_PARAMETER_HELPER(PARAMETER_NAME, MEMBER_NAME) if(name == PARAMETER_NAME) MEMBER_NAME = value

class ScriptObject {
public:
  ScriptObject() {
    using namespace std::placeholders;

    std::cout << this_node << ": ScriptObject()" << std::endl;
    MpiCallbacks::Function f;
    f = std::bind(&ScriptObject::callback, this, _1);
    m_callback_id = MpiCallbacks::add(f);
    std::cout << this_node << ": m_callback_id = " << m_callback_id << ")\n";
  }
  virtual void callback(int par) {
    std::cout << this_node << ": ScriptObject::callback(" << par << ")\n";
  }
  void call_slaves() {
    std::cout << this_node << ": ScriptObject::call_slaves(), m_callback_id = " << m_callback_id << ")\n";
    MpiCallbacks::call(m_callback_id, m_callback_id);
  }
  virtual const std::string name() const = 0;
  virtual void set_parameters(Parameters &parameters) {
    for(Parameters::iterator it = parameters.begin(); it != parameters.end(); ++it) {
      if(it->second.set != 0)
        set_parameter(it->first, it->second.value);
    }
  }
  
  virtual Parameters get_parameters() = 0;
  virtual Parameters all_parameters() const = 0;
  virtual const Variant &get_parameter(std::string name) { return get_parameters()[name].value; };
  virtual void set_parameter(const std::string &name, const Variant &value) = 0;
 private:
  int m_callback_id;
};

#endif
