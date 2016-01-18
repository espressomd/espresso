#ifndef __SCRIPTOBJECT_HPP
#define __SCRIPTOBJECT_HPP

#include "Parameters.hpp"

#include <iostream>

#define SET_PARAMETER_HELPER(PARAMETER_NAME, MEMBER_NAME) if(name == PARAMETER_NAME) MEMBER_NAME = value

class ScriptObject {
public:
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
};

#endif
