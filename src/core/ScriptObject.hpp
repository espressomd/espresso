#ifndef __SCRIPTOBJECT_HPP
#define __SCRIPTOBJECT_HPP

#include "Parameters.hpp"

class ScriptObject {
public:
  virtual const std::string name() const = 0;
  virtual void set_parameters(Parameters &parameters) = 0;
  virtual Parameters get_parameters() = 0;
  virtual Parameters &all_parameters() const = 0;
  virtual const Variant &get_parameter(std::string name) { return get_parameters()[name].value; };
  virtual void set_parameter(std::string name, const Variant &value) {
    Parameter p = all_parameters()[name];
    p.value = value;
    Parameters P;
    P[name] = p;
    set_parameters(P);
  }
};

#endif
