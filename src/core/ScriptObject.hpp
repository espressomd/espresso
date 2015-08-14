#ifndef __SCRIPTOBJECT_HPP
#define __SCRIPTOBJECT_HPP

#include "Parameters.hpp"

class ScriptObject {
public:
  virtual const std::string name() const = 0;
  virtual void set_parameters(const Parameters &parameters) = 0;
  virtual Parameters get_parameters() = 0;
};

#endif
