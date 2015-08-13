#ifndef __SCRIPTOBJECT_HPP
#define __SCRIPTOBJECT_HPP

#include "Parameters.hpp"

class ScriptObject {
public:
  virtual const std::string &name() const = 0;
  virtual void set_parameters(const Parameters &parameters) = 0;
  virtual const Parameters &get_parameters(void) = 0;
};

#endif
