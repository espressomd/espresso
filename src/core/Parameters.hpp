#include "Variant.hpp"
#include <map>
#include <string>

struct Parameter {
  Parameter() : type(Variant::NONE), n_elements(0), set(false), required(false) {}
  Parameter(Variant::Type _type, Variant _value, int _n_elements, bool _set, bool _required) :
    type(_type), value(_value), n_elements(_n_elements), set(_set), required(_required) {}
  Parameter(Variant::Type _type, bool _req) : type(_type), value(_type), n_elements(0), set(false), required(_req) {}
  Parameter(Variant::Type _type, int _n, bool _req) : type(_type), value(_type), n_elements(_n), set(false), required(_req) {}

  /** @TODO: Should throw if types do not match */
  Parameter &operator=(const Variant& rhs)  {
    value = rhs;
    set = true;
    return *this;
  }
  
  Variant::Type type;
  Variant value;
  int n_elements;
  bool set, required;
};

class Parameters : public std::map<std::string, Parameter> {
public:
  virtual bool check() {
    bool ret = true;
    for(iterator it = begin(); it != end(); ++it) {
      ret &= (it->second.set | !it->second.required);
    }
    return ret;
  }
  virtual void unset() {
    for(iterator it = begin(); it != end(); ++it) {
      it->second.set = false;
    }
  }
};
