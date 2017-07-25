#include "observables.hpp"

namespace Observables {
std::vector<std::shared_ptr<Observables::Observable>> auto_update_observables;

void auto_update() {
  for (auto& o : auto_update_observables) {
    o->calculate();
  }
}

void auto_write() {
  for (auto& o : auto_update_observables)
  {
    if ( o->writable() )
      o->write();
  }
}

bool auto_write_enabled() {
  for (auto &p : auto_update_observables)
  {
    if ( p->writable() )
      return true;
  }
  return false;
}

}

