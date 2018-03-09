#include "observables.hpp"
#include "partCfg_global.hpp"

namespace Observables {

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

