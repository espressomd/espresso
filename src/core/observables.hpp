#pragma once

#include "observables/Observable.hpp"
#include <vector>
#include <memory>

namespace Observables {

extern std::vector<std::shared_ptr<Observables::Observable>> auto_update_observables;

void auto_update();
void auto_write();
bool auto_write_enabled();

inline bool auto_update_enabled() {
  return auto_update_observables.size() >0;
}

}
