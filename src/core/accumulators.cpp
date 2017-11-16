
#include "accumulators.hpp"
#include "integrate.hpp"

namespace Accumulators {
std::vector<std::shared_ptr<Accumulators::Accumulator>> auto_update_accumulators;

void auto_update() {
  for (auto& c : auto_update_accumulators) {
    c->update();
  }

}

}
