
#include "correlators.hpp"

namespace Correlators {
std::vector<std::shared_ptr<Correlators::Correlator>> auto_update_correlators;

void auto_update() {
for (auto& c : auto_update_correlators) {
  c->get_data();
}

}

}

