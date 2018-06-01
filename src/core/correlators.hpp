
#include "correlators/Correlator.hpp"
#include <vector>
#include <memory>

namespace Correlators {

extern std::vector<std::shared_ptr<Correlators::Correlator>> auto_update_correlators;


void auto_update();

inline bool auto_update_enabled() {
  return auto_update_correlators.size() >0;
}

}
