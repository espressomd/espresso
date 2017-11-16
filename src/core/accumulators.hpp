#ifndef ESPRESSO_ACCUMULATORS_HPP
#define ESPRESSO_ACCUMULATORS_HPP

#include "accumulators/Accumulator.hpp"
#include <vector>
#include <memory>

namespace Accumulators {

extern std::vector<std::shared_ptr<Accumulators::Accumulator>> auto_update_accumulators;


void auto_update();

inline bool auto_update_enabled() {
  return auto_update_accumulators.size() >0;
}

}


#endif //ESPRESSO_ACCUMULATORS_HPP