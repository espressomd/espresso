#include "lees_edwards_protocol.hpp"
#include <memory>

namespace LeesEdwards {

std::shared_ptr<ActiveProtocol> active_protocol;

double pos_offset_at_last_resort = 0;
} // namespace LeesEdwards
