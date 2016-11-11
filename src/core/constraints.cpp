#include "config.hpp"

#ifdef CONSTRAINTS

#include <memory>
#include <vector>

#include "ObjectRegistry.hpp"
#include "constraints/Constraint.hpp"

namespace Constraints {
ObjectRegistry<std::vector<std::shared_ptr<Constraint>>> constraints;
}
#endif
