#include "config.hpp"
#ifdef CONSTRAINTS

#include "ObjectRegistry.hpp"
#include <vector>
#include <memory>
#include "constraints/Constraint.hpp"

namespace Constraints {
ObjectRegistry<std::vector<std::shared_ptr<Constraint>>> constraints;
}
#endif
