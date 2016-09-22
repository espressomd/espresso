#include "ObjectRegistry.hpp"
#include <vector>
#include <memory>
#include "lbboundaries/LBBoundary.hpp"

namespace LBBoundaries {
ObjectRegistry<std::vector<std::shared_ptr<LBBoundary>>> lbboundaries;
}
