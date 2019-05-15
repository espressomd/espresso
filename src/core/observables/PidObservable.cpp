#include "PidObservable.hpp"

#include "partCfg_global.hpp"

namespace Observables {
std::vector<double> PidObservable::evaluate() const {
  return this->evaluate(partCfg());
}
} // namespace Observables