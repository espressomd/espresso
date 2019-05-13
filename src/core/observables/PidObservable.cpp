#include "PidObservable.hpp"

#include "partCfg_global.hpp"

namespace Observables {
std::vector<double> PidObservable::operator()() const {
  return this->evaluate(partCfg());
}
} // namespace Observables