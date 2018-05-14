#ifndef SCRIPT_INTERFACE_OBSERVABLES_LBPROFILEOBSERVABLE_HPP
#define SCRIPT_INTERFACE_OBSERVABLES_LBPROFILEOBSERVABLE_HPP

#include <memory>

#include "core/observables/LBProfileObservable.hpp"

#include "LBObservable.hpp"
#include "ProfileObservable.hpp"
#include "auto_parameters/AutoParameters.hpp"

namespace ScriptInterface {
namespace Observables {

template <typename CoreObs>
class LBProfileObservable : public ProfileObservable<CoreObs>,
                            public LBObservable<CoreObs> {
public:
  static_assert(
      std::is_base_of<::Observables::LBProfileObservable, CoreObs>::value, "");
  LBProfileObservable() : m_observable(std::make_shared<CoreObs>()) {}

  std::shared_ptr<::Observables::LBObservable> lb_observable() const {
    return m_observable;
  }

  std::shared_ptr<::Observables::ProfileObservable> profile_observable() const {
    return m_observable;
  }

  std::shared_ptr<::Observables::Observable> observable() const override {
    return m_observable;
  }

private:
  std::shared_ptr<CoreObs> m_observable;
};

} /* namespace Observables */
} /* namespace ScriptInterface */

#endif
