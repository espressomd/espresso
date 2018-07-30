#ifndef SCRIPTINTERFACE_ACCUMULATORS_ACCUMULATORBASE_HPP
#define SCRIPTINTERFACE_ACCUMULATORS_ACCUMULATORBASE_HPP

#include "ScriptInterface.hpp"
#include "core/accumulators/AccumulatorBase.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

namespace ScriptInterface {
namespace Accumulators {

class AccumulatorBase : public AutoParameters<AccumulatorBase> {
public:
  AccumulatorBase() {
    add_parameters({{"delta_N", [this](const Variant &v) {accumulator()->delta_N() = get_value<int>(v);}, // NOLINT, bug in clang-tidy-4.0
                                [this]() {return accumulator()->delta_N();}}});} // NOLINT, bug in clang-tidy-4.0
  virtual std::shared_ptr<const ::Accumulators::AccumulatorBase> accumulator() const = 0;
  virtual std::shared_ptr<::Accumulators::AccumulatorBase> accumulator() = 0;
};

}
}

#endif
