#pragma once

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"
#include <string>

namespace ScriptInterface {
namespace BondBreakage {

class BondBreakage : public AutoParameters<BondBreakage> {
public:
  Variant do_call_method(std::string const &name,
                         VariantMap const &parameters) override;
};

} // namespace BondBreakage
} // namespace ScriptInterface
