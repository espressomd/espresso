#include "config.hpp"
#include "core/lees_edwards.hpp"
#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

namespace ScriptInterface {
namespace LeesEdwards {

class Protocol : public AutoParameters<Protocol> {
public:
  virtual std::shared_ptr<::LeesEdwards::ActiveProtocol> protocol() = 0;

}; // Class Protocol;

} // namespace LeesEdwards
} // namespace ScriptInterface
