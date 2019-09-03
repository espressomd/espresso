#include "../ScriptInterfaceBase.hpp"
#include "config.hpp"
#include "core/lees_edwards.hpp"

#ifdef LEES_EDWARDS

namespace ScriptInterface {
namespace LeesEdwards {

class Protocol : public AutoParameters<Protocol> {
public:
  virtual std::shared_ptr<::LeesEdwards::ActiveProtocol> protocol() = 0;

}; // Class Protocol;

} // namespace LeesEdwards
} // namespace ScriptInterface

#endif
