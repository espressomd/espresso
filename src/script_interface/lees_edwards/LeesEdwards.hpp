
#include "../ScriptInterfaceBase.hpp"
#include "core/lees_edwards.hpp"

#ifdef LEES_EDWARDS

namespace ScriptInterface {
namespace LeesEdwards {

class LeesEdwards : public AutoParameters<LeesEdwards> {
public:
  LeesEdwards() {
    add_parameters(
	{{"type", lees_edwards_protocol.type},
	 {"time0", lees_edwards_protocol.time0},
         {"offset", lees_edwards_protocol.offset},
         {"velocity", lees_edwards_protocol.velocity},
         {"amplitude", lees_edwards_protocol.amplitude},
         {"frequency", lees_edwards_protocol.frequency},
         {"sheardir", lees_edwards_protocol.sheardir},
         {"shearplanenormal", lees_edwards_protocol.shearplanenormal}});
  };

}; // Class LeesEdwards

} // namespace LeesEdwards
} // namespace ScriptInterface

#endif
