#include "BondBreakage.hpp"

#include "script_interface/Variant.hpp"
#include "script_interface/get_value.hpp"

#include "core/bond_breakage/bond_breakage.hpp" 
#include "core/system/System.hpp"

#include <string>

namespace ScriptInterface {
namespace BondBreakage {

Variant BondBreakage::do_call_method(std::string const &name, VariantMap const &) {
    if (name == "execute") {
        
        context()->parallel_try_catch([]() {    
        execute_bond_breakage(System::get_system(),*System::get_system().bond_breakage);
        });
        
        return {};
    }
    return {}; 
}

} // namespace BondBreakage
} // namespace ScriptInterface
