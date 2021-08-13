#ifndef SCRIPT_INTERFACE_WALBERLA_EK_CONTAINER_HPP
#define SCRIPT_INTERFACE_WALBERLA_EK_CONTAINER_HPP

#include "EKSpecies.hpp"
#include "core/grid_based_algorithms/ek_container.hpp"

#include "script_interface/ObjectList.hpp"
#include "script_interface/ScriptInterface.hpp"

namespace ScriptInterface::walberla {

class EKContainer : public ObjectList<EKSpecies> {
  void add_in_core(std::shared_ptr<EKSpecies> const &obj_ptr) override {
    EK::ek_container.add(obj_ptr->get_ekinstance());
  }
  void remove_in_core(std::shared_ptr<EKSpecies> const &obj_ptr) override {
    EK::ek_container.remove(obj_ptr->get_ekinstance());
  }

  Variant do_call_method(std::string const &method,
                         VariantMap const &parameters) override {
    if (method == "set_tau") {
      EK::ek_container.set_tau(get_value<double>(parameters, "tau"));
      return none;
    }
    if (method == "get_tau") {
      return EK::ek_container.get_tau();
    }

    return ObjectList<EKSpecies>::do_call_method(method, parameters);
  }
};
} // namespace ScriptInterface::walberla

#endif // SCRIPT_INTERFACE_WALBERLA_EK_CONTAINER_HPP