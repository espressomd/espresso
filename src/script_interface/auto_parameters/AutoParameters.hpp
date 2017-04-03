#ifndef SCRIPT_INTERFACE_AUTO_PARAMETERS_AUTO_PARAMETERS_HPP
#define SCRIPT_INTERFACE_AUTO_PARAMETERS_AUTO_PARAMETERS_HPP

#include <unordered_map>
#include <vector>

#include "AutoParameter.hpp"
#include "script_interface/ScriptInterfaceBase.hpp"

namespace ScriptInterface {

class AutoParameters : public ScriptInterfaceBase {
public:
  AutoParameters() = delete;

protected:
  AutoParameters(std::vector<AutoParameter> &&params) {
    for (auto const &p : params) {
      m_parameters.emplace(std::make_pair(
          p.name, Parameter{p.type, std::move(p.set), std::move(p.get)}));
    }
  }

public:
  ParameterMap valid_parameters() const final {
    ParameterMap valid_params;

    for (auto const &p : m_parameters) {
      valid_params.emplace(std::make_pair(
          p.first, ScriptInterface::Parameter{p.second.type, true}));
    }

    return valid_params;
  }

  Variant get_parameter(const std::string &name) const final {
    return m_parameters.at(name).get();
  }

  void set_parameter(const std::string &name, const Variant &value) final {
    m_parameters.at(name).set(value);
  }

private:
  struct Parameter {
    VariantType type;
    std::function<void(Variant const &)> set;
    std::function<Variant()> get;
  };

  std::unordered_map<std::string, Parameter> m_parameters;
};
}

#endif
