#ifndef SCRIPT_INTERFACE_AUTO_PARAMETERS_AUTO_PARAMETERS_HPP
#define SCRIPT_INTERFACE_AUTO_PARAMETERS_AUTO_PARAMETERS_HPP

#include <unordered_map>
#include <vector>

#include "AutoParameter.hpp"
#include "script_interface/ScriptInterfaceBase.hpp"

namespace ScriptInterface {

/**
 * @brief Bind parameters in the script interface.
 */
class AutoParameters : public ScriptInterfaceBase {
public:
  /* Exceptions */
  struct UnknownParameter : public std::runtime_error {
    UnknownParameter(std::string const &name)
        : runtime_error("Parameter " + name + " is read-only.") {}
  };

  struct WriteError : public std::runtime_error {
    WriteError(std::string const &name)
        : runtime_error("Unknown parameter '" + name + "'.") {}
  };

protected:
  AutoParameters() = default;
  AutoParameters(std::vector<AutoParameter> &&params) {
    add_parameters(std::move(params));
  }

  void add_parameters(std::vector<AutoParameter> &&params) {
    for (auto const &p : params) {
      m_parameters.emplace(std::make_pair(
          std::move(p.name),
          Parameter{p.type, p.length, std::move(p.set), std::move(p.get)}));
    }
  }

public:
  /* ScriptInterfaceBase interface */
  ParameterMap valid_parameters() const final {
    ParameterMap valid_params;

    for (auto const &p : m_parameters) {
      valid_params.emplace(std::make_pair(
          p.first,
          ScriptInterface::Parameter{p.second.type, p.second.length, true}));
    }

    return valid_params;
  }

  Variant get_parameter(const std::string &name) const final {
    try {
      return m_parameters.at(name).get();
    } catch (std::out_of_range const &e) {
      throw UnknownParameter{name};
    }
  }

  void set_parameter(const std::string &name, const Variant &value) final {
    try {
      m_parameters.at(name).set(value);
    } catch (AutoParameter::WriteError const &e) {
      throw WriteError{name};
    } catch (std::out_of_range const &e) {
      throw UnknownParameter{name};
    }
  }

private:
  struct Parameter {
    VariantType type;
    size_t length;
    std::function<void(Variant const &)> set;
    std::function<Variant()> get;
  };

  std::unordered_map<std::string, Parameter> m_parameters;
};
}

#endif
