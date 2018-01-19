#ifndef SCRIPT_INTERFACE_AUTO_PARAMETERS_AUTO_PARAMETERS_HPP
#define SCRIPT_INTERFACE_AUTO_PARAMETERS_AUTO_PARAMETERS_HPP

#include <unordered_map>
#include <vector>

#include "AutoParameter.hpp"
#include "script_interface/ScriptInterfaceBase.hpp"

namespace ScriptInterface {

/**
 * @brief Bind parameters in the script interface.
 *
 * This class implementes @c ScriptInterfaceBase, binding
 * the parameters added by add_parameters or by the constructor.
 * To use, derive from this class and add parameters. For example,
 * given a class A
 * ~~~{.cpp}
 * class A {
 * public:
 *  int i() { return m_i; }
 * private:
 * int m_i;
 * };
 * ~~~
 * that should have i exposed, this can be achieved by extending it
 * like this:
 * ~~~{.cpp}
 * class A : public AutoParameters {
 * public:
 * A() : AutoParameters({"name_for_i", i}) {}
 *  int i() { return m_i; }
 * private:
 * int m_i;
 * };
 * ~~~
 *
 * If there is more complicated logic needed, specific setters and
 * getters can be provided. E.g. given a class B like
 * ~~~{.cpp}
 * class B {
 * public:
 * void set_i(int);
 * int get_i();
 * private:
 * int m_i;
 * };
 * ~~~
 * we can use a lambdas to set and get the parameter like this:
 * ~~~{.cpp}
 * class B : public AutoParameters {
 * public:
 * B() : AutoParameters({"name_for_i",
 * [this](Variant const& v) {
 *   set_i(get_value<int>(v));
 *  },
 * [this]() {
 *   return get_i();
 * }}) {}
 * void set_i(int);
 * int get_i();
 * private:
 * int m_i;
 * };
 * ~~~
 * (this has to be caputerd in the lambdas to have acces to the member functions
 * of the class).
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
          p.name,
          Parameter{p.type, p.length, p.set, p.get}));
    }
  }

public:
  /* ScriptInterfaceBase implementation */
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
