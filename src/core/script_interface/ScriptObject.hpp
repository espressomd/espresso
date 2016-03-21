#ifndef __SCRIPTOBJECT_HPP
#define __SCRIPTOBJECT_HPP

#include "Parameters.hpp"

#define SET_PARAMETER_HELPER(PARAMETER_NAME, MEMBER_NAME) if(name == PARAMETER_NAME) MEMBER_NAME = value

namespace ScriptInterface {

/**
 * @brief Base class for generic script interface.
 *
 * @TODO Add extensive documentation.
 *
 */
class ScriptObject {
 public:
  /**
   * @brief Human-readable name of the object.
   *
   * @return Name of the object.
   */   
  virtual const std::string name() const = 0;

  /**
   * @brief get current parameters.
   *
   * @return Parameters set in class.
   */
  virtual Parameters get_parameters() = 0;

  /**
   * @brief Get requiered and optional parameters for class
   *
   * Get requiered and optional parameters for class.
   * Should set default values for parameters if any.
   *
   * @return Expected parameters.
   */  
  virtual Parameters all_parameters() const = 0;
  /**
   * @brief Get single parameter.
   *
   * @param name Name of the parameter
   * @return Value of parameter @param name.
   */
  virtual const Variant get_parameter(const std::string &name) {
    return get_parameters()[name].value;
  };
  /**
   * @brief Set single parameter.
   *
   * @param name Name of the parameter
   * @param value Set parameter to this value.
   */  
  virtual void set_parameter(const std::string &name, const Variant &value) = 0;
  /**
   * @brief Set multiple parameters.
   *
   * @param Paramters Parameters to set.
   */  
  virtual void set_parameters(Parameters &parameters) {
    for(Parameters::iterator it = parameters.begin(); it != parameters.end(); ++it) {
      if(it->second.set != 0)
        set_parameter(it->first, it->second.value);
    }
  }
  /**
   * @brief Broadcast parameters to slave nodes.
   *
   * Broadcast parameters to slave nodes, only defined here
   * as NOOP so that ScriptObject and ParallelScriptObject
   * have the same interface.
   * Typcially the implementation from @class ParallelScriptobject
   * is used.
   */
  virtual void broadcast_parameters() {};
};

}

#endif
