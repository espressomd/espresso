#ifndef ESPRESSO_CLASSTYPE_HPP
#define ESPRESSO_CLASSTYPE_HPP

namespace ScriptInterface {
/**
 * @brief Binding an interface name to C++ type.
 */
template <class T> struct ClassName {
  using class_type = T;
  const char *name;
};
} // namespace ScriptInterface

#endif // ESPRESSO_CLASSTYPE_HPP
