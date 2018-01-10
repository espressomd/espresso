#ifndef SCRIPT_INTERFACE_SERIALIZER_HPP
#define SCRIPT_INTERFACE_SERIALIZER_HPP

#include "get_value.hpp"

#include <boost/variant/static_visitor.hpp>

namespace ScriptInterface {
/**
 * @brief Serialize a Variant into a Variant with type info.
 *
 * ObjectId values are flattend by the get_state function of
 * the ScriptObject they refer to.
 */
class Serializer : public boost::static_visitor<Variant> {
public:
  template <typename T> Variant operator()(T const &val) const {
    return std::vector<Variant>{{static_cast<int>(infer_type<T>()), val}};
  }

  Variant operator()(ObjectId const &oid) const {
    auto so_ptr = get_value<std::shared_ptr<ScriptInterfaceBase>>(oid);
    if (so_ptr) {
      return std::vector<Variant>{
          {static_cast<int>(VariantType::OBJECTID), so_ptr->name(),
           static_cast<int>(so_ptr->policy()), so_ptr->get_state()}};
    } else {
      return std::vector<Variant>{{static_cast<int>(VariantType::NONE)},
                                  None{}};
    }
  }
};

/**
 * @brief Serialize a Variant into a Variant with type info.
 *
 * ObjectId values are flattend by the get_state function of
 * the ScriptObject they refer to.
 */
class UnSerializer : public boost::static_visitor<Variant> {
  std::vector<std::shared_ptr<ScriptInterfaceBase>> m_created_objects;

public:
  std::vector<std::shared_ptr<ScriptInterfaceBase>> const &
  created_objects() const {
    return m_created_objects;
  }

  template <typename T> Variant operator()(T const &val) {
    throw std::runtime_error("Invalid format.");
  }

  Variant operator()(std::vector<Variant> const &val) {
    using boost::get;
    switch (val.size()) {
    case 2: /* Normal value */
      return val[1];
      break;
    case 4: /* Object value */
    {
      auto so_ptr = ScriptInterfaceBase::make_shared(
          get<std::string>(val[1]),
          ScriptInterfaceBase::CreationPolicy(get<int>(val[2])), val[3]);
      /* Store a copy to keep the so alive. */
      m_created_objects.push_back(so_ptr);

      return so_ptr->id();
    } break;
    default: /* Error */
      throw std::runtime_error("Invalid format.");
    }
  }
};
}
#endif
