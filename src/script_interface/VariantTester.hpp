#include "ScriptInterface.hpp"

namespace ScriptInterface {
namespace Testing {
class VariantTester : public ScriptInterfaceBase {
public:
  Variant call_method(std::string const &method,
                      VariantMap const &par) override {
    if(method == "default") {
      return Variant{};
    }

    if (method == "true") {
      return true;
    }

    if (method == "false") {
      return false;
    }

    if (method == "flat") {
      return std::vector<Variant>{true,
                                  std::string("a string"),
                                  3.14159,
                                  std::vector<int>{3, 1, 4, 1, 5},
                                  std::vector<double>{1.1, 2.2, 3.3},
                                  ObjectId(),
                                  this->id()};
    }

    if (method == "recursive") {
      return make_recursive_variant(0, boost::get<int>(par.at("max_level")));
    }

    if (method == "mixed") {
      return std::vector<Variant>{1, std::string("another string"),
                                  make_recursive_variant(0, 5)};
    }

    if (method == "check_parameter_type") {
      auto const type = boost::get<std::string>(par.at("type"));

      if(type == "none")
        return is_none(par.at("value"));

      if (type == "bool") {
        return is_bool(par.at("value"));
      }

      if (type == "int") {
        return is_int(par.at("value"));
      }

      if (type == "double") {
        return is_double(par.at("value"));
      }

      if (type == "string") {
        return is_string(par.at("value"));
      }

      if (type == "objectid") {
        return is_objectid(par.at("value"));
      }

      if (type == "double_vector") {
        return is_double_vector(par.at("value"));
      }

      if (type == "int_vector") {
        return is_int_vector(par.at("value"));
      }

      if (type == "vector") {
        return is_vector(par.at("value"));
      }
    }

    throw std::runtime_error("Unknown method");
  }

private:
  Variant make_recursive_variant(int i, int max_level) const {
    if (i < max_level)
      return std::vector<Variant>{i, make_recursive_variant(i + 1, max_level)};
    else
      return std::string("end");
  }
};
}
}
