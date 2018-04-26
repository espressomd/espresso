#include <boost/mpi.hpp>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE ParallelScriptInterface test
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "core/MpiCallbacks.hpp"
#include "core/utils/make_unique.hpp"

#include "../../script_interface/ParallelScriptInterface.hpp"

namespace mpi = boost::mpi;
std::unique_ptr<Communication::MpiCallbacks> callbacks;

using namespace ScriptInterface;

struct TestClass : public ScriptInterfaceBase {
  TestClass() { constructed = true; }
  ~TestClass() { destructed = true; }

  void set_parameter(const std::string &name, const Variant &value) override {
    last_parameter = make_pair(name, value);

    if (name == "obj_param") {
      obj_param = get_instance(boost::get<ObjectId>(value)).lock();
    }
  }

  Variant get_parameter(std::string const &name) const override {
    if (name == "obj_param") {
      return obj_param->id();
    } else {
      return last_parameter.second;
    }
  }

  Variant call_method(const std::string &method,
                      const VariantMap &params) override {
    last_method_parameters = make_pair(method, params);

    return std::string("TestResult");
  }

  static std::pair<std::string, VariantMap> last_method_parameters;
  static std::pair<std::string, Variant> last_parameter;

  std::shared_ptr<ScriptInterfaceBase> obj_param;

  static bool constructed;
  static bool destructed;
};

bool TestClass::constructed = false;
bool TestClass::destructed = false;
std::pair<std::string, VariantMap> TestClass::last_method_parameters;
std::pair<std::string, Variant> TestClass::last_parameter;

/**
 * Check that instances are created and correctly destroyed on
 * the slave nodes.
 */
BOOST_AUTO_TEST_CASE(ctor_dtor) {
  /* Reset */
  TestClass::constructed = false;
  TestClass::destructed = false;

  if (callbacks->comm().rank() == 0) {
    /* Create an instance everywhere */
    auto so = std::make_shared<ParallelScriptInterface>("TestClass");
    /* Force destruction */
    so = nullptr;

    callbacks->abort_loop();
  } else {
    callbacks->loop();
  }

  /* Check that ctor and dtor were run on all nodes */
  BOOST_CHECK(TestClass::constructed);
  BOOST_CHECK(TestClass::destructed);
}

/**
 * Check that parameters are forwarded correctly.
 */
BOOST_AUTO_TEST_CASE(set_parameter) {
  if (callbacks->comm().rank() == 0) {
    auto so = std::make_shared<ParallelScriptInterface>("TestClass");

    so->set_parameter("TestParam", std::string("TestValue"));

    callbacks->abort_loop();
  } else {
    callbacks->loop();
  }

  auto const &last_parameter = TestClass::last_parameter;
  BOOST_CHECK(last_parameter.first == "TestParam");
  BOOST_CHECK(boost::get<std::string>(last_parameter.second) == "TestValue");
}

/*
 * Check that the method name and parameters are forwarded correctly
 * to the payload object, and check that the return value is
 * propagated correctly.
 */
BOOST_AUTO_TEST_CASE(call_method) {
  const VariantMap params{{"TestParam", std::string("TestValue")}};
  const std::string method{"TestMethod"};

  if (callbacks->comm().rank() == 0) {
    auto so = std::make_shared<ParallelScriptInterface>("TestClass");

    auto result = so->call_method(method, params);

    /* Check return value */
    BOOST_CHECK(boost::get<std::string>(result) == "TestResult");

    callbacks->abort_loop();
  } else {
    callbacks->loop();
  }

  auto const &last_parameters = TestClass::last_method_parameters;
  BOOST_CHECK(last_parameters.first == method);
  BOOST_CHECK(last_parameters.second == params);
}

BOOST_AUTO_TEST_CASE(parameter_lifetime) {
  if (callbacks->comm().rank() == 0) {
    auto host = std::make_shared<ParallelScriptInterface>("TestClass");
    ScriptInterfaceBase *bare_ptr;

    {
      auto parameter = ScriptInterfaceBase::make_shared(
          "TestClass", ScriptInterfaceBase::CreationPolicy::GLOBAL);
      bare_ptr = parameter.get();

      BOOST_CHECK(get_instance(parameter->id()) == parameter);

      host->set_parameter("obj_param", parameter->id());
    }

    auto param_id = host->get_parameter("obj_param");
    auto parameter = get_instance(param_id);

    /* Check that we got the original instance back */
    BOOST_CHECK(parameter.get() == bare_ptr);

    callbacks->abort_loop();
  } else {
    callbacks->loop();
  }
}

int main(int argc, char **argv) {
  mpi::environment mpi_env;
  mpi::communicator world;
  callbacks = Utils::make_unique<Communication::MpiCallbacks>(
      world, /* abort_on_exit */ false);

  ParallelScriptInterface::initialize(*callbacks);
  register_new<TestClass>("TestClass");

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
