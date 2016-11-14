#include <iostream>

#include <boost/mpi.hpp>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE MpiCallbacks test
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "MpiCallbacks.hpp"

#include "../../script_interface/ParallelScriptInterface.hpp"

using Communication::mpiCallbacks;
using namespace ScriptInterface;

namespace Testing {

void reduce_and_check(const boost::mpi::communicator &comm, bool local_value) {
  if (comm.rank() == 0) {
    bool total;
    boost::mpi::reduce(comm, local_value, total, std::logical_and<bool>(), 0);
    BOOST_CHECK(total);
  } else {
    boost::mpi::reduce(comm, local_value, std::logical_and<bool>(), 0);
  }
}
}

struct TestClass : public ScriptInterfaceBase {
  TestClass() {
    constructed = true;
    last_instance = this;
  }
  ~TestClass() { destructed = true; }

  void set_parameter(const std::string &name, const Variant &value) override {
    last_parameter = make_pair(name, value);
  }

  std::pair<std::string, Variant> last_parameter;
  static TestClass *last_instance;

  const std::string name() const override { return "TestClass"; }
  static bool constructed;
  static bool destructed;
};

bool TestClass::constructed = false;
bool TestClass::destructed = false;
TestClass *TestClass::last_instance = nullptr;

BOOST_AUTO_TEST_CASE(ctor_dtor) {
  if (mpiCallbacks().comm().rank() == 0) {
    auto so = std::make_shared<ParallelScriptInterface>("TestClass");
    /* Force destruction */
    so = nullptr;

    mpiCallbacks().abort_loop();
  } else {
    mpiCallbacks().loop();
  }

  Testing::reduce_and_check(mpiCallbacks().comm(), TestClass::constructed);
  Testing::reduce_and_check(mpiCallbacks().comm(), TestClass::destructed);
}

BOOST_AUTO_TEST_CASE(set_parmeter) {
  if (mpiCallbacks().comm().rank() == 0) {
    auto so = std::make_shared<ParallelScriptInterface>("TestClass");

    so->set_parameter("TestParam", std::string("TestValue"));

    mpiCallbacks().abort_loop();
    Testing::reduce_and_check(mpiCallbacks().comm(),
                              TestClass::last_instance != nullptr);

    auto const &last_parameter = TestClass::last_instance->last_parameter;
    Testing::reduce_and_check(mpiCallbacks().comm(),
                              last_parameter.first == "TestParam");
    Testing::reduce_and_check(mpiCallbacks().comm(),
                              boost::get<std::string>(last_parameter.second) ==
                                  "TestValue");
  } else {
    mpiCallbacks().loop();
    Testing::reduce_and_check(mpiCallbacks().comm(),
                              TestClass::last_instance != nullptr);
    auto const &last_parameter = TestClass::last_instance->last_parameter;
    Testing::reduce_and_check(mpiCallbacks().comm(),
                              last_parameter.first == "TestParam");
    Testing::reduce_and_check(mpiCallbacks().comm(),
                              boost::get<std::string>(last_parameter.second) ==
                                  "TestValue");
  }
}

/* TODO: test call_method */

int main(int argc, char **argv) {
  boost::mpi::environment mpi_env(argc, argv);
  ParallelScriptInterface::initialize();
  register_new<TestClass>("TestClass");

  boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
