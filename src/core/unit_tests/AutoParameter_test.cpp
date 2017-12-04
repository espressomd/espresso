#define BOOST_TEST_MODULE AutoParameter test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "script_interface/auto_parameters/AutoParameter.hpp"

BOOST_AUTO_TEST_CASE(infer_length) {
  using ScriptInterface::infer_length;
  static_assert(infer_length<Vector<11, int>>() == 11, "");
  static_assert(infer_length<int>() == 0, "");
}

BOOST_AUTO_TEST_CASE(direct_binding) {
  using namespace ScriptInterface;
  int i{19};

  auto p = AutoParameter("i", i);

  BOOST_CHECK(p.type == VariantType::INT);

  BOOST_CHECK(boost::get<int>(p.get()) == 19);
  p.set(42);
  BOOST_CHECK(boost::get<int>(p.get()) == 42);
  BOOST_CHECK(i == 42);
}

BOOST_AUTO_TEST_CASE(method_pointer) {
  using namespace ScriptInterface;

  class A {
    int m_i;

  public:
    int i() const { return m_i; }
    void set_i(int const &i) { m_i = i; }
  };

  A a;

  auto p = AutoParameter("i", &a, &A::set_i, &A::i);

  BOOST_CHECK(p.type == VariantType::INT);

  p.set(42);
  BOOST_CHECK(boost::get<int>(p.get()) == 42);
}

BOOST_AUTO_TEST_CASE(read_only) {
  using namespace ScriptInterface;
  const int i = 12;

  auto p = AutoParameter("i", i);

  BOOST_CHECK(p.type == VariantType::INT);

  /* Getting should work as usual */
  BOOST_CHECK(boost::get<int>(p.get()) == i);

  /* Setting should throw */
  BOOST_CHECK_EXCEPTION(
      p.set(2), AutoParameter::WriteError,
      [](AutoParameter::WriteError const &e) { return true; });
}

BOOST_AUTO_TEST_CASE(user_provided) {
  using namespace ScriptInterface;
  int i{12};

  auto setter = [&i](Variant const &j) { i = boost::get<int>(j); };
  auto getter = [&i]() { return i; };

  auto p = AutoParameter("i", setter, getter);

  BOOST_CHECK(p.type == VariantType::INT);

  BOOST_CHECK(boost::get<int>(p.get()) == 12);
  p.set(42);
  BOOST_CHECK(boost::get<int>(p.get()) == 42);
  BOOST_CHECK(i == 42);
}
