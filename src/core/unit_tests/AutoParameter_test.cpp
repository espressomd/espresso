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

BOOST_AUTO_TEST_CASE(pointer_to_method) {
  using namespace ScriptInterface;
  struct C {
    C() = default;
    C(int i) : m_i(i) {}
    int m_i;

    void setter(int const &i) { m_i = i; }
    int &setter_getter() { return m_i; }
    int value_getter() const { return m_i; }
    int const &ref_getter() const { return m_i; }
  };

  {
    auto c_ptr = std::make_shared<C>();

    auto p = AutoParameter("name", c_ptr, &C::setter, &C::value_getter);
    p.set(5);
    BOOST_CHECK(5 == boost::get<int>(p.get()));
  }

  {
    auto c_ptr = std::make_shared<C>();

    auto p = AutoParameter("name", c_ptr, &C::setter, &C::value_getter);
    p.set(5);
    BOOST_CHECK(5 == boost::get<int>(p.get()));
  }

  {
    auto c_ptr = std::make_shared<C>();

    auto p = AutoParameter("name", c_ptr, &C::setter, &C::ref_getter);
    p.set(5);
    BOOST_CHECK(5 == boost::get<int>(p.get()));
  }

  {
    auto c_ptr = std::make_shared<C>();

    auto p_setter_getter = AutoParameter("name", c_ptr, &C::setter_getter);
    p_setter_getter.set(5);
    BOOST_CHECK(5 == boost::get<int>(p_setter_getter.get()));
  }

  {
    auto c_ptr = std::make_shared<C>(5);

    auto p = AutoParameter("name", c_ptr, &C::value_getter);
    BOOST_CHECK_THROW(p.set(5), AutoParameter::WriteError);
    BOOST_CHECK(5 == boost::get<int>(p.get()));
  }

  {
    auto c_ptr = std::make_shared<C>(5);

    auto p = AutoParameter("name", c_ptr, &C::ref_getter);
    BOOST_CHECK_THROW(p.set(5), AutoParameter::WriteError);
    BOOST_CHECK(5 == boost::get<int>(p.get()));
  }
}
