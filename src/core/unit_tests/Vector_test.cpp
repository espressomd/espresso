
#include <iostream>
#include "../Vector.hpp"

#define BOOST_TEST_MODULE Vector test
#include <boost/test/included/unit_test.hpp>

/** Number of nontrivial Baxter permutations of length 2n-1. (A001185) */
#define TEST_NUMBERS { 0, 1, 1, 7, 21, 112, 456, 2603, 13203 }
#define TEST_NUMBERS_PARTIAL_NORM2 { 0, 1, 2, 51, 492, 13036 }
const int test_numbers[] = TEST_NUMBERS;
const int test_numbers_partial_norm2[] = TEST_NUMBERS_PARTIAL_NORM2;
const int n_test_numbers = sizeof(test_numbers) / sizeof(int);

template <int n>
bool il_constructor() {
  bool pass = true;
  Vector<n, int> v(TEST_NUMBERS);

  pass &= v.size() == n; 
  
  for(int i = 0; i < std::min(n_test_numbers, n); i++)
    pass &= v[i] == test_numbers[i];
    
  return pass;
}

template<int n>
bool default_constructor() {
  bool pass = true;
  Vector<n, int> v;

  for(int i = 0; i < n; i++) {
    v[i] = i;
  }
  
  return pass;
}

template<int n>
bool norm2() {
  Vector<n, int> v(TEST_NUMBERS);

  std::cout << v.norm2() << " == " << test_numbers_partial_norm2[n] << std::endl;
  
  return v.norm2() == test_numbers_partial_norm2[n-1];
}

BOOST_AUTO_TEST_CASE( test_constructor )
{
  BOOST_CHECK(il_constructor<1>());
  BOOST_CHECK(il_constructor<2>());
  BOOST_CHECK(il_constructor<3>());
  BOOST_CHECK(il_constructor<4>());
  BOOST_CHECK(il_constructor<5>());
  BOOST_CHECK(il_constructor<6>());
  BOOST_CHECK(il_constructor<7>());
  BOOST_CHECK(il_constructor<8>());
  BOOST_CHECK(il_constructor<9>());
  BOOST_CHECK(il_constructor<10>());

  BOOST_CHECK(default_constructor<1>());
  BOOST_CHECK(default_constructor<2>());
  BOOST_CHECK(default_constructor<3>());
  BOOST_CHECK(default_constructor<4>());
  BOOST_CHECK(default_constructor<5>());
  BOOST_CHECK(default_constructor<6>());
  BOOST_CHECK(default_constructor<7>());
  BOOST_CHECK(default_constructor<8>());
  BOOST_CHECK(default_constructor<9>());
  BOOST_CHECK(default_constructor<10>());  
}

BOOST_AUTO_TEST_CASE( test_norm2 ) {
  BOOST_CHECK(norm2<1>());
  BOOST_CHECK(norm2<2>());
  BOOST_CHECK(norm2<3>());
  BOOST_CHECK(norm2<4>());
  BOOST_CHECK(norm2<5>());
  BOOST_CHECK(norm2<6>()); 
}
