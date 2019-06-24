#define BOOST_TEST_MODULE tensor test
#define BOOST_TEST_DYN_LINK

#include <array>
#include <boost/test/unit_test.hpp>

#include <utils/tensor.hpp>

template <typename T> using Tensor = Utils::Tensor<T>;

BOOST_AUTO_TEST_CASE(ctor) {
  auto const extents = std::initializer_list<std::size_t>{4, 5, 6, 7};
  auto test_tensor = Tensor<double>(extents);
  for (int i = 0; i < 4; ++i) {
    BOOST_CHECK_EQUAL(*(std::begin(test_tensor.extents()) + i),
                      *(std::begin(extents) + i));
  }
}

BOOST_AUTO_TEST_CASE(ctor_from_container) {
  auto const extents = std::vector<std::size_t>{4, 5, 6, 7};
  auto test_tensor = Tensor<double>(extents);
  for (int i = 0; i < 4; ++i) {
    BOOST_CHECK_EQUAL(*(std::begin(test_tensor.extents()) + i),
                      *(std::begin(extents) + i));
  }
}

BOOST_AUTO_TEST_CASE(ctor_iterator) {
  auto const extents = std::array<std::size_t, 4>{4, 5, 6, 7};
  auto test_tensor = Tensor<double>(extents.begin(), extents.end());
  for (int i = 0; i < 4; ++i) {
    BOOST_CHECK_EQUAL(*(std::begin(test_tensor.extents()) + i),
                      *(std::begin(extents) + i));
  }
}

BOOST_AUTO_TEST_CASE(rank) {
  auto test_tensor = Tensor<double>({4, 6, 7, 8});
  BOOST_CHECK_EQUAL(test_tensor.rank(), 4);
}

BOOST_AUTO_TEST_CASE(data) {
  auto test_tensor = Tensor<double>({4, 6, 7, 8});
  test_tensor(0, 0, 0, 0) = 0.3;
  BOOST_CHECK_EQUAL(*test_tensor.data(), 0.3);
}

BOOST_AUTO_TEST_CASE(assign_initializer_list) {
  auto test_tensor = Tensor<double>({4, 6, 7, 8});
  test_tensor({0, 1, 3, 2}) = 2.3;
  BOOST_CHECK_EQUAL(test_tensor({0, 1, 3, 2}), 2.3);
}

BOOST_AUTO_TEST_CASE(assign_variadic) {
  auto test_tensor = Tensor<double>({4, 6, 7, 8});
  test_tensor(0, 1, 3, 2) = 2.3;
  BOOST_CHECK_EQUAL(test_tensor(0, 1, 3, 2), 2.3);
}

BOOST_AUTO_TEST_CASE(wrong_number_of_indices) {
  auto test_tensor = Tensor<double>({4, 6, 7, 8});
  BOOST_CHECK_THROW(test_tensor.at(0, 1, 3), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(out_of_bounds) {
  auto test_tensor = Tensor<double>({4, 6, 7, 8});
  BOOST_CHECK_THROW(test_tensor.at(4, 1, 3, 2), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(iterator) {
  auto test_tensor = Tensor<double>({4, 6, 7, 8});
  for (auto &t : test_tensor) {
    t = 1.4;
  }
  for (auto &t : test_tensor) {
    BOOST_CHECK_EQUAL(t, 1.4);
  }
}

BOOST_AUTO_TEST_CASE(front) {
  auto test_tensor = Tensor<double>({4, 6, 7, 8});
  test_tensor(0, 0, 0, 0) = 1.6;
  BOOST_CHECK_EQUAL(test_tensor.front(), 1.6);
}

BOOST_AUTO_TEST_CASE(back) {
  auto test_tensor = Tensor<double>({1, 3, 4, 5});
  test_tensor(0, 2, 3, 4) = 1.6;
  BOOST_CHECK_EQUAL(test_tensor.back(), 1.6);
}

BOOST_AUTO_TEST_CASE(row_major) {
  auto test_tensor = Tensor<double>({3, 4});
  test_tensor(1, 0) = 1.6;
  BOOST_CHECK_EQUAL(*(test_tensor.begin() + 1), 1.6);
}
