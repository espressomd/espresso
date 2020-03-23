/** \file
 *  This file provides a wrapper around THRUST functions and types. In case
 *  that THRUST is not present, equivalent standard C++ types and functions
 *  are used.
 */

#pragma once



#ifdef __GNUC__
#define MAYBE_UNUSED __attribute__((unused))
#else
#define MAYBE_UNUSED
#endif





// Dependencies with THRUST
//#if __has_include(<thrust/device_vector.h>) // looks if the file exists
#if __has_include(<thrust/device_vectorBLAH.h>) // looks if the file exists // for testing non-thrust build
#  define SD_TRUST
#  include <thrust/device_vector.h>
#  include <thrust/execution_policy.h>
#  include <thrust/tabulate.h>
#  include <thrust/tuple.h>

namespace thrust_wrapper {
  // types
  using thrust::counting_iterator;
  using thrust::device_vector;
  using thrust::host_vector;
  using thrust::minus;
  using thrust::negate;
  using thrust::plus;
  using thrust::tuple;
  //using thrust::host;
  //using thrust::device;

  // routines
  using thrust::copy;
  using thrust::equal;
  using thrust::fill;
  using thrust::for_each;
  using thrust::get;
  using thrust::make_tuple;
  using thrust::raw_pointer_cast;
  using thrust::tabulate;
  using thrust::tie;
  using thrust::transform;
}


// vector addition
template <typename T>
thrust::device_vector<T> operator+(thrust::device_vector<T> const &x,
                                   thrust::device_vector<T> const &y) {
    assert(x.size() == y.size());
    thrust::device_vector<T> z(x.size());
    thrust::transform(x.begin(), x.end(), y.begin(), z.begin(),
                      thrust::plus<T>{});
    return z;
}

template <typename T>
thrust::host_vector<T> operator+(thrust::host_vector<T> const &x,
                                 thrust::host_vector<T> const &y) {
    assert(x.size() == y.size());
    thrust::host_vector<T> z(x.size());
    thrust::transform(x.begin(), x.end(), y.begin(), z.begin(),
                      thrust::plus<T>{});
    return z;
}








#else // Dependencies without THRUST

#  include <algorithm>
#  include <cassert>
#  include <iterator>
#  include <tuple>
#  include <vector>

#  include <boost/iterator/counting_iterator.hpp>


namespace thrust_wrapper {
  // types
  using boost::counting_iterator;
  using std::tuple;
  using std::minus;
  using std::negate;
  using std::plus;
  template<typename T> using device_vector = std::vector<T>;
  template<typename T> using host_vector   = std::vector<T>;
  // using host = void;
  // using device = void;


  // routines
  using std::for_each;
  using std::tie;
  using std::make_tuple;
  using std::get;
  using std::copy;
  // using std::fill; // uncomment if needed
  // using std::equal; // uncomment if needed


  template <typename T, typename ForwardIterator, typename UnaryOperation>
  void tabulate(MAYBE_UNUSED const T &exec,
                ForwardIterator first,
                ForwardIterator last,
                UnaryOperation unary_op) {
    while (first != last) {
      *first = unary_op(*first);
      ++first;
    }
  }

  // This wrapper function is needed because the THRUST variant gets the
  // execution policy as argument, which the STD library one doesn't.
  template <typename DerivedPolicy, typename InputIterator, typename UnaryFunction>
  void for_each(MAYBE_UNUSED const DerivedPolicy &exec,
                InputIterator first,
                InputIterator last,
                UnaryFunction f) {
    std::for_each(first, last, f);
  }

  // literally does nothing with its argument, since data isn't on a "device"
  template <typename T>
  T *raw_pointer_cast(T *ptr) {
    return ptr;
  }
  
  template <typename DerivedPolicy, typename ForwardIterator, typename T>
  void fill(MAYBE_UNUSED const DerivedPolicy &exec,
            ForwardIterator first,
            ForwardIterator last,
            const T &value) {
    std::fill(first, last, value);
  }

  // unary transform
  template <typename DerivedPolicy, typename ForwardIterator, 
            typename OutputIterator, typename UnaryFunction>
  void transform(MAYBE_UNUSED const DerivedPolicy &exec,
                 ForwardIterator first,
                 ForwardIterator last,
                 OutputIterator result,
                 UnaryFunction op) {
    std::transform(first, last, result, op);
  }

  // binary transform
  template <typename DerivedPolicy, typename InputIterator1, 
            typename InputIterator2, typename OutputIterator,
            typename BinaryFunction>
  bool transform(MAYBE_UNUSED const DerivedPolicy &exec,
                 InputIterator1 first1,
                 InputIterator1 last1,
                 InputIterator2 first2,
                 OutputIterator result,
                 BinaryFunction op) {
    return std::transform(first1, last1, first2, result, op);
  }
    
  template <typename DerivedPolicy, typename InputIterator1, 
            typename InputIterator2>
  bool equal(MAYBE_UNUSED const DerivedPolicy &exec,
             InputIterator1 first1,
             InputIterator1 last1,
             InputIterator2 first2) {
    return std::equal(first1, last1, first2);
  }

  // Wait, this following stuff has been way too complicated?!

  // Aliasing a template function with template argument deduction doesn't work
  //     https://www.fluentcpp.com/2017/10/27/function-aliases-cpp/
  //     Section "type deduction didn't follow"
  // Therefore we provide just the templates needed for Stokesian Dynamics.
  // If more are needed, add more, e.g. cases with three or more arguments.
  //template <typename T1, typename T2>
  //auto tie = std::template tie<T1, T2>;

  //template <typename T1, typename T2>
  //auto make_tuple = std::make_tuple<T1, T2>;

  //template <std::size_t I>
  //auto get = std::get<I>;

  //template <typename T, typename InputIterator, typename UnaryFunction>
  //void for_each(MAYBE_UNUSED const T &exec,
  //              InputIterator first,
  //              InputIterator last,
  //              UnaryFunction f) {
  //  std::for_each(first, last, f);
  //}
  
  // There are more parameters to this template, but Stokesian Dynamics is only
  // using this variant so far.
  //template<typename Incrementable>
  //using boost::counting_iterator<Incrementable>;

/*
  using device_vector = std::vector;

  // routines
  using std::tie;
*/
}


// vector addition
template <typename T>
std::vector<T> operator+(std::vector<T> const &x,
                         std::vector<T> const &y) {
    assert(x.size() == y.size());
    std::vector<T> z(x.size());
    std::transform(x.begin(), x.end(), y.begin(), z.begin(),
                   std::plus<T>{});
    return z;
}

#endif

#undef MAYBE_UNUSED
