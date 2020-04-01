/** \file
 *  This file provides a wrapper around THRUST functions and types. In case
 *  that CUDA/THRUST is not present, equivalent standard C++ types and 
 *  functions are used.
 */

#pragma once



// Dependencies with THRUST
#ifdef STOKESIAN_DYNAMICS_GPU
#  include <thrust/device_vector.h>
#  include <thrust/execution_policy.h>
#  include <thrust/tabulate.h>
#  include <thrust/tuple.h>

namespace thrust_wrapper {
  // types and constants
  using thrust::counting_iterator;
  using thrust::device_vector;
  using thrust::host_vector;
  using thrust::minus;
  using thrust::negate;
  using thrust::plus;
  using thrust::tuple;

  using thrust::host;
  using thrust::device;

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
  template<typename T> using device_vector = std::vector<T>;
  template<typename T> using host_vector   = std::vector<T>;
  using std::minus;
  using std::negate;
  using std::plus;
  using std::tuple;

  // dummy constants, they could be anything
  const int host = 0; 
  const int device = 0; 

  // routines
  using std::copy;
  using std::for_each;
  using std::get;
  using std::make_tuple;
  using std::tie;

  // tabulate doesn't exist in the standard library
  template <typename DerivedPolicy, typename ForwardIterator, 
            typename UnaryOperation>
  void tabulate(const DerivedPolicy &,
                ForwardIterator first,
                ForwardIterator last,
                UnaryOperation unary_op) {
    ForwardIterator i = first;
    while (i != last) {
      *i = unary_op(i - first);
      ++i;
    }
  }

  // literally does nothing with its argument
  template <typename T>
  T *raw_pointer_cast(T *ptr) {
    return ptr;
  }

  // These wrapper functions are needed because the THRUST variant gets the
  // execution policy as argument, which the STD library one doesn't.
  template <typename DerivedPolicy, typename InputIterator1, 
            typename InputIterator2>
  bool equal(const DerivedPolicy &,
             InputIterator1 first1,
             InputIterator1 last1,
             InputIterator2 first2) {
    return std::equal(first1, last1, first2);
  }

  template <typename DerivedPolicy, typename ForwardIterator, typename T>
  void fill(const DerivedPolicy &,
            ForwardIterator first,
            ForwardIterator last,
            const T &value) {
    std::fill(first, last, value);
  }

  template <typename DerivedPolicy, typename InputIterator, 
            typename UnaryFunction>
  void for_each(const DerivedPolicy &,
                InputIterator first,
                InputIterator last,
                UnaryFunction f) {
    std::for_each(first, last, f);
  }
  
  // unary transform
  template <typename DerivedPolicy, typename ForwardIterator, 
            typename OutputIterator, typename UnaryFunction>
  void transform(const DerivedPolicy &,
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
  bool transform(const DerivedPolicy &,
                 InputIterator1 first1,
                 InputIterator1 last1,
                 InputIterator2 first2,
                 OutputIterator result,
                 BinaryFunction op) {
    return std::transform(first1, last1, first2, result, op);
  }
    
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
