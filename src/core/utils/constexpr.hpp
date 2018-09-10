#ifndef UTILS_CONSTEXPR_HPP
#define UTILS_CONSTEXPR_HPP

/**
 * This exports the macro CXX14_CONSTEXPR
 * that is empty for c++11, and constexpr
 * for c++14, to be used in places that are
 * only allowed to be constexpr in c++14.
 */
#if __cplusplus < 201402L
#define CXX14_CONSTEXPR
#else
#define CXX14_CONSTEXPR constexpr
#endif
#endif
