#ifndef UTILS_MATH_SQR_HPP
#define UTILS_MATH_SQR_HPP

/** Calculates the SQuaRe of 'double' x, returning 'double'. */
template <typename T> inline T SQR(T x) { return x * x; }

namespace Utils {
  /** Calculates the SQuaRe of 'double' x, returning 'double'. */
  template <typename T> inline T sqr(T x) { return x * x; }

}

#endif
