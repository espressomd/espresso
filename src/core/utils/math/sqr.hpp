#ifndef UTILS_MATH_SQR_HPP
#define UTILS_MATH_SQR_HPP

namespace Utils {
/** Calculates the SQuaRe of 'double' x, returning 'double'. */
template <typename T> inline T sqr(T x) { return x * x; }

} // namespace Utils

#endif
