#ifndef ESPRESSO_DISTANCE_HPP
#define ESPRESSO_DISTANCE_HPP

/**
 * @brief Distance vector and length handed to pair kernels.
 */
struct Distance {
  explicit Distance(Utils::Vector3d const &vec21)
      : vec21(vec21), dist2(vec21.norm2()) {}

  Utils::Vector3d vec21;
  double dist2;
};

#endif // ESPRESSO_DISTANCE_HPP
