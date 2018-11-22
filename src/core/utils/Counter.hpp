#ifndef UTILS_COUNTER_HPP
#define UTILS_COUNTER_HPP

namespace Utils {
template <typename T> class Counter {
private:
  T m_val;
  T m_initial;

public:
  explicit Counter(T initial_value = T(0)) noexcept
      : m_val(initial_value), m_initial(initial_value) {}
  Counter(T initial_value, T value) noexcept
      : m_val(value), m_initial(initial_value) {}

  void increment() { ++m_val; }

  T value() const { return m_val; }
  T initial_value() const { return m_initial; }
};
} // namespace Utils

#endif // UTILS_COUNTER_HPP
