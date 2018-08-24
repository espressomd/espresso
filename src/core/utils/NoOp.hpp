#ifndef CORE_UTILS_NO_OP_HPP
#define CORE_UTILS_NO_OP_HPP

namespace Utils {

/**
 * @brief A NoOp functor that does nothing.
 */
class NoOp {
public:
  template <typename... Args> void operator()(Args...) const { return; }
};

}
#endif
