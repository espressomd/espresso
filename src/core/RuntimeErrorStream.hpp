#ifndef ERROR_HANDLING_RUNTIME_ERROR_STREAM_HPP
#define ERROR_HANDLING_RUNTIME_ERROR_STREAM_HPP

#include <sstream>
#include <string>

#include "RuntimeError.hpp"

namespace ErrorHandling {

class RuntimeErrorCollector;

/** \brief Allows creating a runtime error messing by using the streaming
 * operator */

class RuntimeErrorStream {
public:
  RuntimeErrorStream(RuntimeErrorStream &&) = default;
  RuntimeErrorStream(const RuntimeErrorStream &rhs);
  RuntimeErrorStream(RuntimeErrorCollector &ec, RuntimeError::ErrorLevel level,
                     const std::string &file, const int line,
                     const std::string &function);
  ~RuntimeErrorStream();
  template <typename T> RuntimeErrorStream &operator<<(T const &value) {
    m_buff << value;

    return *this;
  }

private:
  RuntimeErrorCollector &m_ec;
  RuntimeError::ErrorLevel m_level;
  const int m_line;
  const std::string m_file, m_function;
  std::ostringstream m_buff;
};

} /* ErrorHandling */

#endif
