#ifndef __ERROR_HANDLING_RUNTIME_ERROR_HPP
#define __ERROR_HANDLING_RUNTIME_ERROR_HPP

#include <sstream>
#include <string>

class RuntimeErrorCollector;

namespace ErrorHandling {

/** \brief Allows creating a runtime error messing by using the streaming operator */

class RuntimeErrorStream {
 public:
  RuntimeErrorStream(const RuntimeErrorStream &rhs);
  RuntimeErrorStream(RuntimeErrorCollector &ec, const std::string &file, const int line, const std::string &function);
  ~RuntimeErrorStream();
  template<typename T>
  RuntimeErrorStream &operator<<(T const &value) {
    m_buff << value;

    return *this;
  }
 private:
  RuntimeErrorCollector &m_ec;
  const int m_line;
  const std::string m_file, m_function;
  std::ostringstream m_buff;
};

}

#endif
