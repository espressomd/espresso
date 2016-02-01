#include "RuntimeErrorStream.hpp"
#include "RuntimeErrorCollector.hpp"

namespace ErrorHandling {
/** ostringstream is not copyable, but it is fine here to copy just the content. */
RuntimeErrorStream::RuntimeErrorStream(const RuntimeErrorStream &rhs) : m_ec(rhs.m_ec), m_line(rhs.m_line), m_file(rhs.m_file), m_function(rhs.m_function) {
  m_buff << rhs.m_buff.rdbuf();
}


RuntimeErrorStream::RuntimeErrorStream(RuntimeErrorCollector &ec, const std::string &file, const int line, const std::string &function) :
    m_ec(ec), m_line(line), m_file(file), m_function(function) {}
RuntimeErrorStream::~RuntimeErrorStream() {
  m_ec.error(m_buff, m_function.c_str(), m_file.c_str(), m_line);
}

}
