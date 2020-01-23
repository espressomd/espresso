#define MATHEVAL_IMPLEMENTATION

#include "parser_def.hpp"

namespace matheval {

template struct parser::grammar<std::string::const_iterator>;

} // namespace matheval
