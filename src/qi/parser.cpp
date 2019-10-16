#define MATHEVAL_IMPLEMENTATION

#include "parser_def.hpp"

namespace matheval {

typedef std::string::const_iterator iterator_type;
template struct parser::grammar<iterator_type>;

} // namespace matheval
