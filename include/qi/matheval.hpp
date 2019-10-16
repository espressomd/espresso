#pragma once

#include <map>
#include <memory>
#include <string>

namespace matheval {

/// @brief Parse a mathematical expression
///
/// This can parse and evaluate a mathematical expression for a given
/// symbol table using Boost.Spirit X3.  The templates of Boost.Spirit
/// are very expensive to parse and instantiate, which is why we hide
/// it behind an opaque pointer.
///
/// The drawback of this approach is that calls can no longer be
/// inlined and because the pointer crosses translation unit
/// boundaries, dereferencing it can also not be optimized out at
/// compile time.  We have to rely entirely on link-time optimization
/// which might be not as good.
///
/// The pointer to the implementation is a std::unique_ptr which makes
/// the class not copyable but only moveable.  Copying shouldn't be
/// required but is easy to implement.
class Parser {
    class impl;
    std::auto_ptr<impl> pimpl;

public:
    /// @brief Constructor
    Parser();

    /// @brief Destructor
    ~Parser();

    /// @brief Parse the mathematical expression into an abstract syntax tree
    ///
    /// @param[in] expr The expression given as a std::string
    void parse(std::string const &expr);

    /// @brief Perform constant folding onto the abstract syntax tree
    void optimize();

    /// @brief Evaluate the abstract syntax tree for a given symbol table
    ///
    /// @param[in] st The symbol table
    double evaluate(std::map<std::string, double> const &st);
};

/// @brief Convenience function
///
/// This function builds the grammar, parses the iterator to an AST,
/// evaluates it, and returns the result.
///
/// @param[in] expr  mathematical expression
/// @param[in] st    the symbol table for variables
inline double parse(std::string const &expr,
                    std::map<std::string, double> const &st) {
    Parser parser;
    parser.parse(expr);
    return parser.evaluate(st);
}

} // namespace matheval
