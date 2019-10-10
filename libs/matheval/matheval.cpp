#define MATHEVAL_IMPLEMENTATION

#include "matheval.hpp"
#include "ast.hpp"
#include "evaluator.hpp"
#include "parser.hpp"

#include <memory>
#include <stdexcept>
#include <string>

namespace matheval {

class Parser::impl {
    ast::operand ast;

public:
    void parse(std::string const &expr) {
        auto ast_ = ast::expression{};

        auto first = expr.begin();
        auto last = expr.end();

        boost::spirit::x3::ascii::space_type space;
        bool r = phrase_parse(first, last, grammar(), space, ast_);

        if (!r || first != last) {
            std::string rest(first, last);
            throw std::runtime_error("Parsing failed at " + rest); // NOLINT
        }

        ast = ast_;
    }

    void optimize() { ast = boost::apply_visitor(ast::ConstantFolder{}, ast); }

    double evaluate(std::map<std::string, double> const &st) {
        return boost::apply_visitor(ast::eval{st}, ast);
    }
};

Parser::Parser() : pimpl{std::make_unique<Parser::impl>()} {}

Parser::~Parser() {}

void Parser::parse(std::string const &expr) { pimpl->parse(expr); }

void Parser::optimize() { pimpl->optimize(); }

double Parser::evaluate(std::map<std::string, double> const &st) {
    return pimpl->evaluate(st);
}

} // namespace matheval
