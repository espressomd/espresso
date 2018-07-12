// Copyright (C) 2011-2013 Rhys Ulerich
// Copyright (C) ??? Martin Bauer
// Copyright (C) 2017 Henri Menke
//
// This code borrows heavily from code written by Rhys Ulerich and
// Martin Bauer.  They licensed it under the Mozilla Public License,
// v. 2.0 and the GNU General Public License (no version info),
// respectively.  I believe that I have made enough contributions and
// altered this code far enough from the originals that I can
// relicense it under the Boost Software License.
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#pragma once

#if defined __clang__
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wunsequenced"
#elif defined __GNUC__
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wsequence-point"
#endif

#include <algorithm>
#include <cmath>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

#define BOOST_RESULT_OF_USE_DECLTYPE
#define BOOST_SPIRIT_USE_PHOENIX_V3
#include <boost/math/constants/constants.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/variant.hpp>


namespace matheval {

namespace detail {

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

namespace math {

/** @brief Sign function
 *
 * Missing function in the STL.  This calculates the mathematical sign
 * function.
 *
 * @f[
 *   \mathop{\mathrm{sgn}}(x) =
 *   \begin{cases}
 *      1 & x > 0 \\
 *      0 & x = 0 \\
 *     -1 & x < 0 \\
 *   \end{cases}
 * @f]
 *
 * @param[in] x number
 * @returns the sign of x
 */
template < typename T >
T sgn(T x) { return (T{0} < x) - (x < T{0}); }

/** @brief isnan function with adjusted return type */
template < typename T >
T isnan(T x) { return std::isnan(x); }

/** @brief isinf function with adjusted return type */
template < typename T >
T isinf(T x) { return std::isinf(x); }

/** @brief Convert radians to degrees */
template < typename T >
T deg(T x) { return x*boost::math::constants::radian<T>(); }

/** @brief Convert degrees to radians */
template < typename T >
T rad(T x) { return x*boost::math::constants::degree<T>(); }

}

// AST

template < typename real_t > struct unary_op;
template < typename real_t > struct binary_op;

struct nil {};

/** @brief Abstract Syntax Tree
 *
 * Stores the abstract syntax tree (AST) of the parsed mathematical
 * expression.
 */
template < typename real_t >
struct expr_ast
{
    using tree_t = boost::variant<
        nil // can't happen!
        , real_t
        , std::string
        , boost::recursive_wrapper<expr_ast<real_t>>
        , boost::recursive_wrapper<binary_op<real_t>>
        , boost::recursive_wrapper<unary_op<real_t>>
        >;
public:
    /** @brief AST storage
     *
     * The syntax tree can hold various types.  Numbers (`real_t`),
     * variables (`std::string`), the recursive tree itself
     * (`expr_ast`), binary operators (`binary_op`), and unary
     * operators (`unary_op`).
     */
    tree_t tree;

    /** @brief Default constructor
     *
     * Initializes the tree to a nil value to indicate inconsistent
     * state.
     */
    expr_ast() : tree(nil{}) {}

    /** @brief Copy constructor
     *
     * Deep copies the syntax tree.
     */
    template <typename Expr>
    expr_ast(Expr other) : tree(std::move(other)) {} // NOLINT

    /** @brief Add a tree */
    expr_ast& operator+=(expr_ast const &rhs);
    /** @brief subtract a tree */
    expr_ast& operator-=(expr_ast const &rhs);
    /** @brief Multiply by a tree */
    expr_ast& operator*=(expr_ast const &rhs);
    /** @brief Divide by a tree */
    expr_ast& operator/=(expr_ast const &rhs);
};

/** @brief Store a unary operator and its argument tree */
template < typename real_t >
struct unary_op
{
    /** @brief Signature of a unary operator: op(x) */
    using op_t = std::function<real_t(real_t)>;

    /** @brief Save the operator and the argument tree */
    unary_op(op_t op, expr_ast<real_t> rhs)
        : op(std::move(op)), rhs(std::move(rhs))
    {}

    /** @brief Stored operator */
    op_t op;
    /** @brief Stored argument tree */
    expr_ast<real_t> rhs;
};

/** @brief Store a binary operator and its argument trees */
template < typename real_t >
struct binary_op
{
    /** @brief Signature of a binary operator: op(x,y) */
    using op_t = std::function<real_t(real_t,real_t)>;

    /** @brief Save the operator and the argument trees */
    binary_op(op_t op, expr_ast<real_t> lhs, expr_ast<real_t> rhs)
        : op(std::move(op)), lhs(std::move(lhs)), rhs(std::move(rhs))
    {}

    /** @brief Stored operator */
    op_t op;
    /** @brief Stored argument tree of first argument */
    expr_ast<real_t> lhs;
    /** @brief Stored argument tree of second argument */
    expr_ast<real_t> rhs;
};

template < typename real_t >
expr_ast<real_t>& expr_ast<real_t>::operator+=(expr_ast<real_t> const &rhs)
{
    tree = binary_op<real_t>(std::plus<real_t>{}, tree, rhs);
    return *this;
}
template < typename real_t >
expr_ast<real_t>& expr_ast<real_t>::operator-=(expr_ast<real_t> const &rhs)
{
    tree = binary_op<real_t>(std::minus<real_t>{}, tree, rhs);
    return *this;
}
template < typename real_t >
expr_ast<real_t>& expr_ast<real_t>::operator*=(expr_ast<real_t> const &rhs)
{
    tree = binary_op<real_t>(std::multiplies<real_t>{}, tree, rhs);
    return *this;
}
template < typename real_t >
expr_ast<real_t>& expr_ast<real_t>::operator/=(expr_ast<real_t> const &rhs)
{
    tree = binary_op<real_t>(std::divides<real_t>{}, tree, rhs);
    return *this;
}

/** @brief Evaluate the Abstract Syntax Tree
 *
 * This visits all the variants of the AST and applies the stored
 * operators.
 */
template < typename real_t >
class eval_ast
{
public:
    /** @brief Necessary typedef for `boost::apply_visitor` */
    using result_type = real_t;

    /** @brief Type of the symbol table */
    using symbol_table_t = std::map<std::string, result_type>;

    /** @brief Constructor
     *
     * Saves the symbol table to apply variables.
     */
    explicit eval_ast(symbol_table_t sym) : st(std::move(sym)) {}

    /** @brief Empty nodes in the tree evaluate to 0 */
    result_type operator()(nil /*unused*/) const { return 0; }

    /** @brief Numbers evaluate to themselves */
    result_type operator()(result_type n)  const { return n; }

    /** @brief Variables evaluate to their value in the symbol table */
    result_type operator()(std::string const &c) const
    {
        auto it = st.find(c);
        if (it == st.end()) {
            throw std::invalid_argument("Unknown variable " + c); // NOLINT
        }
        return it->second;
    }

    /** @brief Recursively evaluate the AST */
    result_type operator()(expr_ast<real_t> const& ast) const
    {
        return boost::apply_visitor(*this, ast.tree);
    }

    /** @brief Evaluate a binary operator and optionally recurse its operands */
    result_type operator()(binary_op<real_t> const& tree) const
    {
        return tree.op(
            boost::apply_visitor(*this, tree.lhs.tree),
            boost::apply_visitor(*this, tree.rhs.tree)
            );
    }

    /** @brief Evaluate a unary operator and optionally recurse its operand */
    result_type operator()(unary_op<real_t> const& tree) const
    {
        return tree.op(
            boost::apply_visitor(*this, tree.rhs.tree)
            );
    }

private:
    symbol_table_t st;
};

template <typename T> struct holds_alternative_impl {
  using result_type = bool;

  template <typename U> bool operator()(U const & /*unused*/) const {
    return std::is_same<U, T>::value;
  }
};

template <typename T, typename... Ts>
bool holds_alternative(boost::variant<Ts...> const &v) {
  return boost::apply_visitor(holds_alternative_impl<T>(), v);
}

template <typename real_t> struct ConstantFolder {
  /** @brief Necessary typedef for `boost::apply_visitor` */
  using result_type = typename expr_ast<real_t>::tree_t;

  /** @brief Empty nodes in the tree evaluate to 0 */
  result_type operator()(nil /*unused*/) const { return 0; }

  /** @brief Numbers evaluate to themselves */
  result_type operator()(real_t n) const { return n; }

  /** @brief Variables do not evaluate */
  result_type operator()(std::string const &c) const { return c; }

  /** @brief Recursively evaluate the AST */
  result_type operator()(expr_ast<real_t> const &ast) const {
    return boost::apply_visitor(*this, ast.tree);
  }

  /** @brief Evaluate a binary operator and optionally recurse its operands */
  result_type operator()(binary_op<real_t> const &tree) const {
    auto lhs = boost::apply_visitor(*this, tree.lhs.tree);
    auto rhs = boost::apply_visitor(*this, tree.rhs.tree);

    /* If both operands are known, we can directly evaluate the function,
     * else we just update the children with the new expressions. */
    if (holds_alternative<real_t>(lhs) && holds_alternative<real_t>(rhs)) {
      return tree.op(boost::get<real_t>(lhs), boost::get<real_t>(rhs));
    }
    return binary_op<real_t>(tree.op, lhs, rhs);
  }

  /** @brief Evaluate a unary operator and optionally recurse its operand */
  result_type operator()(unary_op<real_t> const &tree) const {
    auto rhs = boost::apply_visitor(*this, tree.rhs.tree);
    /* If the operand is known, we can directly evaluate the function. */
    if (holds_alternative<real_t>(rhs)) {
      return tree.op(boost::get<real_t>(rhs));
    }
    return unary_op<real_t>(tree.op, rhs);
  }
};

// Expressions

/** @brief Unary expression functor */
template < typename real_t >
struct unary_expr_ {
    /** @brief Make boost::phoenix::function happy */
    template < typename T > struct result { using type = T; };

    /** @brief Create a new AST containing the unary function */
    expr_ast<real_t> operator()(typename unary_op<real_t>::op_t op,
                                expr_ast<real_t> const &rhs) const {
        return expr_ast<real_t>(unary_op<real_t>(op, rhs));
    }
};

/** @brief Binary expression functor */
template < typename real_t >
struct binary_expr_ {
    /** @brief Make boost::phoenix::function happy */
    template < typename T > struct result { using type = T; };

    /** @brief Create a new AST containing the binary function */
    expr_ast<real_t> operator()(typename binary_op<real_t>::op_t op,
                                expr_ast<real_t> const &lhs,
                                expr_ast<real_t> const &rhs) const {
        return expr_ast<real_t>(binary_op<real_t>(op, lhs, rhs));
    }
};

/** @brief Error handler for expectation errors */
struct expectation_handler {
    /** @brief Throw an exception saying where and why parsing failed */
    template < typename Iterator >
    void operator()(Iterator first, Iterator last,
                    boost::spirit::info const& info) const {
       std::stringstream msg;
       msg << "Expected "
           << info
           << " at \""
           << std::string{first, last}
           << "\"";

       throw std::runtime_error(msg.str()); // NOLINT
    }
};


// Grammar

/** @brief Expression Grammar */
template < typename real_t, typename Iterator >
struct grammar
    : qi::grammar<
            Iterator, expr_ast<real_t>(), ascii::space_type
        >
{
private:
    expectation_handler err_handler;
    qi::rule<Iterator, expr_ast<real_t>(), ascii::space_type> expression;
    qi::rule<Iterator, expr_ast<real_t>(), ascii::space_type> term;
    qi::rule<Iterator, expr_ast<real_t>(), ascii::space_type> factor;
    qi::rule<Iterator, expr_ast<real_t>(), ascii::space_type> primary;
    qi::rule<Iterator, std::string()> variable;
public:
    /** @brief symbol table for constants like "pi" */
    struct constant_
        : boost::spirit::qi::symbols<
                typename std::iterator_traits<Iterator>::value_type,
                real_t
            >
    {
        constant_()
        {
            this->add
                ("e"      , boost::math::constants::e<real_t>()   )
                ("epsilon", std::numeric_limits<real_t>::epsilon())
                ("phi"    , boost::math::constants::phi<real_t>() )
                ("pi"     , boost::math::constants::pi<real_t>()  )
            ;
        }
    } constant;

    /** @brief symbol table for unary functions like "abs" */
    struct ufunc_
        : boost::spirit::qi::symbols<
                typename std::iterator_traits<Iterator>::value_type,
                typename unary_op<real_t>::op_t
            >
    {
        ufunc_()
        {
            this->add
                ("abs"   , static_cast<real_t(*)(real_t)>(&std::abs   ))
                ("acos"  , static_cast<real_t(*)(real_t)>(&std::acos  ))
                ("acosh" , static_cast<real_t(*)(real_t)>(&std::acosh ))
                ("asin"  , static_cast<real_t(*)(real_t)>(&std::asin  ))
                ("asinh" , static_cast<real_t(*)(real_t)>(&std::asinh ))
                ("atan"  , static_cast<real_t(*)(real_t)>(&std::atan  ))
                ("atanh" , static_cast<real_t(*)(real_t)>(&std::atanh ))
                ("cbrt"  , static_cast<real_t(*)(real_t)>(&std::cbrt  ))
                ("ceil"  , static_cast<real_t(*)(real_t)>(&std::ceil  ))
                ("cos"   , static_cast<real_t(*)(real_t)>(&std::cos   ))
                ("cosh"  , static_cast<real_t(*)(real_t)>(&std::cosh  ))
                ("deg"   , static_cast<real_t(*)(real_t)>(&math::deg  ))
                ("erf"   , static_cast<real_t(*)(real_t)>(&std::erf   ))
                ("erfc"  , static_cast<real_t(*)(real_t)>(&std::erfc  ))
                ("exp"   , static_cast<real_t(*)(real_t)>(&std::exp   ))
                ("exp2"  , static_cast<real_t(*)(real_t)>(&std::exp2  ))
                ("floor" , static_cast<real_t(*)(real_t)>(&std::floor ))
                ("isinf" , static_cast<real_t(*)(real_t)>(&math::isinf))
                ("isnan" , static_cast<real_t(*)(real_t)>(&math::isnan))
                ("log"   , static_cast<real_t(*)(real_t)>(&std::log   ))
                ("log2"  , static_cast<real_t(*)(real_t)>(&std::log2  ))
                ("log10" , static_cast<real_t(*)(real_t)>(&std::log10 ))
                ("rad"   , static_cast<real_t(*)(real_t)>(&math::rad  ))
                ("round" , static_cast<real_t(*)(real_t)>(&std::round ))
                ("sgn"   , static_cast<real_t(*)(real_t)>(&math::sgn  ))
                ("sin"   , static_cast<real_t(*)(real_t)>(&std::sin   ))
                ("sinh"  , static_cast<real_t(*)(real_t)>(&std::sinh  ))
                ("sqrt"  , static_cast<real_t(*)(real_t)>(&std::sqrt  ))
                ("tan"   , static_cast<real_t(*)(real_t)>(&std::tan   ))
                ("tanh"  , static_cast<real_t(*)(real_t)>(&std::tanh  ))
                ("tgamma", static_cast<real_t(*)(real_t)>(&std::tgamma))
            ;
        }
    } ufunc;

    /** @brief symbol table for binary functions like "pow" */
    struct bfunc_
        : boost::spirit::qi::symbols<
                typename std::iterator_traits<Iterator>::value_type,
                typename binary_op<real_t>::op_t
            >
    {
        bfunc_()
        {
            this->add
                ("atan2", static_cast<real_t(*)(real_t,real_t)>(&std::atan2))
                ("max"  , static_cast<real_t(*)(real_t,real_t)>(&std::fmax ))
                ("min"  , static_cast<real_t(*)(real_t,real_t)>(&std::fmin ))
                ("pow"  , static_cast<real_t(*)(real_t,real_t)>(&std::pow  ))
            ;
        }
    } bfunc;

    /** @brief Constructor builds the grammar */
    grammar() : grammar::base_type(expression)
    {
        using boost::spirit::qi::real_parser;
        using boost::spirit::qi::real_policies;
        real_parser<real_t,real_policies<real_t>> real;

        using boost::spirit::lexeme;
        using boost::spirit::qi::_1;
        using boost::spirit::qi::_2;
        using boost::spirit::qi::_3;
        using boost::spirit::qi::_4;
        using boost::spirit::qi::_val;
        using boost::spirit::qi::alpha;
        using boost::spirit::qi::alnum;
        using boost::spirit::qi::raw;

        boost::phoenix::function<unary_expr_<real_t>> unary_expr;
        boost::phoenix::function<binary_expr_<real_t>> binary_expr;

        auto fmod = static_cast<real_t(*)(real_t,real_t)>(&std::fmod);
        auto pow = static_cast<real_t(*)(real_t,real_t)>(&std::pow);

        expression =
            term                  [_val =  _1]
            >> *(  ('+' > term    [_val += _1])
                |  ('-' > term    [_val -= _1])
                )
            ;

        term =
            factor                [_val =  _1]
            >> *(  ('*' > factor  [_val *= _1])
                |  ('/' > factor  [_val /= _1])
                |  ('%' > factor  [_val = binary_expr(fmod, _val, _1)])
                )
            ;

        factor =
            primary               [_val =  _1]
            >> *(  ("**" > factor [_val = binary_expr(pow, _val, _1)])
                )
            ;

        variable =
            raw[lexeme[alpha >> *(alnum | '_')]];

        primary =
            real                   [_val =  _1]
            |   ('(' > expression  [_val =  _1] > ')')
            |   ('-' > primary     [_val = unary_expr(std::negate<real_t>{}, _1)])
            |   ('+' > primary     [_val =  _1])
            |   (bfunc > '(' > expression > ',' > expression > ')')
                                   [_val = binary_expr(_1, _2, _3)]
            |   (ufunc > '(' > expression > ')')
                                   [_val = unary_expr(_1, _2)]
            |   constant           [_val =  _1]
            |   variable           [_val =  _1]
            ;

        expression.name("expression");
        term.name("term");
        factor.name("factor");
        variable.name("variable");
        primary.name("primary");

        using boost::spirit::qi::fail;
        using boost::spirit::qi::on_error;
        using boost::phoenix::bind;
        using boost::phoenix::ref;

        on_error<fail>
        (
            expression,
            bind(ref(err_handler), _3, _2, _4)
        );
    }
};

/** @brief Parse an expression
 *
 * This function builds the grammar and parses the iterator into
 * an AST.
 *
 * @param[in] first iterator to the start of the input sequence
 * @param[in] last  iterator to the end of the input sequence
 */
template <typename real_t, typename Iterator>
detail::expr_ast<real_t> parse(Iterator first, Iterator last) {
  static detail::grammar<real_t, Iterator> const g;

  auto ast = detail::expr_ast<real_t>{};

  bool r = qi::phrase_parse(first, last, g, ascii::space, ast);

  if (!r || first != last) {
    std::string rest(first, last);
    throw std::runtime_error("Parsing failed at " + rest); // NOLINT
  }

  return ast;
}

} // namespace detail


/** @brief Class interface
 *
 * This class hides the grammar, AST, and AST traversal behind some
 * member functions.
 *
 * @tparam real_t datatype of the result
 */
template < typename real_t >
class Parser
{
    detail::expr_ast<real_t> ast;
public:
    /** @brief Parse an expression
     *
     * This function builds the grammar and parses the iterator into
     * an AST.
     *
     * @param[in] first iterator to the start of the input sequence
     * @param[in] last  iterator to the end of the input sequence
     */
    template < typename Iterator >
    void parse(Iterator first, Iterator last)
    {
        ast = detail::parse<real_t>(first, last);
    }

    /** @overload parse(Iterator first, Iterator last) */
    void parse(std::string const &str)
    {
        parse(str.begin(), str.end());
    }

    void optimize() {
      ast.tree = boost::apply_visitor(detail::ConstantFolder<real_t>{}, ast.tree);
    }

    /** @brief Evaluate the AST with a given symbol table
     *
     * @param[in] st the symbol table for variables
     */
    real_t evaluate(typename detail::eval_ast<real_t>::symbol_table_t const &st)
    {
        detail::eval_ast<real_t> solver(st);
        return solver(ast);
    }
};


/** @brief Convenience function
 *
 * This function builds the grammar, parses the iterator to an AST,
 * evaluates it, and returns the result.
 *
 * @param[in] first iterator to the start of the input sequence
 * @param[in] last  iterator to the end of the input sequence
 * @param[in] st    the symbol table for variables
 */
template < typename real_t, typename Iterator >
real_t parse(Iterator first, Iterator last,
             typename detail::eval_ast<real_t>::symbol_table_t const &st)
{
    Parser<real_t> parser;
    parser.parse(first, last);
    return parser.evaluate(st);
}

/** @overload parse(Iterator first, Iterator last, typename detail::eval_ast<real_t>::symbol_table_t const &st) */
template < typename real_t >
real_t parse(std::string const &str,
             typename detail::eval_ast<real_t>::symbol_table_t const &st)
{
    return parse<real_t>(str.begin(), str.end(), st);
}

} // namespace expression

#if defined __clang__
#  pragma clang diagnostic pop
#elif defined __GNUC__
#  pragma GCC diagnostic pop
#endif
