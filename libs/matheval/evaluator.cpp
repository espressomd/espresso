#define MATHEVAL_IMPLEMENTATION

#include "evaluator.hpp"
#include "ast.hpp"

namespace matheval {

namespace ast {

// Optimizer

template <typename T> struct holds_alternative_impl {
    using result_type = bool;

    template <typename U> bool operator()(U const &) const {
        return std::is_same<U, T>::value;
    }
};

template <typename T, typename... Ts>
bool holds_alternative(x3::variant<Ts...> const &v) {
    return boost::apply_visitor(holds_alternative_impl<T>(), v);
}

ConstantFolder::result_type ConstantFolder::operator()(nil) const {
    return result_type{0};
}

ConstantFolder::result_type ConstantFolder::operator()(double n) const {
    return result_type{n};
}

ConstantFolder::result_type ConstantFolder::
operator()(std::string const &c) const {
    return result_type{c};
}

ConstantFolder::result_type ConstantFolder::
operator()(operation const &x, operand const &lhs) const {
    auto rhs = boost::apply_visitor(*this, x.rhs);

    if (holds_alternative<double>(lhs) && holds_alternative<double>(rhs)) {
        return result_type{
            x.op(boost::get<double>(lhs), boost::get<double>(rhs))};
    }
    return result_type{binary_op{x.op, lhs, rhs}};
}

ConstantFolder::result_type ConstantFolder::
operator()(unary_op const &x) const {
    auto rhs = boost::apply_visitor(*this, x.rhs);

    /// If the operand is known, we can directly evaluate the function.
    if (holds_alternative<double>(rhs)) {
        return result_type{x.op(boost::get<double>(rhs))};
    }
    return result_type{unary_op{x.op, rhs}};
}

ConstantFolder::result_type ConstantFolder::
operator()(binary_op const &x) const {
    auto lhs = boost::apply_visitor(*this, x.lhs);
    auto rhs = boost::apply_visitor(*this, x.rhs);

    /// If both operands are known, we can directly evaluate the function,
    /// else we just update the children with the new expressions.
    if (holds_alternative<double>(lhs) && holds_alternative<double>(rhs)) {
        return result_type{
            x.op(boost::get<double>(lhs), boost::get<double>(rhs))};
    }
    return result_type{binary_op{x.op, lhs, rhs}};
}

ConstantFolder::result_type ConstantFolder::
operator()(expression const &x) const {
    auto state = boost::apply_visitor(*this, x.lhs);
    for (operation const &oper : x.rhs) {
        state = (*this)(oper, state);
    }
    return result_type{state};
}

// Evaluator

double eval::operator()(nil) const {
    BOOST_ASSERT(0);
    return 0;
}

double eval::operator()(double n) const { return n; }

double eval::operator()(std::string const &c) const {
    auto it = st.find(c);
    if (it == st.end()) {
        throw std::invalid_argument("Unknown variable " + c); // NOLINT
    }
    return it->second;
}

double eval::operator()(operation const &x, double lhs) const {
    double rhs = boost::apply_visitor(*this, x.rhs);
    return x.op(lhs, rhs);
}

double eval::operator()(unary_op const &x) const {
    double rhs = boost::apply_visitor(*this, x.rhs);
    return x.op(rhs);
}

double eval::operator()(binary_op const &x) const {
    double lhs = boost::apply_visitor(*this, x.lhs);
    double rhs = boost::apply_visitor(*this, x.rhs);
    return x.op(lhs, rhs);
}

double eval::operator()(expression const &x) const {
    double state = boost::apply_visitor(*this, x.lhs);
    for (operation const &oper : x.rhs) {
        state = (*this)(oper, state);
    }
    return state;
}

} // namespace ast

} // namespace matheval
