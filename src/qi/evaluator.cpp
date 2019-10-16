#define MATHEVAL_IMPLEMENTATION

#include "evaluator.hpp"
#include "ast.hpp"

namespace matheval {

namespace ast {

// Optimizer

template <typename T, typename U>
struct is_same {
    static const bool value = false;
};

template <typename T>
struct is_same<T, T> {
    static const bool value = true;
};

template <typename T>
struct holds_alternative_impl {
    typedef bool result_type;

    template <typename U>
    bool operator()(U const & /*unused*/) const {
        return is_same<U, T>::value;
    }
};

template <typename T>
bool holds_alternative(operand const &v) {
    return boost::apply_visitor(holds_alternative_impl<T>(), v);
}

ConstantFolder::result_type ConstantFolder::operator()(nil) const { return 0; }

ConstantFolder::result_type ConstantFolder::operator()(double n) const {
    return n;
}

ConstantFolder::result_type ConstantFolder::
operator()(std::string const &c) const {
    return c;
}

ConstantFolder::result_type ConstantFolder::
operator()(operation const &x, operand const &lhs) const {
    operand rhs = boost::apply_visitor(*this, x.rhs);

    if (holds_alternative<double>(lhs) && holds_alternative<double>(rhs)) {
        return x.op(boost::get<double>(lhs), boost::get<double>(rhs));
    }
    return binary_op(x.op, lhs, rhs);
}

ConstantFolder::result_type ConstantFolder::
operator()(unary_op const &x) const {
    operand rhs = boost::apply_visitor(*this, x.rhs);

    /// If the operand is known, we can directly evaluate the function.
    if (holds_alternative<double>(rhs)) {
        return x.op(boost::get<double>(rhs));
    }
    return unary_op(x.op, rhs);
}

ConstantFolder::result_type ConstantFolder::
operator()(binary_op const &x) const {
    operand lhs = boost::apply_visitor(*this, x.lhs);
    operand rhs = boost::apply_visitor(*this, x.rhs);

    /// If both operands are known, we can directly evaluate the function,
    /// else we just update the children with the new expressions.
    if (holds_alternative<double>(lhs) && holds_alternative<double>(rhs)) {
        return x.op(boost::get<double>(lhs), boost::get<double>(rhs));
    }
    return binary_op(x.op, lhs, rhs);
}

ConstantFolder::result_type ConstantFolder::
operator()(expression const &x) const {
    operand state = boost::apply_visitor(*this, x.lhs);
    for (std::list<operation>::const_iterator it = x.rhs.begin();
         it != x.rhs.end(); ++it) {
        state = (*this)(*it, state);
    }
    return state;
}

// Evaluator

double eval::operator()(nil) const {
    BOOST_ASSERT(0);
    return 0;
}

double eval::operator()(double n) const { return n; }

double eval::operator()(std::string const &c) const {
    std::map<std::string, double>::const_iterator it = st.find(c);
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
    for (std::list<operation>::const_iterator it = x.rhs.begin();
         it != x.rhs.end(); ++it) {
        state = (*this)(*it, state);
    }
    return state;
}

} // namespace ast

} // namespace matheval
