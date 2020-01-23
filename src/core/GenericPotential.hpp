#ifndef CORE_GENERIC_POTENTIAL_HPP
#define CORE_GENERIC_POTENTIAL_HPP

#include "config.hpp"
#ifdef MATHEVAL

#include <cassert>
#include <memory>
#include <string>

#include <boost/serialization/access.hpp>
#include <boost/serialization/string.hpp>

#include "integrate.hpp"
#include "matheval.hpp"

struct GenericPotential {
  double maxval = -1.0;
  std::string force_expr;
  std::string energy_expr;
  std::shared_ptr<matheval::Parser> force_parser{nullptr};
  std::shared_ptr<matheval::Parser> energy_parser{nullptr};
  bool is_parsed{false};

  void parse() {
    if (!force_parser || !energy_parser) {
      throw std::runtime_error("nullptr dereference");
    }

    if (!force_expr.empty() && !energy_expr.empty()) {
      force_parser->parse(force_expr);
      force_parser->optimize();
      energy_parser->parse(energy_expr);
      energy_parser->optimize();

      is_parsed = true;
    }
  }

  double force(double x) const {
    assert(x <= maxval);
    assert(is_parsed);
    return force_parser->evaluate(
        {std::make_pair("x", x), std::make_pair("t", sim_time)});
  }

  double energy(double x) const {
    assert(x <= maxval);
    assert(is_parsed);
    return energy_parser->evaluate(
        {std::make_pair("x", x), std::make_pair("t", sim_time)});
  }

  double cutoff() const { return maxval; }

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, long int /* version */) {
    ar &maxval;
    ar &force_expr;
    ar &energy_expr;
  }
};

#endif
#endif
