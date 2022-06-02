#ifndef SCRIPT_INTERFACE_WALBERLA_NONE_HPP
#define SCRIPT_INTERFACE_WALBERLA_NONE_HPP

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameter.hpp"
#include "script_interface/walberla/EKPoissonSolver.hpp"

#include "walberla_bridge/electrokinetics/ek_walberla_init.hpp"

#include <memory>

namespace ScriptInterface::walberla {

class EKNone : public EKPoissonSolver {
public:
  void do_construct(VariantMap const &args) override {
    m_single_precision = get_value_or<bool>(args, "single_precision", false);
    auto lattice =
        get_value<std::shared_ptr<LatticeWalberla>>(args, "lattice")->lattice();

    m_noneinstance = new_ek_poisson_none(lattice, m_single_precision);
  }

  [[nodiscard]] std::shared_ptr<::walberla::PoissonSolver> get_instance() const
      noexcept override {
    return m_noneinstance;
  }

private:
  /* The actual instance */
  std::shared_ptr<::walberla::PoissonSolver> m_noneinstance;

  bool m_single_precision;
};
} // namespace ScriptInterface::walberla

#endif // SCRIPT_INTERFACE_WALBERLA_NONE_HPP
