#ifndef SCRIPT_INTERFACE_WALBERLA_NONE_HPP
#define SCRIPT_INTERFACE_WALBERLA_NONE_HPP

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameter.hpp"
#include "script_interface/walberla/EKPoissonSolver.hpp"

#include "LatticeWalberla.hpp"
#include "walberla_bridge/PoissonSolver/None.hpp"

#include <memory>

namespace ScriptInterface::walberla {

class EKNone : public EKPoissonSolver {
public:
  void do_construct(VariantMap const &args) override {
    m_noneinstance = std::make_shared<::walberla::None<double>>(
        get_value<std::shared_ptr<LatticeWalberla>>(args, "lattice")
            ->lattice());
  }

  [[nodiscard]] std::shared_ptr<::walberla::PoissonSolver<double>>
  get_instance() override {
    return m_noneinstance;
  }

private:
  /* The actual instance */
  std::shared_ptr<::walberla::None<double>> m_noneinstance;
};
} // namespace ScriptInterface::walberla

#endif // SCRIPT_INTERFACE_WALBERLA_NONE_HPP
