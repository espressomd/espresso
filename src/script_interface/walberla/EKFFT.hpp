#ifndef SCRIPT_INTERFACE_WALBERLA_FFT_HPP
#define SCRIPT_INTERFACE_WALBERLA_FFT_HPP

#include "EKPoissonSolver.hpp"
#include "LatticeWalberla.hpp"
#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameter.hpp"

#include "grid_based_algorithms/walberla_blockforest.hpp"
#include "walberla_bridge/PoissonSolver/FFT.hpp"

#include <memory>

namespace ScriptInterface::walberla {

class EKFFT : public EKPoissonSolver {
public:
  void do_construct(VariantMap const &args) override {
    m_fftinstance = std::make_shared<::walberla::FFT<double>>(
        get_value<std::shared_ptr<LatticeWalberla>>(args, "lattice")->lattice(),
        get_value<double>(args, "permittivity"));

    add_parameters({{"permittivity",
                     [this](Variant const &v) {
                       m_fftinstance->set_permittivity(get_value<double>(v));
                     },
                     [this]() { return m_fftinstance->get_permittivity(); }}});
  }

  [[nodiscard]] std::shared_ptr<::walberla::PoissonSolver<double>>
  get_instance() override {
    return m_fftinstance;
  }

  [[nodiscard]] std::shared_ptr<::walberla::FFT<double>> get_fftinstance() {
    return m_fftinstance;
  }

private:
  /* The actual instance */
  std::shared_ptr<::walberla::FFT<double>> m_fftinstance;
};
} // namespace ScriptInterface::walberla

#endif // SCRIPT_INTERFACE_WALBERLA_FFT_HPP
