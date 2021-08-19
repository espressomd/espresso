#ifndef SCRIPT_INTERFACE_WALBERLA_FFT_HPP
#define SCRIPT_INTERFACE_WALBERLA_FFT_HPP

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameter.hpp"

#include "grid_based_algorithms/walberla_blockforest.hpp"
#include "walberla_bridge/PoissonSolver/FFT.hpp"

#include <memory>

namespace ScriptInterface::walberla {

class EKFFT : public AutoParameters<::walberla::FFT<double>> {
public:
  void do_construct(VariantMap const &args) override {
    m_fftinstance =
        std::make_shared<::walberla::FFT<double>>(get_walberla_blockforest());
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
