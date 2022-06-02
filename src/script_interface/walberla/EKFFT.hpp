#ifndef SCRIPT_INTERFACE_WALBERLA_FFT_HPP
#define SCRIPT_INTERFACE_WALBERLA_FFT_HPP

#include "EKPoissonSolver.hpp"
#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameter.hpp"

#include "walberla_bridge/electrokinetics/ek_walberla_init.hpp"

#include <memory>

namespace ScriptInterface::walberla {

class EKFFT : public EKPoissonSolver {
public:
  void do_construct(VariantMap const &args) override {
    m_single_precision = get_value_or<bool>(args, "single_precision", false);
    auto lattice =
        get_value<std::shared_ptr<LatticeWalberla>>(args, "lattice")->lattice();
    auto permittivity = get_value<double>(args, "permittivity");

    m_fftinstance =
        new_ek_poisson_fft(lattice, permittivity, m_single_precision);

    add_parameters({{"permittivity",
                     [this](Variant const &v) {
                       m_fftinstance->set_permittivity(get_value<double>(v));
                     },
                     [this]() { return m_fftinstance->get_permittivity(); }},
                    {"single_precision", AutoParameter::read_only,
                     [this]() { return m_single_precision; }}});
  }

  [[nodiscard]] std::shared_ptr<::walberla::PoissonSolver> get_instance() const
      noexcept override {
    return m_fftinstance;
  }

private:
  /* The actual instance */
  std::shared_ptr<::walberla::PoissonSolver> m_fftinstance;

  bool m_single_precision;
};
} // namespace ScriptInterface::walberla

#endif // SCRIPT_INTERFACE_WALBERLA_FFT_HPP
