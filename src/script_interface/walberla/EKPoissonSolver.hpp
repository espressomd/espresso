#ifndef SCRIPT_INTERFACE_WALBERLA_EKPOISSONSOLVER_HPP
#define SCRIPT_INTERFACE_WALBERLA_EKPOISSONSOLVER_HPP

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameter.hpp"

#include "walberla_bridge/electrokinetics/PoissonSolver/PoissonSolver.hpp"

#include <memory>

namespace ScriptInterface {
namespace walberla {

class EKPoissonSolver
    : public AutoParameters<::walberla::PoissonSolver<double>> {
public:
  virtual std::shared_ptr<::walberla::PoissonSolver<double>> get_instance() = 0;
};
} // namespace walberla
} // namespace ScriptInterface

#endif // SCRIPT_INTERFACE_WALBERLA_EKPOISSONSOLVER_HPP
