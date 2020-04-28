/*
 * Copyright (C) 2015-2019 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "initialize.hpp"
#include "script_interface/ScriptInterface.hpp"

#include "Constraints.hpp"

#include "HomogeneousMagneticField.hpp"
#include "ShapeBasedConstraint.hpp"

#include "ExternalField.hpp"
#include "ExternalPotential.hpp"

#include "couplings.hpp"
#include "fields.hpp"

namespace ScriptInterface {
namespace Constraints {

using namespace FieldCoupling::Coupling;
using namespace FieldCoupling::Fields;

/* Generic Fields */
using TabulatedForceField = ExternalField<Scaled, Interpolated<double, 3>>;
using TabulatedPotentialField =
    ExternalPotential<Scaled, Interpolated<double, 1>>;

/* Physical Fields */
using Gravity = ExternalField<Mass, Constant<double, 3>>;

using FlowField = ExternalField<Viscous, Interpolated<double, 3>>;
using LinearFlowField = ExternalField<Viscous, AffineMap<double, 3>>;
using HomogeneousFlowField = ExternalField<Viscous, Constant<double, 3>>;

using ElectricPotential = ExternalPotential<Charge, Interpolated<double, 1>>;
using LinearElectricPotential = ExternalPotential<Charge, AffineMap<double, 1>>;
using ElectricPlaneWave = ExternalField<Charge, PlaneWave<double, 3>>;

void initialize() {
  ScriptInterface::register_new<ScriptInterface::Constraints::Constraints>(
      "Constraints::Constraints");

  ScriptInterface::register_new<
      ScriptInterface::Constraints::ShapeBasedConstraint>(
      "Constraints::ShapeBasedConstraint");

  ScriptInterface::register_new<
      ScriptInterface::Constraints::HomogeneousMagneticField>(
      "Constraints::HomogeneousMagneticField");

  ScriptInterface::register_new<TabulatedForceField>("Constraints::ForceField");
  ScriptInterface::register_new<TabulatedPotentialField>(
      "Constraints::PotentialField");

  ScriptInterface::register_new<Gravity>("Constraints::Gravity");
  ScriptInterface::register_new<FlowField>("Constraints::FlowField");
  ScriptInterface::register_new<HomogeneousFlowField>(
      "Constraints::HomogeneousFlowField");

#ifdef ELECTROSTATICS
  ScriptInterface::register_new<ElectricPotential>(
      "Constraints::ElectricPotential");
  ScriptInterface::register_new<LinearElectricPotential>(
      "Constraints::LinearElectricPotential");
  ScriptInterface::register_new<ElectricPlaneWave>(
      "Constraints::ElectricPlaneWave");
#endif
}
} /* namespace Constraints */
} /* namespace ScriptInterface */
