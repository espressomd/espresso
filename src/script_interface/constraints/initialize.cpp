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

#include "script_interface/ClassName.hpp"

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

constexpr auto class_names() {
  return std::make_tuple(
      ClassName<Constraints>{"Constraints::Constraints"},
      ClassName<ShapeBasedConstraint>{"Constraints::ShapeBasedConstraint"},
      ClassName<HomogeneousMagneticField>{
          "Constraints::HomogeneousMagneticField"},
      ClassName<TabulatedForceField>{"Constraints::ForceField"},
      ClassName<TabulatedPotentialField>{"Constraints::PotentialField"},
      ClassName<Gravity>{"Constraints::Gravity"},
      ClassName<FlowField>{"Constraints::FlowField"},
      ClassName<HomogeneousFlowField> { "Constraints::HomogeneousFlowField" }
#ifdef ELECTROSTATICS
      ,
      ClassName<ElectricPotential>{"Constraints::ElectricPotential"},
      ClassName<LinearElectricPotential>{
          "Constraints::LinearElectricPotential"},
      ClassName<ElectricPlaneWave> { "Constraints::ElectricPlaneWave" }
#endif
  );
}

void initialize(Utils::Factory<ObjectHandle> *om) {
  Utils::for_each(
      [om](auto name) {
        om->register_new<typename decltype(name)::class_type>(name.name);
      },
      class_names());
}
} /* namespace Constraints */
} /* namespace ScriptInterface */
