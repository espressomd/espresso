/*
 * Copyright (C) 2015-2022 The ESPResSo project
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

void initialize(Utils::Factory<ObjectHandle> *om) {
  om->register_new<Constraints>("Constraints::Constraints");
  om->register_new<ShapeBasedConstraint>("Constraints::ShapeBasedConstraint");
  om->register_new<HomogeneousMagneticField>(
      "Constraints::HomogeneousMagneticField");
  om->register_new<TabulatedForceField>("Constraints::ForceField");
  om->register_new<TabulatedPotentialField>("Constraints::PotentialField");
  om->register_new<Gravity>("Constraints::Gravity");
  om->register_new<FlowField>("Constraints::FlowField");
  om->register_new<HomogeneousFlowField>("Constraints::HomogeneousFlowField");
#ifdef ELECTROSTATICS
  om->register_new<ElectricPotential>("Constraints::ElectricPotential");
  om->register_new<LinearElectricPotential>(
      "Constraints::LinearElectricPotential");
  om->register_new<ElectricPlaneWave>("Constraints::ElectricPlaneWave");
#endif
}
} /* namespace Constraints */
} /* namespace ScriptInterface */
