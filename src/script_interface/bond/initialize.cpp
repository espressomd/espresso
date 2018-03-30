/*
  Copyright (C) 2015,2016 The ESPResSo project

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "initialize.hpp"

#include "Bonds.hpp"

//pairbonds
#include "Fene.hpp"
#include "Harmonic.hpp"
#include "HarmonicDumbbell.hpp"
#include "BondedCoulomb.hpp"
#include "BondedCoulombP3MSR.hpp"
#include "Quartic.hpp"
#include "SubtLj.hpp"
#include "Umbrella.hpp"
#include "TabulatedBondLength.hpp"
#include "OverlapBondLength.hpp"
#include "ThermalizedBond.hpp"

//3 particle bonds
#include "OverlapBondAngle.hpp"
#include "TabulatedBondAngle.hpp"
#include "AngleHarmonic.hpp"
#include "AngleCosine.hpp"
#include "AngleCosSquare.hpp"
#include "AngleDist.hpp"
#include "IbmTriel.hpp"

//4 particle bonds
#include "OverlapBondDihedral.hpp"
#include "TabulatedBondDihedral.hpp"
#include "HydrogenBond.hpp"
#include "Dihedral.hpp"
#include "MembraneCollision.hpp"
#include "OifLocalForces.hpp"

namespace ScriptInterface {
  namespace Bond {
    void initialize() {

      //pair bonds
      ScriptInterface::register_new<ScriptInterface::Bond::Bonds>("Bond::Bonds");
      ScriptInterface::register_new<ScriptInterface::Bond::Fene>("Bond::Fene");
      ScriptInterface::register_new<ScriptInterface::Bond::Harmonic>("Bond::Harmonic");
      ScriptInterface::register_new<ScriptInterface::Bond::HarmonicDumbbell>
	("Bond::HarmonicDumbbell");
      ScriptInterface::register_new<ScriptInterface::Bond::BondedCoulomb>
	("Bond::BondedCoulomb");
      ScriptInterface::register_new<ScriptInterface::Bond::BondedCoulombP3MSR>
	("Bond::BondedCoulombP3MSR");
      ScriptInterface::register_new<ScriptInterface::Bond::Quartic>("Bond::Quartic");
      ScriptInterface::register_new<ScriptInterface::Bond::SubtLj>("Bond::SubtLj");
      ScriptInterface::register_new<ScriptInterface::Bond::Umbrella>("Bond::Umbrella");
      ScriptInterface::register_new<ScriptInterface::Bond::TabulatedBondLength>
	("Bond::TabulatedBondLength");
      ScriptInterface::register_new<ScriptInterface::Bond::OverlapBondLength>
	("Bond::OverlapBondLength");
      ScriptInterface::register_new<ScriptInterface::Bond::ThermalizedBond>
	("Bond::ThermalizedBond");

      //3 partricle bonds
      ScriptInterface::register_new<ScriptInterface::Bond::AngleHarmonic>
	("Bond::AngleHarmonic");
      ScriptInterface::register_new<ScriptInterface::Bond::AngleCosine>
	("Bond::AngleCosine");
      ScriptInterface::register_new<ScriptInterface::Bond::AngleCosSquare>
	("Bond::AngleCosSquare");
      ScriptInterface::register_new<ScriptInterface::Bond::OverlapBondAngle>
	("Bond::OverlapBondAngle");
      ScriptInterface::register_new<ScriptInterface::Bond::TabulatedBondAngle>
	("Bond::TabulatedBondAngle");
      ScriptInterface::register_new<ScriptInterface::Bond::AngleDist>
	("Bond::AngleDist");
      ScriptInterface::register_new<ScriptInterface::Bond::IbmTriel>
	("Bond::IbmTriel");

      //4 partricle bonds
      ScriptInterface::register_new<ScriptInterface::Bond::HydrogenBond>
	("Bond::HydrogenBond");
      ScriptInterface::register_new<ScriptInterface::Bond::OverlapBondDihedral>
	("Bond::OverlapBondDihedral");
      ScriptInterface::register_new<ScriptInterface::Bond::TabulatedBondDihedral>
	("Bond::TabulatedBondDihedral");
      ScriptInterface::register_new<ScriptInterface::Bond::Dihedral>
	("Bond::Dihedral");
      ScriptInterface::register_new<ScriptInterface::Bond::MembraneCollision>
	("Bond::MembraneCollision");
      ScriptInterface::register_new<ScriptInterface::Bond::OifLocalForces>
	("Bond::OifLocalForces");
    }
    
  } /* namespace Shapes */
} /* namespace ScriptInterface */
