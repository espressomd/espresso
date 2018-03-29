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

#include "Fene.hpp"
#include "Harmonic.hpp"
#include "HarmonicDumbbell.hpp"
#include "BondedCoulomb.hpp"
#include "BondedCoulombP3MSR.hpp"
#include "Quartic.hpp"
#include "SubtLj.hpp"
#include "Umbrella.hpp"
#include "TabulatedBondLength.hpp"
//#include "OverlapBondLength.hpp"
#include "ThermalizedBond.hpp"

namespace ScriptInterface {
  namespace Bond {
    void initialize() {
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
      //      ScriptInterface::register_new<ScriptInterface::Bond::OverlapBondLength>
      //("Bond::OverlapBondLength");
      ScriptInterface::register_new<ScriptInterface::Bond::ThermalizedBond>
	("Bond::ThermalizedBond");
    }
    
  } /* namespace Shapes */
} /* namespace ScriptInterface */
