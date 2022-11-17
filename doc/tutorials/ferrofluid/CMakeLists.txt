#
# Copyright (C) 2018-2022 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

configure_tutorial_target(
  TARGET tutorial_fe DEPENDS ferrofluid_part1.ipynb ferrofluid_part2.ipynb
  ferrofluid_part3.ipynb figures/Electro-Steric_Stabilization.jpg
  figures/Ferrofluid_Magnet_under_glass_edit.jpg figures/headtotailconf.png)

nb_export(
  TARGET
  tutorial_fe
  SUFFIX
  "1"
  FILE
  "ferrofluid_part1.ipynb"
  HTML_RUN
  VAR_SUBST
  "EQUIL_STEPS=100;EQUIL_ROUNDS=10;\"CI_DP3M_PARAMS={'cao':3,'r_cut':8.34,'mesh':[8,8,8],'alpha':0.2115,'tune':False}\""
)
nb_export(TARGET tutorial_fe SUFFIX "2" FILE "ferrofluid_part2.ipynb" HTML_RUN
          VAR_SUBST
          "equil_steps=100;equil_rounds=10;alphas=[0,1,2,3,4.5,8];loops=100")
nb_export(
  TARGET tutorial_fe SUFFIX "3" FILE "ferrofluid_part3.ipynb" HTML_RUN
  VAR_SUBST
  "equil_steps=100;equil_rounds=10;alphas=[0,0.5,1,2,4,8];loops=200;loops_m=100"
)
