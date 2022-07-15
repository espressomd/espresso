#!/bin/sh
#
# Copyright (C) 2021-2022 The ESPResSo project
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

pypresso=../../build/pypresso

if [ ! -e "${pypresso}" ]
then
    echo "Invalid path to pypresso script: ${pypresso}" 1>&2
    exit 1
fi

# Run the simulations for the temperatures given as parameters
for kT in $*
do
    for seed in 1 2 3 4 5
    do

        test -s "temp_${kT}_seed_${seed}.dat.gz" ||
            "${pypresso}" run_sim.py ${kT} --log --steps 1000000 --seed ${seed} > "temp_${kT}_seed_${seed}.out" 2>&1 &
    done
done
wait

# Plot the results, create fits
"${pypresso}" create_fits.py temp_*.dat.gz
