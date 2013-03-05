#!/bin/sh
#
# Copyright (C) 2012,2013 The ESPResSo project
# Copyright (C) 2012,2013 Olaf Lenz
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
# Generate the animated logo.
#

if ! test -d cup-animation; then
    echo "Directory cup-animation doesn't exist!"
    exit 1
fi

echo "GENERATING ANIMATED GIF..."

echo "  Combining logo box with cup in all frames..."
for i in `seq 24 49`; do
    convert -page 500x500 -gravity NorthWest logo-template.png cup-animation/$i.png -geometry +113+0 -composite cup-animation/logo-$i.png
done

echo "  Generating animation..."
for i in `seq 25 49; seq 48 -1 25`; do
    sequence="$sequence cup-animation/logo-$i.png"
done
convert \
    -dispose None \
    -delay 0 \
    cup-animation/logo-24.png \
    -dispose previous \
    -delay 5 \
    $sequence \
    -loop 0 \
    -coalesce \
    cup-animation/logo-animated-500.miff

echo "  Generating different sizes..."
convert -layers Optimize cup-animation/logo-animated-500.miff logo-animated-500.gif
convert -size 500x500 -layers Optimize cup-animation/logo-animated-500.miff -resize 200x200 logo-animated-200.gif
convert -size 500x500 -layers Optimize cup-animation/logo-animated-500.miff -resize 100x100 logo-animated-100.gif
ls -l logo-animated-500.gif logo-animated-200.gif logo-animated-100.gif

echo "Finished."
