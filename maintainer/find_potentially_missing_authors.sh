# Copyright (C) 2019 The ESPResSo project
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

# This gets a tentative list of authors that are missing from the AUTHORS file
# based on git commits. The output has to be checked manually, because
# users have committed under different spellings and abbreviations.

# This has to be run from the root directory of the source tree

git log|
grep -i ^author|
cut -f2- -d\ |
sed -e 's/ <.*//'|
sort -u|
while read author
do 
  grep -i "$author" AUTHORS >/dev/null || echo Missing: $author
done

