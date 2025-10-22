#
# Copyright (C) 2019-2024 The ESPResSo project
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

file(READ "${INPUT}" FILE_CONTENTS)

# transform BibTeX DOI fields into URL fields (bibliographic styles available
# to Doxygen do not process the DOI field)
string(REGEX REPLACE
       "([\r\n])[\t ]*doi *= *([\\{\"]+)(10\\.[0-9]+)"
       "\\1url=\\2https://doi.org/\\3"
       FILE_CONTENTS "${FILE_CONTENTS}")

file(WRITE "${OUTPUT}" "${FILE_CONTENTS}")
