#
# Copyright (C) 2023-2024 The ESPResSo project
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

if(EXISTS "${BIBLIOGRAPHY}")
  file(READ "${BIBLIOGRAPHY}" FILE_CONTENTS)
  # convert unsupported LaTeX formatting commands to HTML
  string(REGEX REPLACE " textsuperscript ([^ ]+) " "<sup>\\1</sup> " FILE_CONTENTS "${FILE_CONTENTS}")
  string(REGEX REPLACE " textsubscript ([^ ]+) " "<sub>\\1</sub> " FILE_CONTENTS "${FILE_CONTENTS}")
  file(WRITE "${BIBLIOGRAPHY}" "${FILE_CONTENTS}")
endif()
