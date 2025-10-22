/*
 * Copyright (C) 2010-2025 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#pragma once

#include <filesystem>
#include <vector>

namespace Writer {
namespace H5md {

struct Dataset;

/**
 * @brief Layout information for H5MD files.
 * In order to add a new particle property you have to add an entry to the
 * H5MD_Specification::DATASETS member and extend the File::write() and the
 * File::write_units() functions accordingly.
 */
struct Specification {
  Specification(unsigned int fields);

  auto const &get_datasets() const { return m_datasets; }

  bool is_compliant(std::filesystem::path const &file) const;

private:
  std::vector<Dataset> m_datasets;
};

} // namespace H5md
} // namespace Writer
