/*
 * Copyright (C) 2010-2020 The ESPResSo project
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
#ifndef CORE_IO_WRITER_H5MD_SPECIFICATION_HPP
#define CORE_IO_WRITER_H5MD_SPECIFICATION_HPP

#include <h5xx/h5xx.hpp>

#include <algorithm>
#include <string>
#include <vector>

namespace Writer {
namespace H5md {

/**
 * @brief Layout information for H5MD files.
 * In order to add a new particle property you have to add an entry to the
 * H5MD_Specification::DATASETS member and extend the File::write() and the
 * File::write_units() functions accordingly.
 */
struct H5MD_Specification {

  struct Dataset {
    std::string path() const { return group + "/" + name; }

    std::string group;
    std::string name;
    hsize_t rank;
    hid_t type;
    hsize_t data_dim;
    bool is_link;
  };

  H5MD_Specification(unsigned int fields);

  auto const &get_datasets() const { return m_datasets; }

  bool is_compliant(std::string const &filename) const {
    h5xx::file h5md_file(filename, h5xx::file::in);

    auto const all_groups_exist = std::all_of(
        m_datasets.begin(), m_datasets.end(), [&h5md_file](auto const &d) {
          return h5xx::exists_group(h5md_file, d.group);
        });
    auto const all_datasets_exist = std::all_of(
        m_datasets.begin(), m_datasets.end(), [&h5md_file](auto const &d) {
          return h5xx::exists_dataset(h5md_file, d.path());
        });
    return all_groups_exist and all_datasets_exist;
  }

private:
  std::vector<Dataset> m_datasets;
};

} // namespace H5md
} // namespace Writer

#endif
