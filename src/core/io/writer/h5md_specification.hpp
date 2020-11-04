#ifndef CORE_IO_WRITER_H5MD_HPP
#define CORE_IO_WRITER_H5MD_HPP

#include <array>

#include <h5xx/h5xx.hpp>

#include <algorithm>
#include <string>

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

  static std::array<Dataset, 30> DATASETS;

  static bool is_compliant(std::string const &filename) {
    h5xx::file h5md_file(filename, h5xx::file::in);

    auto all_groups_exist = std::all_of(
        DATASETS.begin(), DATASETS.end(), [&h5md_file](auto const &dataset) {
          return h5xx::exists_group(h5md_file, dataset.group);
        });
    auto all_datasets_exist = std::all_of(
        DATASETS.begin(), DATASETS.end(), [&h5md_file](auto const &dataset) {
          return h5xx::exists_dataset(h5md_file, dataset.path());
        });
    return all_groups_exist and all_datasets_exist;
  }
};

} // namespace H5md
} // namespace Writer

#endif
