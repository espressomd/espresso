#define BOOST_TEST_MODULE Histogram test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils.hpp"
#include "utils/Histogram.hpp"

BOOST_AUTO_TEST_CASE(histogram) {
  std::vector<size_t> n_bins{10, 10};
  std::vector<std::pair<double, double>> limits{
      {std::make_pair(1.0, 20.0), std::make_pair(5.0, 10.0)}};
  size_t n_dims_data = 2;
  auto hist = Utils::Histogram<double>(n_bins, n_dims_data, limits);
  // Check getters.
  BOOST_CHECK(hist.get_limits() == limits);
  BOOST_CHECK(hist.get_n_bins() == n_bins);
  BOOST_CHECK(
      (hist.get_bin_sizes() == std::vector<double>{{19.0 / 10.0, 5.0 / 10.0}}));
  // Check that histogram is initialized to zero.
  BOOST_CHECK(hist.get_histogram() ==
              std::vector<double>(n_dims_data * n_bins[0] * n_bins[1], 0.0));
  // Check that histogram still empty if data is out of bounds.
  hist.update(std::vector<double>{{1.0, 4.0}});
  BOOST_CHECK(hist.get_histogram() ==
              std::vector<double>(n_dims_data * n_bins[0] * n_bins[1], 0.0));
  // Check if putting in data at the first bin is set correctly.
  hist.update(std::vector<double>{{limits[0].first, limits[1].first}});
  BOOST_CHECK((hist.get_histogram())[0] == 1.0);
  BOOST_CHECK((hist.get_histogram())[1] == 1.0);
  // Check if weights are correctly set.
  hist.update(std::vector<double>{{limits[0].first, limits[1].first}},
              std::vector<double>{{10.0, 10.0}});
  BOOST_CHECK((hist.get_histogram())[0] == 11.0);
  BOOST_CHECK((hist.get_histogram())[1] == 11.0);
  BOOST_CHECK_THROW(hist.update(std::vector<double>{{1.0, 5.0, 3.0}}),
                    std::invalid_argument);
}
