#include "ek_boundaries.hpp"
#include "config.hpp"
#include "ek_container.hpp"
#include "ek_interface.hpp"
#include "event.hpp"
#include "walberla_blockforest.hpp"

#include <vector>

namespace EKBoundaries {

std::vector<std::shared_ptr<EKBoundary>> ekboundaries;
void ek_init_boundaries() {
#ifdef EK_WALBERLA
  std::for_each(
      EK::ek_container.begin(), EK::ek_container.end(), [](const auto ek) {
        auto const agrid = get_walberla_blockforest_params()->get_agrid();
        ek->clear_boundaries();

        for (auto index_and_pos : ek->node_indices_positions(true)) {
          // Convert to MD units
          auto const index = index_and_pos.first;
          auto const pos = index_and_pos.second * agrid;

          for (auto const &ekboundary : ekboundaries) {
            if (not ekboundary->shape().is_inside(pos)) {
              ek->set_node_noflux_boundary(index);
            }
          }
        }
      });
#endif // EK_WALBERLA
}

void add(const std::shared_ptr<EKBoundary> &b) {
  auto &ekb = ekboundaries;
  assert(std::find(ekb.begin(), ekb.end(), b) == ekb.end());
  ekb.emplace_back(b);

  on_ekboundary_change();
}

void remove(const std::shared_ptr<EKBoundary> &b) {
  auto &ekb = ekboundaries;
  assert(std::find(ekb.begin(), ekb.end(), b) != ekb.end());
  ekb.erase(std::remove(ekb.begin(), ekb.end(), b), ekb.end());

  on_ekboundary_change();
}
} // namespace EKBoundaries