#include "PartCfg.hpp"

#include "grid.hpp"
#include "particle_data.hpp"

#include <boost/algorithm/clamp.hpp>

void PartCfg::update() {
  if (m_valid)
    return;

  remote_parts.clear();

  auto const ids = get_particle_ids();
  auto const chunk_size = fetch_cache_max_size();

  for (size_t offset = 0; offset < ids.size();) {
    auto const this_size =
        boost::algorithm::clamp(chunk_size, 0, ids.size() - offset);
    auto const chunk_ids =
        Utils::make_const_span(ids.data() + offset, this_size);

    prefetch_particle_data(chunk_ids);

    for (auto id : chunk_ids) {
      remote_parts.push_back(get_particle_data(id));

      auto &p = remote_parts.back();
      p.r.p += image_shift(p.l.i, box_geo.length());
      p.l.i = {};
    }

    offset += this_size;
  }

  m_valid = true;
}
