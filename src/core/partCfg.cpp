#include "partCfg.hpp"

#include <cassert>

PartCfg &partCfg(std::unique_ptr<PartCfg> init) {
  static std::unique_ptr<PartCfg> m_partCfg;

  if (init) {
    m_partCfg = std::move(init);
  }

  assert(m_partCfg);
  return *m_partCfg;
}
