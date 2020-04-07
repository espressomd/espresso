#include "bond_error.hpp"

#include "errorhandling.hpp"

void bond_broken_error(int id, Utils::Span<const int> partner_ids) {
  auto error_msg = runtimeErrorMsg();

  error_msg << "bond broken between particles " << id;
  for (auto partner_id : partner_ids) {
    error_msg << ", " << partner_id;
  }
}
