#ifndef ESPRESSO_BOND_ERROR_HPP
#define ESPRESSO_BOND_ERROR_HPP

#include <utils/Span.hpp>

#include <stdexcept>

void bond_broken_error(int id, Utils::Span<const int> partner_ids);

/**
 * Exception indicating that a particle id
 * could not be resolved.
 */
struct BondResolutionError : std::exception {};

/**
 * Exception indicating that a bond type
 * was unknown.
 */
struct BondUnknownTypeError : std::exception {
  explicit BondUnknownTypeError(int type) : type(type) {}

  int type;
};

/**
 * Exception indicating that a bond with an
 * unexpected number of partners was encountered.
 */
struct BondInvalidSizeError : std::exception {
  explicit BondInvalidSizeError(int size) : size(size) {}

  int size;
};

#endif // ESPRESSO_BOND_ERROR_HPP
