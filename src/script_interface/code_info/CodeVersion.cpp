/*
 * Copyright (C) 2022-2026 The ESPResSo project
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

#include "CodeVersion.hpp"

#include "config/version.hpp"

#include <string>
#include <vector>

namespace ScriptInterface {
namespace CodeInfo {

Variant CodeVersion::do_call_method(std::string const &name,
                                    VariantMap const &) {
  if (name == "version_major") {
    return ESPRESSO_VERSION_MAJOR;
  }
  if (name == "version_minor") {
    return ESPRESSO_VERSION_MINOR;
  }
  if (name == "version_friendly") {
    return std::string(ESPRESSO_VERSION);
  }
  if (name == "version") {
    return std::vector<int>{ESPRESSO_VERSION_MAJOR, ESPRESSO_VERSION_MINOR,
                            ESPRESSO_VERSION_PATCH};
  }
  if (name == "git_branch") {
    return std::string{GIT_BRANCH};
  }
  if (name == "git_commit") {
    return std::string{GIT_COMMIT_HASH};
  }
  if (name == "git_state") {
    return std::string{GIT_STATE};
  }
  return {};
}

} // namespace CodeInfo
} // namespace ScriptInterface
