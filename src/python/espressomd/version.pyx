# Copyright (C) 2010-2019 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from .utils import to_str


def major():
    """Prints the major version of ESPResSo.
    """
    return ESPRESSO_VERSION_MAJOR  # pylint: disable=undefined-variable


def minor():
    """Prints the minor version of ESPResSo.
    """
    return ESPRESSO_VERSION_MINOR  # pylint: disable=undefined-variable


def friendly():
    """Dot version of the version.
    """
    return "{}.{}".format(major(), minor())


def git_branch():
    """Git branch of the build if known, otherwise
       empty.
    """
    return to_str(GIT_BRANCH)  # pylint: disable=undefined-variable


def git_commit():
    """Git commit of the build if known, otherwise
       empty.
    """
    return to_str(GIT_COMMIT_HASH)  # pylint: disable=undefined-variable


def git_state():
    """Git state of the build if known, otherwise
       empty. State is "CLEAN" if the repository
       was not changed from :meth:`git_commit()`,
       "DIRTY" otherwise.
    """
    return to_str(GIT_STATE)  # pylint: disable=undefined-variable
