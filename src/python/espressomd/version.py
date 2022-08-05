#
# Copyright (C) 2010-2022 The ESPResSo project
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
#

from .script_interface import script_interface_register, ScriptInterfaceHelper


@script_interface_register
class _Version(ScriptInterfaceHelper):
    _so_name = "CodeInfo::CodeVersion"
    _so_creation_policy = "LOCAL"
    _so_bind_methods = (
        "version",
        "version_major",
        "version_minor",
        "version_friendly",
        "git_branch",
        "git_commit",
        "git_state",
    )


def version():
    """Get the version of ESPResSo."""
    return tuple(_Version().version())


def major():
    """Get the major version of ESPResSo."""
    return _Version().version_major()


def minor():
    """Get the minor version of ESPResSo."""
    return _Version().version_minor()


def friendly():
    """Dot version of the version."""
    return _Version().version_friendly()


def git_branch():
    """Git branch of the build if known, otherwise empty."""
    return _Version().git_branch()


def git_commit():
    """Git commit of the build if known, otherwise empty."""
    return _Version().git_commit()


def git_state():
    """
    Git state of the build if known, otherwise empty. State is "CLEAN" if the
    repository was not changed from :meth:`git_commit()`, "DIRTY" otherwise.
    """
    return _Version().git_state()
