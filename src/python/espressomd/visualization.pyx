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
import sys
import os

try:
    if sys.platform == "darwin" or os.name != "posix" or "DISPLAY" in os.environ:
        from .visualization_mayavi import mayaviLive
    else:
        raise ImportError("Cannot connect to X server")
except BaseException as e:
    if isinstance(e, ImportError) or isinstance(e, RuntimeError) and \
            e.args[0] == 'No pyface.toolkits plugin could be loaded for wx':
        class mayaviLive:
            deferred_ImportError = e

            def __init__(self, *args, **kwargs):
                raise self.deferred_ImportError
    else:
        raise e

try:
    from .visualization_opengl import openGLLive
except ImportError as e:
    class openGLLive:
        deferred_ImportError = e

        def __init__(self, *args, **kwargs):
            raise self.deferred_ImportError

__all__ = ['mayaviLive', 'openGLLive']
