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
        class mayaviLive(object):
            deferred_ImportError = e

            def __init__(self, *args, **kwargs):
                raise self.deferred_ImportError
    else:
        raise e

try:
    from .visualization_opengl import openGLLive
except ImportError as e:
    class openGLLive(object):
        deferred_ImportError = e

        def __init__(self, *args, **kwargs):
            raise self.deferred_ImportError

__all__ = ['mayaviLive', 'openGLLive']
