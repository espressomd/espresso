import sys, os

try:
    if sys.platform == "darwin" or os.name != "posix" or "DISPLAY" in os.environ:
        from .visualization_mayavi import mayaviLive
    else:
        raise ImportError("Cannot connect to X server")
except ImportError as e:
    class mayaviLive(object):
        def __init__(*args, **kwargs):
            raise e

try:
    from .visualization_opengl import openGLLive
except ImportError as e:
    class openGLLive(object):
        def __init__(*args, **kwargs):
            raise e

__all__ = ['mayaviLive', 'openGLLive']
