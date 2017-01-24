try:
    from .visualizationMayavi import mayaviLive
except ImportError as e:
    class mayaviLive:
        def __init__(*args, **kwargs):
            raise e

try:
    from .visualizationOpenGL import openGLLive
except ImportError as e:
    class openGLLive:
        def __init__(*args, **kwargs):
            raise e

__all__ = ['mayaviLive', 'openGLLive']
