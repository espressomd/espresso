from functools import update_wrapper

class ThereCanOnlyBeOne(BaseException):
    def __init__(self, cls):
        self._cls = cls
    def __str__(self):
        return "There can only be one instance of '{}' at any time.".format(self._cls)

def highlander(klass):
    klass.highlander_created = False

    def cls_init(self, *args, **kwargs):
        "__init__ method by the highlander decorator."
        if self.__class__.highlander_created:
            raise ThereCanOnlyBeOne(self.__class__)
        self.__class__.highlander_created = True
    def cls_init_call_orig(self, *args, **kwargs):
        if self.__class__.highlander_created:
            raise ThereCanOnlyBeOne(self.__class__)
        self.__class__.highlander_created = True
        self.__class__.__init_orig__(self, *args, **kwargs)

    # override the __init__ method of the class to store the bool
    # "highlander_created"
    if hasattr(klass, '__init__'):
        klass.__init_orig__ = klass.__init__
        klass.__init__ = cls_init_call_orig
        update_wrapper(cls_init_call_orig, klass.__init_orig__)
    else:
        klass.__init__ = cls_init
        
    # override the __del__ method of the class
    def cls_del(self):
        "__del__ method by the highlander decorator."
        self.__class__.highlander_created = False
    def cls_del_call_orig(self):
        cls_del(self)
        self.__class__.__del_orig__(self)

    if hasattr(klass, '__del__'):
        klass.__del_orig__ = klass.__del__
        klass.__del__ = cls_del_call_orig
        update_wrapper(cls_del_call_orig, klass.__del_orig__)
    else:
        klass.__del__ = cls_del

    return klass

@highlander
class A:
    def __init__(self, x):
        "Create it."
        self.x = x
    def __del__(self):
        print("A.del called")

@highlander
class B:
    pass

