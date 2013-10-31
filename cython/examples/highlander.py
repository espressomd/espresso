import gc

class ThereCanOnlyBeOne():
    def __init__(self, cls):
        self._cls = cls
    def __str__(self):
        return "Can create only one instance of '%s'." % self._cls
        
class Highlander(object):
    _instance_count = 0
    def __init__(self):
        # if there is more than one instance, it might not have been
        # garbage collected yet, so start collection now
        if self.__class__._instance_count > 0: gc.collect()
        # afterwards, check once more
        if self.__class__._instance_count > 0: 
            raise ThereCanOnlyBeOne(self.__class__)
        #
        self.__class__._instance_count += 1

    def __del__(self):
        self.__class__._instance_count -= 1

class A(Highlander):
    def __init__(self):
        pass
    def __del__(self):
        Highlander.__del__(self)

h1 = Highlander()
print A._instance_count
del h1
h2 = Highlander()
try:
    h3 = Highlander()
except:
    pass

print A._instance_count
a1 = A()
print A._instance_count
del a1
print A._instance_count
a2 = A()
print A._instance_count

