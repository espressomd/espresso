
cdef class PShape(PScriptInterface):
    def __init__(self, name):
        self.si = factory_make[Shape](name)
        self.sip = <ScriptInterfaceBase *> self.si.get()
        super(PShape, self).__init__()


cdef class PWall(PShape):
    def __init__(self):
        super(PWall, self).__init__("wall")

