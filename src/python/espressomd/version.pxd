cdef extern from "version.hpp":
    cdef int ESPRESSO_VERSION_MAJOR
    cdef int ESPRESSO_VERSION_MINOR
    cdef const char * GIT_BRANCH
    cdef const char * GIT_COMMIT_HASH
    cdef const char * GIT_STATE
