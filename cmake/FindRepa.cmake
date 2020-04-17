
set(REPA_DIR "" CACHE PATH "Directory of librepa")

find_path(REPA_INCLUDE_DIR
          repa/repa.hpp
          HINTS ${REPA_DIR}
          ENV C_INCLUDE_PATH
          PATH_SUFFIXES include)

find_library(REPA_LIBRARIES
             repa
             HINTS ${REPA_DIR}
             ENV LIBRARY_PATH
             PATH_SUFFIXES lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(REPA
                                  DEFAULT_MSG
                                  REPA_LIBRARIES
                                  REPA_INCLUDE_DIR)
mark_as_advanced(REPA_LIBRARIES REPA_INCLUDE_DIR)