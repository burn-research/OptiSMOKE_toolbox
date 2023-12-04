# ------------------------------------------------------------------------- #
#  Found somewhere don't remember where and slighlty modified               #
# ------------------------------------------------------------------------- #

if(NOT Config_INCLUDE_DIR)
    find_path(Config_INCLUDE_DIR
        NAMES
        libconfig.h++
        HINTS
        ENV Config_ROOT/include)
endif()

if(NOT Config_LIBRARIES)
    find_library(Config_LIBRARIES
        NAMES 
        libconfig++.a
        HINTS
        ENV Config_ROOT/lib)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Config DEFAULT_MSG Config_INCLUDE_DIR)
find_package_handle_standard_args(Config DEFAULT_MSG Config_LIBRARIES)

mark_as_advanced(Config_INCLUDE_DIR)
mark_as_advanced(Config_LIBRARIES)
