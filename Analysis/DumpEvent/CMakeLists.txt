gaudi_subdir(DumpEvent v0r0)

find_package(EDM4HEP REQUIRED) 
find_package(AIDA REQUIRED)
find_package(CLHEP REQUIRED;CONFIG)

gaudi_depends_on_subdirs(
    Utilities/DataHelper
)

set(DumpEvent_srcs src/*.cpp)

# Modules
gaudi_add_module(DumpEvent ${DumpEvent_srcs}
    INCLUDE_DIRS AIDA
    LINK_LIBRARIES DataHelperLib GaudiKernel EDM4HEP ${CLHEP_LIBRARIES}
)
