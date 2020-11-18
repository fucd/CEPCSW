gaudi_subdir(CompareToMarlin v0r0)

find_package(EDM4HEP REQUIRED) 
find_package(AIDA REQUIRED)
find_package(CLHEP REQUIRED;CONFIG)

gaudi_depends_on_subdirs(
    Utilities/DataHelper
)

set(CompareToMarlin_srcs src/*.cpp)

# Modules
gaudi_add_module(CompareToMarlin ${CompareToMarlin_srcs}
    INCLUDE_DIRS k4FWCore AIDA
    LINK_LIBRARIES k4FWCore DataHelperLib GaudiKernel EDM4HEP ${CLHEP_LIBRARIES}
)
