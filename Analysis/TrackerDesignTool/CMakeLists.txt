include(DD4hep)

message("fucd ${ROOT_LIBRARIES}")
# Modules
gaudi_add_module(TrackerDesignTool
                 SOURCES src/TrackerDesignToolAlg.cpp src/TDTHit.cpp
                 LINK DataHelperLib 
                      Gaudi::GaudiKernel 
                      EDM4HEP::edm4hep 
                      ${ROOT_LIBRARIES}
                      ${CLHEP_LIBRARIES}
		      DetInterface
)

install(TARGETS TrackerDesignTool
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)
