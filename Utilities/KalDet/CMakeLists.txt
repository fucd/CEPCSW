##############################################################################
# Package: KalDet
#    Desc: import from ILCSoft
##############################################################################

gaudi_subdir(KalDet v0r0)

find_package(CLHEP REQUIRED;CONFIG)
find_package(LCIO)
find_package(GEAR)
find_package(ROOT COMPONENTS MathCore)
find_package(EDM4HEP REQUIRED)
find_package(DD4hep COMPONENTS DDCore DDRec REQUIRED)

gaudi_depends_on_subdirs(
  Detector/DetInterface
  Utilities/KalTest
)

get_target_property(to_incl KalTestLib SOURCE_DIR)
if (to_incl)
  LIST( APPEND DICT_INCLUDE_DIRS ${to_incl})
else()
  message(FATAL_ERROR "Failed to get the source dir for package KalTestLib")
endif()

set( DICT_CINT_DEFINITIONS "HANDLE_DICT_EXCEPTIONS=IGNORED_FOR_CINT" )
set( DICT_INPUT_DIRS gen kern lctpc/gearTPC )
set( DICT_OUTPUT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/rootdict )

foreach( DICT_DIR ${DICT_INPUT_DIRS} )
  list( APPEND DICT_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src/${DICT_DIR} )
  aux_source_directory( src/${DICT_DIR} LIB_SOURCES )
endforeach()

foreach( DICT_DIR ${DICT_INPUT_DIRS} )
  file( GLOB DICT_HEADER_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/${DICT_DIR}/*.h" )
  file( GLOB LINKDEFS "${CMAKE_CURRENT_SOURCE_DIR}/src/${DICT_DIR}/LinkDef.h" )

  if( LINKDEFS )
    list( REMOVE_ITEM DICT_HEADER_FILES "${LINKDEFS}" )
  endif()

  string(REGEX REPLACE "/" "_" DICT ${DICT_DIR})    

  add_custom_command(OUTPUT ${DICT_OUTPUT_DIR}/${DICT}.cxx ${DICT_OUTPUT_DIR}/${DICT}_rdict.pcm
  	             COMMAND mkdir -p ${DICT_OUTPUT_DIR}   
                     COMMAND ${ROOT_rootcling_CMD}
                     ARGS -f ${DICT_OUTPUT_DIR}/${DICT}.cxx -s ${DICT_OUTPUT_DIR}/${DICT} -inlineInputHeader
                     "$<$<BOOL:$<JOIN:${DICT_CINT_DEFINITIONS},>>:-D$<JOIN:${DICT_CINT_DEFINITIONS},;-D>>"
                     "$<$<BOOL:$<JOIN:${DICT_INCLUDE_DIRS},>>:-I$<JOIN:${DICT_INCLUDE_DIRS},;-I>>"
                     "$<JOIN:${DICT_HEADER_FILES},;>" "$<JOIN:${LINKDEFS},;>"
                     DEPENDS ${DICT_HEADER_FILES} ${LINKDEFS}
                     COMMAND_EXPAND_LISTS
                     )

  add_custom_target( ${DICT} DEPENDS ${DICT_OUTPUT_DIR}/${DICT}.cxx ${DICT_OUTPUT_DIR}/${DICT}_rdict.pcm )

  set_source_files_properties(${DICT_OUTPUT_DIR}/${DICT}.cxx ${DICT_OUTPUT_DIR}/${DICT}_rdict.pcm
                              PROPERTIES
                              GENERATED TRUE
                              COMPILE_FLAGS "-Wno-unused-function -Wno-overlength-strings"
                              )

  install(FILES ${DICT_OUTPUT_DIR}/${DICT}_rdict.pcm DESTINATION lib)

  list( APPEND LIB_SOURCES ${DICT_OUTPUT_DIR}/${DICT}.cxx )
endforeach()

FILE( GLOB_RECURSE COMMON_SOURCES "src/ild/*.cc" )
FILE( GLOB_RECURSE lib_headers "src/ild/*.h" )

include_directories( src/ild/common )
include_directories( ${DICT_INCLUDE_DIRS} )

set( KalDetLib_srcs ${LIB_SOURCES} ${COMMON_SOURCES} )

gaudi_add_library(KalDetLib ${KalDetLib_srcs}
		 PUBLIC_HEADERS kaldet
                 LINK_LIBRARIES GaudiKernel ROOT ${CLHEP_LIBRARIES} LCIO ${GEAR_LIBRARIES} KalTestLib EDM4HEP::edm4hep EDM4HEP::edm4hepDict ${DD4hep_COMPONENT_LIBRARIES}
)
