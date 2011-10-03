#---------------------
#    OS STUFF
#---------------------
IF (WIN32)

   SET(BUILD_OS WIN32)

ELSE(WIN32)	
	
	IF (UNIX)
		IF (APPLE)
			SET(BUILD_OS MAC)
		ELSEIF(APPLE)
			SET(BUILD_OS LINUX)
		ENDIF(APPLE)
	ENDIF(UNIX)

ENDIF(WIN32)

#-----------------------------------
#    INCLUDE DIRECTORIES
#-----------------------------------

INCLUDE_DIRECTORIES( 
	BEFORE
	${LIBCOPASI_SOURCE_DIR}/libSBML/include
	${LIBCOPASI_SOURCE_DIR}/libSBML/src
	${LIBCOPASI_SOURCE_DIR}/copasi
	${LIBCOPASI_SOURCE_DIR}/raptor/src
	${LIBCOPASI_SOURCE_DIR}/raptor/librdfa
	${LIBCOPASI_SOURCE_DIR}/raptor/util
	${LIBCOPASI_SOURCE_DIR}/antimony
	${LIBCOPASI_SOURCE_DIR}/galib/ga
	${LIBCOPASI_SOURCE_DIR}/galib
	${LIBCOPASI_SOURCE_DIR}/libstructural/include/
)

#---------------------
#    DEFINITIONS
#---------------------

#ADD_DEFINITIONS(-DCOPASI_MAIN)
ADD_DEFINITIONS(-DWITH_CSPMETHOD)
#ADD_DEFINITIONS(-D__cplusplus)
ADD_DEFINITIONS(-DUSE_CLAPACK)
#ADD_DEFINITIONS(-DUSE_LAPACK)
#ADD_DEFINITIONS(-DCOPASI_DEBUG)
#ADD_DEFINITIONS(-DDEBUG_MATRIX)
ADD_DEFINITIONS(-DTC_EXPORTS)

#------------------------------------------
#  DEFINITIONS FOR WINDOWS
#------------------------------------------

IF (WIN32)
	ADD_DEFINITIONS(-DWIN32)
	ADD_DEFINITIONS(-D_WIN32)
	ADD_DEFINITIONS(-D_WINDOWS)
	ADD_DEFINITIONS(-D_USRDLL)
	ADD_DEFINITIONS(-D_EXPORTS)
	ADD_DEFINITIONS(-DLIBSBML_EXPORTS)
	#ADD_DEFINITIONS(-DIN_LIBXML)
	IF (MINGW)
		ADD_DEFINITIONS(-DCYGWIN)
		ADD_DEFINITIONS(-D__CYGWIN__)
		IF (USE_STDCALL)
			ADD_DEFINITIONS(-DUSE_STDCALL)
		 ENDIF(USE_STDCALL)
	ELSE(MINGW)
		ADD_DEFINITIONS(-Dstrcasecmp=_stricmp)
		ADD_DEFINITIONS(-Dstrncasecmp=_strnicmp)
	ENDIF(MINGW)
ENDIF (WIN32)
 IF (APPLE)
  	 ADD_DEFINITIONS(-D_APPLE)
ENDIF(APPLE)
#--------------------------
#  COPASI SOURCE FILES
#--------------------------

FILE( GLOB copasi_SRC_FILES
	${LIBCOPASI_SOURCE_DIR}/copasi/model/*.cpp
	${LIBCOPASI_SOURCE_DIR}/copasi/compareExpressions/*.cpp
	${LIBCOPASI_SOURCE_DIR}/copasi/CopasiDataModel/*.cpp
	${LIBCOPASI_SOURCE_DIR}/copasi/crosssection/*.cpp
	${LIBCOPASI_SOURCE_DIR}/copasi/elementaryFluxModes/*.cpp
	${LIBCOPASI_SOURCE_DIR}/copasi/function/*.cpp
	${LIBCOPASI_SOURCE_DIR}/copasi/layout/*.cpp
	${LIBCOPASI_SOURCE_DIR}/copasi/libs/*.cpp
	${LIBCOPASI_SOURCE_DIR}/copasi/lyap/*.cpp
	${LIBCOPASI_SOURCE_DIR}/copasi/moieties/*.cpp
	${LIBCOPASI_SOURCE_DIR}/copasi/odepack++/*.cpp
	${LIBCOPASI_SOURCE_DIR}/copasi/optimization/*.cpp
	${LIBCOPASI_SOURCE_DIR}/copasi/parameterFitting/*.cpp
	${LIBCOPASI_SOURCE_DIR}/copasi/plot/*.cpp
	${LIBCOPASI_SOURCE_DIR}/copasi/randomGenerator/*.cpp	
	${LIBCOPASI_SOURCE_DIR}/copasi/sbml/*.cpp
	${LIBCOPASI_SOURCE_DIR}/copasi/sbmlunit/*.cpp
	${LIBCOPASI_SOURCE_DIR}/copasi/scan/*.cpp
	${LIBCOPASI_SOURCE_DIR}/copasi/sensitivities/*.cpp
	${LIBCOPASI_SOURCE_DIR}/copasi/steadystate/*.cpp
	${LIBCOPASI_SOURCE_DIR}/copasi/report/*.cpp
	${LIBCOPASI_SOURCE_DIR}/copasi/trajectory/*.cpp
	${LIBCOPASI_SOURCE_DIR}/copasi/tss/*.cpp
	${LIBCOPASI_SOURCE_DIR}/copasi/tssanalysis/*.cpp
	${LIBCOPASI_SOURCE_DIR}/copasi/utilities/*.cpp
	${LIBCOPASI_SOURCE_DIR}/copasi/xml/*.cpp
	${LIBCOPASI_SOURCE_DIR}/copasi/commandline/*.cpp
	${LIBCOPASI_SOURCE_DIR}/copasi/MIRIAM/*.cpp
	${LIBCOPASI_SOURCE_DIR}/copasi/MIRIAM/WebServicesIssues/*.cpp
	${LIBCOPASI_SOURCE_DIR}/copasi/copasi_api.cpp
	${LIBCOPASI_SOURCE_DIR}/TC_structs.c
	${LIBCOPASI_SOURCE_DIR}/mtrand.c
)

ADD_DEFINITIONS(-DHAVE_EXPAT_CONFIG_H)
FILE( GLOB  EXPAT_SRC_FILES   ${LIBCOPASI_SOURCE_DIR}/expat/lib/*.c)

#--------------------------
#  MAKE libcopasi
#--------------------------

ADD_LIBRARY( copasi
#  STATIC
  SHARED
  ${copasi_SRC_FILES}
  ${EXPAT_SRC_FILES}
)

IF (WIN32 AND NOT MINGW)
	IF(CMAKE_BUILG_TYPE STREQUAL "Debug")
		SET( boost_library_name libboost_regex-vc100-mt-1_48)
	ELSEIF(CMAKE_BUILG_TYPE STREQUAL "Debug")
		SET( boost_library_name libboost_regex-vc100-mt-gd-1_48)
	ENDIF(CMAKE_BUILG_TYPE STREQUAL "Debug")
ELSE (WIN32 AND NOT MINGW)
	SET( boost_library_name boost_regex)
ENDIF (WIN32 AND NOT MINGW)

TARGET_LINK_LIBRARIES( copasi
  copasi
  sbml
  raptor
  lapack
  blas
  f2c
  antimony
  galib
  structural
  muparser
  ${boost_library_name}
)

#----------------------------------
#  INSTALL libcopasi in . (CPack)
#----------------------------------

IF ( WIN32 )
   INSTALL(TARGETS copasi DESTINATION bin/../)
ELSE ( WIN32 )
   INSTALL(TARGETS copasi LIBRARY DESTINATION bin/../)
ENDIF( WIN32 )

#----------------------------------
#  COPASI api
#----------------------------------

OPTION( TEST_COPASI_API "test copasi api" ON )

IF (TEST_COPASI_API)
  ADD_EXECUTABLE(copasi_api_test test_copasi.c)
  TARGET_LINK_LIBRARIES( copasi_api_test copasi)
ENDIF (TEST_COPASI_API)

