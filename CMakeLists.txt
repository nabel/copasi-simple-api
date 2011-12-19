PROJECT( LIBCOPASI )
CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )

#--------------------------------------------------------------------------
#    just some setup for compilers etc.
#--------------------------------------------------------------------------
#just a few linker frags needed for Linux
IF (UNIX AND NOT APPLE)
	SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -lrt")
	SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -lrt")
	SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lrt")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
ENDIF(UNIX AND NOT APPLE)

SET (LIBRARY_OUTPUT_PATH ${LIBCOPASI_BINARY_DIR})
SET (EXECUTABLE_OUTPUT_PATH ${LIBCOPASI_BINARY_DIR})

# CMake not supporting Win32 MinGW?
IF( WIN32 AND MINGW)
  ADD_DEFINITIONS(-DMINGW)
  SET(CMAKE_RC_COMPILER_INIT windres)
  ENABLE_LANGUAGE(RC)
  SET(CMAKE_RC_COMPILE_OBJECT
    "<CMAKE_RC_COMPILER> <FLAGS> -O coff <DEFINES> -i <SOURCE> -o <OBJECT>")
ENDIF( WIN32 AND MINGW)

#Build type
OPTION( USE_DEBUG_MODE "Release vs Debug mode" OFF)
IF (USE_DEBUG_MODE)
	SET(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type")
	SET(CMAKE_RC_COMPILE_OBJECT Debug)
ELSE(USE_DEBUG_MODE)
	SET(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type")
	SET(CMAKE_RC_COMPILE_OBJECT Release)
ENDIF(USE_DEBUG_MODE)

#------------------------------------------------------------------------------
#  cmake definitions for curl, sbml, and some other packages
#------------------------------------------------------------------------------
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


#--------------------------------------------------------------------------
# include folders
#--------------------------------------------------------------------------
INCLUDE_DIRECTORIES( BEFORE
    ${LIBCOPASI_SOURCE_DIR}
    ${LIBCOPASI_SOURCE_DIR}/curl/include
    ${LIBCOPASI_SOURCE_DIR}/expat
	${LIBCOPASI_SOURCE_DIR}/expat/lib
    ${LIBCOPASI_SOURCE_DIR}/muparser
    ${LIBCOPASI_SOURCE_DIR}/lapack/INCLUDE
	${LIBCOPASI_SOURCE_DIR}/libstructural/include
)

#------------------------------------------------------------------------------
#build individual subdirectories using CMakeLists.txt in each directory
#------------------------------------------------------------------------------

#if the system already has lapack, then don't build it
OPTION( BUILD_LAPACK "check this if your system does not have lapack library" ON ) #use cmake gui to set this
IF (BUILD_LAPACK)
	ADD_SUBDIRECTORY( lapack )     #linear algebra
ENDIF(BUILD_LAPACK)

ADD_SUBDIRECTORY( muparser )  #math parser
ADD_SUBDIRECTORY( curl )      #needed by raptor
ADD_SUBDIRECTORY( raptor )    #RDF parser
ADD_SUBDIRECTORY( sbml )    #systems biology markup language parser
ADD_SUBDIRECTORY( antimony )   #human-readable test-based reaction networks
ADD_SUBDIRECTORY( boost )      #regular expressions
ADD_SUBDIRECTORY( copasi )     #simulation and analysis package
ADD_SUBDIRECTORY( libstructural )   #network structural analysis based on stoichiometry matrix
ADD_SUBDIRECTORY( galib )   #genetic algorithms library

OPTION( TEST_COPASI_API "test copasi api" ON )

IF (TEST_COPASI_API)
	ADD_SUBDIRECTORY( testSuite )     #linear algebra
ENDIF(TEST_COPASI_API)

ADD_SUBDIRECTORY(python)  #python wrapper
ADD_SUBDIRECTORY(octave)  #octave wrapper

