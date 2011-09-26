SET(LIBRARY_OUTPUT_PATH ${EXECUTABLE_OUTPUT_PATH})

#----------- include folders ----------------------

INCLUDE_DIRECTORIES( 
	BEFORE
	${LIBCOPASI_SOURCE_DIR}/expat
	${LIBCOPASI_SOURCE_DIR}/expat/lib
	${LIBCOPASI_SOURCE_DIR}/libSBML/include
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/xml 
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/util 
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/validator 
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/validator/constraints/
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/units 
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/sbml
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/sbml/layout 
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/math 
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/common 
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/compress 
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/annotation
)

#------ DEFINITIONS --------------

#ADD_DEFINITIONS(-DLIBSBML_STATIC)
ADD_DEFINITIONS(-DLIBSBML_EXPORTS)
ADD_DEFINITIONS(-DLIBLAX_EXPORTS)

ADD_DEFINITIONS(-DPACKAGE_NAME=\"libSBML\")
ADD_DEFINITIONS(-DNDEBUG)
#ADD_DEFINITIONS(-DUSE_LIBXML)
ADD_DEFINITIONS(-DUSE_EXPAT)
ADD_DEFINITIONS(-DUSE_LAYOUT)

IF (WIN32)
	ADD_DEFINITIONS(-DWIN32)
	ADD_DEFINITIONS(-D_WIN32)
	ADD_DEFINITIONS(-D_WINDOWS)
	ADD_DEFINITIONS(-D_USRDLL)
	ADD_DEFINITIONS(-D_EXPORTS)
	#ADD_DEFINITIONS(-DIN_LIBXML)
	IF (MINGW)
		ADD_DEFINITIONS(-DCYGWIN)
	ENDIF(MINGW)
ENDIF (WIN32)

#------ SBML source files --------------

FILE( GLOB SBML_SRC_FILES
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/xml/XML*.cpp 
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/xml/XML*.c
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/xml/Expat*.cpp 
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/xml/Expat*.c
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/util/*.cpp 
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/util/*.c
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/validator/*.cpp 
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/validator/*.c
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/validator/constraints/*.cpp
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/validator/constraints/*.c
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/units/*.cpp 
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/units/*.c
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/math/*.cpp
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/math/*.c
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/annotation/*.cpp 
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/annotation/*.c
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/sbml/*.cpp
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/sbml/*.c
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/sbml/layout/*.cpp
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/sbml/layout/*.c
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/common/*.cpp
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/common/*.c
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/compress/CompressCommon.cpp
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/compress/InputDecompressor.cpp
	${LIBCOPASI_SOURCE_DIR}/libSBML/src/compress/OutputCompressor.cpp
)

ADD_DEFINITIONS(-DHAVE_EXPAT_CONFIG_H)
FILE( GLOB  EXPAT_SRC_FILES   ${LIBCOPASI_SOURCE_DIR}/expat/lib/*.c)

#------ make the library --------------

ADD_LIBRARY( sbml
  SHARED
  #STATIC
  ${SBML_SRC_FILES}
  ${EXPAT_SRC_FILES}
)

TARGET_LINK_LIBRARIES( sbml
  sbml
)

#------ install the library in the main TinkerCell folder --------------

IF ( WIN32 )
   INSTALL(TARGETS sbml DESTINATION bin/../)
ELSE ( WIN32 )
   INSTALL(TARGETS sbml LIBRARY DESTINATION bin/../)
ENDIF( WIN32 )

