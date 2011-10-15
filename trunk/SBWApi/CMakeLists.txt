INCLUDE_DIRECTORIES( 
	BEFORE
	${CMAKE_CURRENT_SOURCE_DIR}
)
ADD_EXECUTABLE(copasi_sbw_api_test test_sbw_api.c)
TARGET_LINK_LIBRARIES( copasi_sbw_api_test copasi)
