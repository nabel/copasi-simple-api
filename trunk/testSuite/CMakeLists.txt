add_executable(testBasicAPI testBasicAPI.c)
target_link_libraries(testBasicAPI copasi_api)

add_executable(testSBWAPI testSBWAPI.c)
target_link_libraries(testSBWAPI copasiSBW)

add_executable(testSBMLread testSBMLread.c)
target_link_libraries(testSBMLread copasi_api)

add_executable(testOptimizer testOptimizer.c)
target_link_libraries(testOptimizer copasi_api)

#add_subdirectory(subdir)

