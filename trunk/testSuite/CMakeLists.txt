add_executable(copasiAPItest copasiAPItest.c)
target_link_libraries(copasiAPItest copasi)

add_executable(sbwAPItest SBWAPItest.c)
target_link_libraries(sbwAPItest copasiSBW)

#add_subdirectory(subdir)

