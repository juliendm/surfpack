include_directories("${Surfpack_SOURCE_DIR}/src" 
  "${Surfpack_SOURCE_DIR}/src/interpreter"
  )

link_directories(${Surfpack_BINARY_DIR}/src/interpreter)
add_executable(surfpack_exe main.cpp)
target_link_libraries(surfpack_exe surfpack_interpreter ${SURFPACK_LIBS} 
  ${SURFPACK_TPL_LIBS} ${SURFPACK_SYSTEM_LIBS}
  )
set_target_properties(surfpack_exe PROPERTIES OUTPUT_NAME surfpack) 

install(TARGETS surfpack_exe DESTINATION bin)
