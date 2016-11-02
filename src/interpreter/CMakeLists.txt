
include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

# In specification maintenance mode only, we build the parser files
# into the source tree to faciliate native Windows builds. Header
# surfparse.hpp gets copied to the build tree during make.
if(ENABLE_SPEC_MAINT)
  bison_target(SurfParser surfparse.ypp 
    ${CMAKE_CURRENT_SOURCE_DIR}/surfparse.cpp)
  flex_target(SurfScanner lexer.l  ${CMAKE_CURRENT_SOURCE_DIR}/lexer.c)
  add_flex_bison_dependency(SurfScanner SurfParser)
endif()

set(surfparse_hpp "${CMAKE_CURRENT_SOURCE_DIR}/surfparse.hpp")
set(surfparse_h "${CMAKE_CURRENT_BINARY_DIR}/surfparse.h")

add_custom_command(
  DEPENDS "${surfparse_hpp}"
  OUTPUT  "${surfparse_h}"
  COMMAND ${CMAKE_COMMAND} -E copy "${surfparse_hpp}" "${surfparse_h}"
)

set(surfpack_interpreter_sources
  ${CMAKE_CURRENT_BINARY_DIR}/surfparse.h
  lexer.c
  surfparse.cpp
  SurfpackParser.cpp
  SurfpackParser.h
  FlexWrapper.cpp
  FlexWrapper.h
  SurfpackInterpreter.cpp
  SurfpackInterpreter.h
)

add_library(surfpack_interpreter ${surfpack_interpreter_sources})
install(TARGETS surfpack_interpreter EXPORT ${ExportTarget} DESTINATION lib)

