
include_directories(${Surfpack_BINARY_DIR}
                    ${Surfpack_SOURCE_DIR}/src)

set(local_library surfpack)

## Top-level Surfpack library sources.
set(lib${local_library}_sources
   AxesBounds.cpp
   AxesBounds.h
   MersenneTwister.h
   ModelFactory.cpp
   ModelFactory.h
   SurfData.cpp
   SurfData.h
   surfpack.cpp
   surfpack.h
   SurfpackMatrix.h
   SurfpackMatrix.cpp
   surfpack_system_headers.h
   SurfpackInterface.h
   SurfpackInterface.cpp
   SurfpackParserArgs.h
   SurfpackParserArgs.cpp
   SurfPoint.cpp
   SurfPoint.h
   Conmin.cpp
   Conmin.h
   ModelFitness.cpp
   ModelFitness.h
)

# Build a single libsurfpack from all subdirectories
set(local_src_subdirs mars surfaces surfaces/nkm)

# Aggregate all the sources into a single library for installation
# BUT split the fortran out into its own library (necessary for native Windows)
foreach(sd_it ${local_src_subdirs})

  add_subdirectory(${sd_it})

  # grab the trailing part of any paths
  string(REGEX MATCH "[^/]+$" last_path_name ${sd_it})

  # extract the variables ${last_path_name}_sources from the subdir's env
  get_directory_property(sd_sources DIRECTORY ${sd_it} 
    DEFINITION ${last_path_name}_sources)

  # add the sources to the library's list and includes
  foreach ( src_file ${sd_sources} )

    if(${src_file} MATCHES "\\.f$" OR ${src_file} MATCHES "\\.f..$")
      list ( APPEND lib${local_library}_sources_fortran "${sd_it}/${src_file}" )
    else()
      list ( APPEND lib${local_library}_sources "${sd_it}/${src_file}" )
    endif()
  endforeach()

  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/${sd_it})

endforeach()

# the interpreter should not be linked into the main library, just
# into executables
if(SURFPACK_STANDALONE)
  add_subdirectory(interpreter)
endif()

# WJB - ToDo: look into CMake template support features: SurfpackMatrix.tpl

add_library(${local_library} ${lib${local_library}_sources})
add_library(${local_library}_fortran ${lib${local_library}_sources_fortran})
target_link_libraries(${local_library} ${local_library}_fortran)
#if (HAVE_PECOS)
#  target_link_libraries(${local_library} pecos ${Pecos_LINK_LIBRARIES} 
#    ${Teuchos_LIBRARIES} ${TEUCHOS_LIBRARY} ${LAPACK_LIBS} ${BLAS_LIBS})
#endif()

#target_link_libraries(surfpack_src mars teuchos conmin ${LAPACK_LIBRARIES})
install(TARGETS ${local_library} EXPORT ${ExportTarget} DESTINATION lib)
install(TARGETS ${local_library}_fortran EXPORT ${ExportTarget} DESTINATION lib)

install(FILES SurfpackMatrix.h surfpack_system_headers.h DESTINATION include)


# Surfpack C interface, compiled into separate library for now, though
# no reason couldn't be directly in libsurfpack.
add_library(surfpack_c_interface surfpack_c_interface.cpp)
install(TARGETS surfpack_c_interface DESTINATION lib)
install(FILES surfpack_c_interface.h DESTINATION include)