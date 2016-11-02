cmake_minimum_required(VERSION 2.8)

project("SurfpackCInterfaceDemo" C)

# Set these to the installed location of a surfpack distro to link against it
#set(Surfpack_INCLUDE_DIR "/apps/surfpack/trunk/include")
#set(Surfpack_LIBRARY_DIR "/apps/surfpack/trunk/lib")
#set(Boost_LIBRARY_DIR "/apps/boost/1.49/lib")
#set(LinAlg_LIBRARY_DIR "/usr/lib64")

set(SURFPACK_LIBS "surfpack;surfpack_fortran")
set(SURFPACK_TPL_LIBS "ncsuopt;conmin")
set(SURFPACK_SYSTEM_LIBS "${Boost_LIBRARIES};${LAPACK_LIBS};${BLAS_LIBS}")

link_directories("${Boost_LIBRARY_DIR}")
list(APPEND SURFPACK_SYSTEM_LIBS "boost_serialization")

link_directories("${LinAlg_LIBRARY_DIR}")
list(APPEND SURFPACK_SYSTEM_LIBS "lapack;blas")

list(APPEND SURFPACK_SYSTEM_LIBS "stdc++")

include_directories("${Surfpack_INCLUDE_DIR}")
link_directories("${Surfpack_LIBRARY_DIR}")

add_executable(eval_model eval_model.c)
target_link_libraries(eval_model surfpack_c_interface 
  ${SURFPACK_LIBS} ${SURFPACK_TPL_LIBS}
  ${SURFPACK_SYSTEM_LIBS}
  )
