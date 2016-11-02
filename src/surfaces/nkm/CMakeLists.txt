## nkm directory sources.
set(nkm_sources
  NKM_KrigingModel.cpp
  NKM_Optimize.cpp
  NKM_SurfData.cpp
  NKM_SurfPack.cpp
  NKM_SurfPackModel.cpp
  NKM_SurfMat.cpp
  NKM_PivotChol.f
)

# Optionally build the ValidateMain binary
if (Surfpack_NKM_Tests)
  link_directories(${Surfpack_BINARY_DIR}/src/interpreter)
  add_executable(ValidateMain NKM_ValidateMain.cpp)
  target_link_libraries(ValidateMain surfpack_interpreter ${SURFPACK_LIBS} ${SURFPACK_TPL_LIBS} 
    ${SURFPACK_SYSTEM_LIBS}
    )
endif()
