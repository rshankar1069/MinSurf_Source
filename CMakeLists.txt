
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

project(minSurf)

find_package(VTK COMPONENTS 
  vtkCommonCore
  vtkCommonDataModel
  vtkIOXML QUIET)

if (NOT VTK_FOUND)
  message("Skipping XMLStructuredGridWriter: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()

#do not use this, because CMAKE_BUILD_TPYE can be set to Release or Debug add_compile_options(-g -O3 -std=c++11)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIERD ON)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")

if (VTK_VERSION VERSION_LESS "8.90.0")
  # old system
  include(${VTK_USE_FILE})
  add_executable(minSurf MACOSX_BUNDLE minSurf.cpp )
  target_link_libraries(minSurf PRIVATE ${VTK_LIBRARIES})
else ()
  # include all components
  add_executable(minSurf MACOSX_BUNDLE minSurf.cpp )
  target_link_libraries(minSurf PRIVATE ${VTK_LIBRARIES})
  # vtk_module_autoinit is needed
  vtk_module_autoinit(
    TARGETS minSurf
    MODULES ${VTK_LIBRARIES}
    )
endif () 
