
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

add_compile_options(-O3)
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
