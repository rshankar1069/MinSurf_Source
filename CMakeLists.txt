
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

project(minSurf)

set(CXX_STANDARD 11)
include_directories(./dco_cpp/include)

## Find VTK
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
set(CMAKE_CXX_STANDARD_REQUIRED ON)

## OpenMP options
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


## Link VTK
message (STATUS "VTK_VERSION: ${VTK_VERSION}")

if (VTK_VERSION VERSION_LESS "8.90.0")
  # old system
  include(${VTK_USE_FILE})
  add_executable(minSurf minSurf.cpp )
  target_link_libraries(minSurf PRIVATE ${VTK_LIBRARIES} ${CMAKE_SOURCE_DIR}/dco_cpp/lib/libdcoc.a)
else ()
  # include all components
  add_executable(minSurf minSurf.cpp )
  target_link_libraries(minSurf PRIVATE ${VTK_LIBRARIES} ${CMAKE_SOURCE_DIR}/dco_cpp/lib/libdcoc.a)
  # vtk_module_autoinit is needed
  vtk_module_autoinit(
    TARGETS minSurf
    MODULES ${VTK_LIBRARIES}
    )
endif () 

set(CMAKE_BUILD_TYPE Release)
include_directories(.)
############ GTEST
# Download and unpack googletest at configure time
configure_file(${CMAKE_SOURCE_DIR}/additionals/gtest.cmake googletest-download/CMakeLists.txt)

execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
    ${CMAKE_BINARY_DIR}/googletest-build
    EXCLUDE_FROM_ALL)

########## END GTEST


include(${CMAKE_SOURCE_DIR}/additionals/cmake_scripts/stce_utils.cmake)


#add tests target for recursively adding subfolders
stce_add_recursively(${CMAKE_SOURCE_DIR}/tests/)









