# set global variable <cmake_debug_messages> to ON for printing
function(debug_message message)
  if (cmake_debug_messages)
    message(${message})
  endif()
endfunction()

# add all sources in current directory for an executable and add
# dependency to target above current (either another executable, or
# just a dummy/custom target); do so for all subdirectories.
function(stce_add_recursively dir)
  debug_message("adding recursively: ${dir}")
  stce_get_name_and_category(name category ${dir})
  
  file(GLOB sources "${dir}/*.cpp")
  if (sources)
    debug_message("  adding executable ${name}")
    add_executable(${name} ${sources})
    target_link_libraries(${name} gtest_main ${DCO4_LIBS})
    if(${name} MATCHES "tests.*")
      add_test(${name} ${name})
    endif()
    #string(FIND ${name} "tests" pos1)
    #message("pos1='${pos1}'")
    #If ( ${pos1} != "-1")
    #  add_test(${name} ${name})
    #endif()
  else()
    debug_message("  adding custom target ${name}")
    add_custom_target(${name})
  endif()
  if (category)
    debug_message("  adding dependency ${category}<-${name}")
    add_dependencies(${category} ${name})
  endif()
  
  stce_get_subdirs(subdirs ${dir})
  foreach(subdir ${subdirs})
    if (EXISTS ${dir}/${subdir}/CMakeLists.txt)
      add_subdirectory(${dir}/${subdir})
    else()
      stce_add_recursively(${dir}/${subdir})
    endif()
  endforeach()

endfunction()

function(stce_get_subdirs result dir)
  file(GLOB children RELATIVE ${dir} ${dir}/*)
  set(dirlist "")
  foreach(child ${children})
    if(IS_DIRECTORY ${dir}/${child})
      list(APPEND dirlist ${child})
    endif()
  endforeach()
  set(${result} ${dirlist} PARENT_SCOPE)
endfunction()

function(stce_get_name_and_category name category dir)
  # remove trailing slash
  string(REGEX REPLACE "\n$" "" dir "${dir}")

  stce_get_name_from_file(local_name     ${dir}/dummy.cpp)
  stce_get_name_from_file(local_category ${dir})

  set(${name}     ${local_name}     PARENT_SCOPE)
  set(${category} ${local_category} PARENT_SCOPE)
endfunction()

# convert filename (file) to a test name by creating a relative path
# and replacing "/" with "."
function(stce_get_name_from_file name file)
  get_filename_component(cur_abs_path ${file} ABSOLUTE)
  get_filename_component(cur_dir ${cur_abs_path} DIRECTORY)

  string(REPLACE ${CMAKE_BINARY_DIR}/ "" cur_dir ${cur_dir})
  string(REPLACE ${CMAKE_BINARY_DIR}  "" cur_dir ${cur_dir})
  string(REPLACE ${CMAKE_SOURCE_DIR}/ "" cur_dir ${cur_dir})
  string(REPLACE ${CMAKE_SOURCE_DIR}  "" cur_dir ${cur_dir})

  if (cur_dir)
    string(REPLACE ${CMAKE_SOURCE_DIR}/ "" cur_name ${cur_dir})
    string(REPLACE "/" "." cur_name ${cur_name})
  endif()

  set(${name} ${cur_name} PARENT_SCOPE)
endfunction()


