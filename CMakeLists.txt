cmake_minimum_required(VERSION 2.8)
  
# Project
project (SearchLab)

# thanks to https://github.com/nathan-osman/CXX11-CMake-Macros
include(${CMAKE_SOURCE_DIR}/CXX11.cmake)

check_for_cxx11_compiler(CXX11_COMPILER)

if(CXX11_COMPILER)
# C++11 compiler available
  include_directories( ${CMAKE_SOURCE_DIR}/include )

  set_source_files_properties("*.c" PROPERTIES LANGUAGE C )
  set_source_files_properties("*.cc" PROPERTIES LANGUAGE CXX )
  set_source_files_properties("*.cpp" PROPERTIES LANGUAGE CXX )
  set_source_files_properties("*.h" PROPERTIES LANGUAGE CXX )
  
  # Kratos
  include_directories("${KRATOS_PATH}/kratos/")
  include_directories("${KRATOS_PATH}/kratos/includes")
  
  # add the executable
  add_executable(SearchLab.exe SearchLab.cpp)
    
  # require a C++11 compiler
  # this seems not to work ...
  set_property(TARGET SearchLab.exe PROPERTY CXX_STANDARD 11)
  set_property(TARGET SearchLab.exe PROPERTY CXX_STANDARD_REQUIRED ON)
  # does not work:
  # target_compile_features(SearchLab PUBLIC cxx_constexpr)
endif()