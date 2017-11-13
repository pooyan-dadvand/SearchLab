cmake_minimum_required(VERSION 3.2)

# Project
project (SearchLab)

# thanks to https://github.com/nathan-osman/CXX11-CMake-Macros
include(${CMAKE_SOURCE_DIR}/CXX11.cmake)

check_for_cxx11_compiler(CXX11_COMPILER)

set( USE_MPI OFF CACHE BOOL "Use MPI implementation")

# for Windows
if (WIN32)
  if ( EXISTS C:/GID64-git/libZ/zlib-1.2.11)
    set( ZLIB_ROOT "C:/GID64-git/libZ/zlib-1.2.11")
  elseif ( EXISTS C:/GID64-git/libZ/windows)
    set( ZLIB_ROOT "C:/GID64-git/libZ/windows")
  else( EXISTS C:/GID64-git/libZ/windows)
    message( "ZLIB: trying to guess Zlib location")
  endif( EXISTS C:/GID64-git/libZ/zlib-1.2.11)
endif (WIN32)

find_package(ZLIB REQUIRED)

if(CXX11_COMPILER)
# C++11 compiler available

  if(USE_MPI MATCHES ON)
    find_package(MPI REQUIRED)
  endif(USE_MPI MATCHES ON)
  
  include_directories( ${CMAKE_SOURCE_DIR}/include )
  include_directories( ${MPI_INCLUDE_PATH} )
  include_directories( ${ZLIB_INCLUDE_DIRS} )
  
  set_source_files_properties("*.c" PROPERTIES LANGUAGE CXX )
  set_source_files_properties("*.h" PROPERTIES LANGUAGE CXX )

  # Warnings
  set ( WARN_FLAGS "")
  if ( CMAKE_COMPILER_IS_GNUCXX)
    set( WARN_FLAGS "-Wreturn-type -Wswitch -Wformat \
	-Wcast-align -Wconversion -Wcomment \
	-Wmissing-braces \
	-Wno-long-long  -Wno-variadic-macros \
	-Wno-char-subscripts -pedantic -Wextra -Wall -Wno-ignored-qualifiers")
  endif ( CMAKE_COMPILER_IS_GNUCXX)
  if ( MSVC)
    set( WARN_FLAGS "/W4")
  endif ( MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${WARN_FLAGS}")
  
  # Kratos
  if(USE_KRATOS MATCHES ON)
    include_directories("${KRATOS_PATH}/kratos/")
    include_directories("${KRATOS_PATH}/kratos/includes")
  endif(USE_KRATOS MATCHES ON)
  
  # add the executable
  add_executable(SearchLab.exe SearchLab.cpp)
  target_link_libraries(SearchLab.exe ${ZLIB_LIBRARIES})

  # require a C++11 compiler
  # this seems not to work ...
  set_property(TARGET SearchLab.exe PROPERTY CXX_STANDARD 11)
  set_property(TARGET SearchLab.exe PROPERTY CXX_STANDARD_REQUIRED ON)
  # does not work:
  # target_compile_features(SearchLab.exe PUBLIC cxx_constexpr)

  if(USE_MPI MATCHES ON)
    message("Mpi used: ${MPI_LIBRARIES}")
    add_executable(SearchLabMpi.exe SearchLabMpi.cpp)
    set_property(TARGET SearchLabMpi.exe PROPERTY CXX_STANDARD 11)
    set_property(TARGET SearchLabMpi.exe PROPERTY CXX_STANDARD_REQUIRED ON)
    target_link_libraries(SearchLabMpi.exe ${MPI_LIBRARIES})
  endif(USE_MPI MATCHES ON)
endif()
