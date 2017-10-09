cmake_minimum_required(VERSION 3.2)

# Project
project (SarchLab)

if(USE_MPI MATCHES ON)
  find_package(MPI REQUIRED)
endif(USE_MPI MATCHES ON)

include_directories( ${CMAKE_SOURCE_DIR}/include )
include_directories( ${MPI_INCLUDE_PATH} )

set_source_files_properties("*.c" PROPERTIES LANGUAGE CXX )
set_source_files_properties("*.h" PROPERTIES LANGUAGE CXX )

# Kratos
if(USE_KRATOS MATCHES ON)
  include_directories("${KRATOS_PATH}/kratos/")
  include_directories("${KRATOS_PATH}/kratos/includes")
endif(USE_KRATOS MATCHES ON)

# add the executable
add_executable(SearchLab SearchLab.cpp)
target_link_libraries(SearchLab)

if(USE_MPI MATCHES ON)
  message("Mpi used: ${MPI_LIBRARIES}")
  add_executable(SearchLabMpi SearchLabMpi.cpp)
  target_link_libraries(SearchLabMpi ${MPI_LIBRARIES})
endif(USE_MPI MATCHES ON)

