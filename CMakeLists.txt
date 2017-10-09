cmake_minimum_required(VERSION 3.2)

# Project
project (SarchLab)

if(DEFINED USE_MPI)
  find_package(MPI REQUIRED)
endif(DEFINED USE_MPI)

include_directories( ${CMAKE_SOURCE_DIR}/include )
include_directories( ${MPI_INCLUDE_PATH} )

set_source_files_properties("*.c" PROPERTIES LANGUAGE CXX )
set_source_files_properties("*.h" PROPERTIES LANGUAGE CXX )

# Kratos
if(DEFINED USE_KRATOS)
  include_directories("${KRATOS_PATH}/kratos/")
  include_directories("${KRATOS_PATH}/kratos/includes")
endif(DEFINED USE_KRATOS)

# add the executable
add_executable(SearchLab SearchLab.cpp)
target_link_libraries(SearchLab)

if(DEFINED USE_MPI)
  message("Mpi used: ${MPI_LIBRARIES}")
  add_executable(SearchLabMpi SearchLabMpi.cpp)
  target_link_libraries(SearchLabMpi ${MPI_LIBRARIES})
endif(DEFINED USE_MPI)

