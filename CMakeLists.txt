cmake_minimum_required(VERSION 3.2)

# Project
project (SarchLab)

find_package(MPI REQUIRED)

include_directories( ${CMAKE_SOURCE_DIR}/include )
include_directories( ${MPI_INCLUDE_PATH} )

set_source_files_properties("*.c" PROPERTIES LANGUAGE CXX )
set_source_files_properties("*.h" PROPERTIES LANGUAGE CXX )

# Kratos
include_directories("${KRATOS_PATH}/kratos/")
include_directories("${KRATOS_PATH}/kratos/includes")

# add the executable
add_executable(SearchLabMpi SearchLabMpi.cpp)
add_executable(SearchLab SearchLab.cpp)

message("Mpi used: ${MPI_LIBRARIES}")

target_link_libraries(SearchLabMpi ${MPI_LIBRARIES})
target_link_libraries(SearchLab)