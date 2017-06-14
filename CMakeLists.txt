cmake_minimum_required(VERSION 3.2)

set_source_files_properties("*.c" PROPERTIES LANGUAGE CXX )
set_source_files_properties("*.h" PROPERTIES LANGUAGE CXX )

# Project
project (SarchLab)

# Kratos
include_directories("${KRATOS_PATH}/kratos/")
include_directories("${KRATOS_PATH}/kratos/includes")

# add the executable
add_executable(SearchLab SearchLab.cpp)
