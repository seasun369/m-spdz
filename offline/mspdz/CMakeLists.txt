cmake_minimum_required (VERSION 3.0)
project (mspdz)
set(NAME "mspdz")

find_path(CMAKE_FOLDER NAMES cmake/emp-tool-config.cmake)
include(${CMAKE_FOLDER}/cmake/emp-base.cmake)

find_package(emp-ot REQUIRED)
include_directories(${EMP-OT_INCLUDE_DIRS})

ENABLE_TESTING()
ADD_SUBDIRECTORY(vole)
