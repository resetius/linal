cmake_minimum_required(VERSION 2.8.0)
project(linal CXX C)

set (LINAL_INCLUDE "" PARENT_SCOPE)
set (LINAL_DEFINES "" PARENT_SCOPE)
set (LINAL_LIBS "" PARENT_SCOPE)

list(APPEND LINAL_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/src)

set (CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(Common)
find_package(OpenMP)
find_package(SelfUmfPack)
find_package(UmfPack)
find_package(SelfSuperLU)
find_package(CUDA)
include(CommonDefinitions)

set (LINAL_LIBS ${LINAL_LIBS} PARENT_SCOPE)
set (LINAL_CUDA_LIBS ${LINAL_CUDA_LIBS} PARENT_SCOPE)
set (LINAL_DEFINES ${LINAL_DEFINES} PARENT_SCOPE)
set (LINAL_INCLUDE ${LINAL_INCLUDE} PARENT_SCOPE)

if (EXISTS "${CMAKE_CURRENT_BINARY_DIR}/settings.cmake")
	include(${CMAKE_CURRENT_BINARY_DIR}/settings.cmake)
endif()

add_subdirectory(contrib)
add_subdirectory(src)
add_subdirectory(tests)

