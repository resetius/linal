cmake_minimum_required(VERSION 2.8.0)
project(linal CXX C)

set (LINAL_HOME ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
set (CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(Common)
find_package(OpenMP)
find_package(SelfUmfPack)
find_package(UmfPack)
find_package(SelfSuperLU)
find_package(CUDA)
include(CommonDefinitions)

if (EXISTS "${CMAKE_CURRENT_BINARY_DIR}/settings.cmake")
	include(${CMAKE_CURRENT_BINARY_DIR}/settings.cmake)
endif()

add_subdirectory(contrib)
add_subdirectory(src)
