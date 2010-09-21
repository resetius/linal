cmake_minimum_required(VERSION 2.8.0)
project(linal CXX C)

set (LINAL_INCLUDE "")
set (LINAL_DEFINES "")
set (LINAL_LIBS "")

list(APPEND LINAL_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/src)

set (CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if (EXISTS "${CMAKE_CURRENT_BINARY_DIR}/settings.cmake")
	include(${CMAKE_CURRENT_BINARY_DIR}/settings.cmake)
endif()

include(Common)
find_package(OpenMP)
find_package(SelfUmfPack)
find_package(UmfPack)
find_package(SelfSuperLU)
find_package(CUDA)
# set(BLA_VENDOR Intel10_64lp)
# set(BLA_VENDOR ACML)
find_package(BLAS)
include(CommonDefinitions)

if (NOT CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR) # <- NOT root project
set (LINAL_LIBS ${LINAL_LIBS} PARENT_SCOPE)
set (LINAL_CUDA_LIBS ${LINAL_CUDA_LIBS} PARENT_SCOPE)
set (LINAL_DEFINES ${LINAL_DEFINES} PARENT_SCOPE)
set (LINAL_INCLUDE ${LINAL_INCLUDE} PARENT_SCOPE)
endif ()

if (BLAS_FOUND)
	message(STATUS "BLAS: ${BLAS_LIBRARIES}")
endif ()

add_subdirectory(contrib)
add_subdirectory(src)
add_subdirectory(tests)

