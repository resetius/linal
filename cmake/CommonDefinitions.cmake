if (NOT _LINAL_COMMON_DEFFINITION_)
set (_LINAL_COMMON_DEFFINITION_ TRUE)

if (MSVC)
        add_definitions(/D_CRT_SECURE_NO_WARNINGS)
endif (MSVC)

if (OPENMP_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
add_definitions(${LINAL_DEFINES})
include_directories(${LINAL_INCLUDE})

if (CUDA_FOUND)
	set (LINAL_CUDA_LIBS ${LINAL_LIBS})
	list(APPEND LINAL_CUDA_LIBS linal_cu)
endif ()

list(APPEND LINAL_LIBS linal)
list(APPEND LINAL_DEFINES -DSPARSE)

if (BLAS_FOUND)
	list(APPEND LINAL_LIBS ${BLAS_LIBRARIES})
endif ()

endif ()

