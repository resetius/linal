if (MSVC)
        add_definitions(/D_CRT_SECURE_NO_WARNINGS)
endif (MSVC)

add_definitions(${OpenMP_C_FLAGS})
add_definitions(${LINAL_DEFINES})
include_directories(${LINAL_INCLUDE})

if (CUDA_FOUND)
	set (LINAL_CUDA_LIBS ${LINAL_LIBS})
	list(APPEND LINAL_CUDA_LIBS linal_cu)
endif ()

list(APPEND LINAL_LIBS linal)
list(APPEND LINAL_DEFINES -DSPARSE)

