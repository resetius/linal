project(linal_test CXX C)

set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
include_directories(${LINAL_INCLUDE})

add_executable(la_test_solver test_solver.cpp)
target_link_libraries(la_test_solver ${LINAL_LIBS})

if (OPENMP_FOUND)
	set_target_properties(la_test_solver PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
endif (OPENMP_FOUND)

if (CUDA_FOUND)
	cuda_add_executable(la_test_gpu_reduction test_gpu_reduction.cpp)
	target_link_libraries(la_test_gpu_reduction ${LINAL_CUDA_LIBS})
	cuda_add_cublas_to_target(la_test_gpu_reduction)
	set_target_properties(la_test_gpu_reduction PROPERTIES COMPILE_FLAGS "-DGPGPU")

	if (OPENMP_FOUND)
		set_target_properties(la_test_gpu_reduction PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
	endif (OPENMP_FOUND)

	add_test(la_test_gpu_reduction ${EXECUTABLE_OUTPUT_PATH}/la_test_gpu_reduction)
endif (CUDA_FOUND)

