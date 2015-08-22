project(linal_test CXX C)

set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
include_directories(${LINAL_INCLUDE})

create_test_sourcelist(la_test_source
	la_test_driver.cpp
	test_solver.cpp
	test_sparse.cpp
	test_blas.cpp
	test_mm_dense.cpp
	test_lu.cpp
	)

add_executable(la_test ${la_test_source})
target_link_libraries(la_test ${LINAL_LIBS})

if (CUDA_FOUND)
	cuda_add_executable(la_test_gpu_reduction test_gpu_reduction.cpp)
	target_link_libraries(la_test_gpu_reduction ${LINAL_CUDA_LIBS})
	cuda_add_cublas_to_target(la_test_gpu_reduction)
	set_target_properties(la_test_gpu_reduction PROPERTIES COMPILE_FLAGS "-DGPGPU")

	cuda_add_executable(la_test_sparse_cu test_sparse.cpp)
	target_link_libraries(la_test_sparse_cu ${LINAL_CUDA_LIBS})
	cuda_add_cublas_to_target(la_test_sparse_cu)
	set_target_properties(la_test_sparse_cu PROPERTIES COMPILE_FLAGS "-DGPGPU")

	cuda_add_executable(la_test_mm_dense_cu test_mm_dense.cpp)
	target_link_libraries(la_test_mm_dense_cu ${LINAL_CUDA_LIBS})
	cuda_add_cublas_to_target(la_test_mm_dense_cu)
	set_target_properties(la_test_mm_dense_cu PROPERTIES COMPILE_FLAGS "-DGPGPU")

	add_test(la_test_gpu_reduction ${EXECUTABLE_OUTPUT_PATH}/la_test_gpu_reduction)
endif (CUDA_FOUND)

add_test(la_test_blas ${EXECUTABLE_OUTPUT_PATH}/la_test test_blas)

