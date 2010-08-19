project(linal_src)

include_directories(.)

set(SOURCE
	linal_base.cpp
	linal.h

	pow.h
	timer.cpp
	timer.h
	fpe.cpp
	fpe.h
	array.h

	gmres.h

	impl/solver_impl.h
	impl/solver_umfpack.h
	impl/solver_superlu.h
	)

add_library(linal
	${SOURCE}
	linal.cpp
	lu_solve.cpp
	lu_solve.h
	mv_dense.cpp
	mv_dense.h
	mv_csr.cpp
	mv_csr.h
	mv_ell.cpp
	mv_ell.h
	)

if (CUDA_FOUND)
	cuda_add_library(linal_cu
		${SOURCE}
		gpu/linal_cuda.cpp
		gpu/linal_cuda.h
		gpu/alloc_cu.h
		gpu/linal.cu
		gpu/mv_csr.cu
		gpu/mv_ell.cu
		gpu/texture.h
		gpu/shmem.h
		gpu/reduction.h
		gpu/gmres.cu
		OPTIONS --generate-code arch=compute_11,code=compute_11 --generate-code arch=compute_13,code=compute_13 --generate-code arch=compute_20,code=compute_20
		)
endif ()

