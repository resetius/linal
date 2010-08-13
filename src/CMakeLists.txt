project(linal_src)

include_directories(.)

set(SOURCE
	linal_base.cpp
	linal.h

	linal_util.h
	util.cpp

	lu_solve.c

	gmres.h

	impl/solver_impl.h
	impl/solver_umfpack.h
	impl/solver_superlu.h
	)

add_library(linal
	${SOURCE}
	linal.cpp
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
		gpu/gmres.cu)
endif ()

