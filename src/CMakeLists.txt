project(linal_src)

include_directories(.)

set(SOURCE
	linal_base.cpp
	linal.h
	gmres.cpp
	gmres.h

	impl/solver_impl.h
	impl/solver_umfpack.h
	impl/solver_superlu.h
	)

add_library(linal
	${SOURCE}
	linal.cpp
	)

if (CUDA_FOUND)
	cuda_add_library(linal_cu
		${SOURCE}
		gpu/linal_cuda.cpp
		gpu/linal_cuda.h
		gpu/alloc_cu.h
		gpu/linal.cu
		gpu/texture.h
		gpu/shmem.h
		gpu/reduction.h
		gpu/gmres.cu)
endif ()

