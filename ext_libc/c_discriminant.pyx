import cython
# # import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free

# declare the interface to the C code
cdef extern from "calc_discriminant.c":
    cdef bint USED_OPENMP
    void c_discriminant (double ref[], double conf[], double n[],
                         int Ni, int Nj, int dim,
                         float ptr_Rvdw_i[], float ptr_Rvdw_j[],
                         double rmaxmin[2], int max_index[2], int min_index[2])

    void c_discriminant_okuwaki (double ref[], double conf[], double n[],
                         int Ni, int Nj, int dim,
                         double dij_base[],
                         double rmaxmin[2], int max_index[2], int min_index[2])

def calc_discriminant(np.ndarray[np.float64_t, ndim=2, mode="c"] ref,
                      np.ndarray[np.float64_t, ndim=2, mode="c"] conf,
                      np.ndarray[np.float64_t, ndim=1, mode="c"] n,
                      int Ni,
                      int Nj,
                      Rvdw_i,
                      Rvdw_j):

    ptr_Rvdw_i = <float *>malloc(Ni*cython.sizeof(float))
    ptr_Rvdw_j = <float *>malloc(Nj*cython.sizeof(float))

    cdef double rmaxmin[2]
    cdef int max_index[2]
    cdef int min_index[2]

    if ptr_Rvdw_i is NULL:
            raise MemoryError()
    for i in range(Ni):
            ptr_Rvdw_i[i] = Rvdw_i[i]

    if ptr_Rvdw_j is NULL:
            raise MemoryError()
    for j in range(Nj):
            ptr_Rvdw_j[j] = Rvdw_j[j]

    cdef int dim1 = ref.shape[1]
    cdef int dim2 = conf.shape[1]
    if dim1 != dim2:
        print ("ERROR. The dimension of points to calculate distances must be equal!!!!")

    c_discriminant(&ref[0,0], &conf[0,0], &n[0], Ni, Nj,
                   dim1, ptr_Rvdw_i, ptr_Rvdw_j, rmaxmin, max_index, min_index)


    free(ptr_Rvdw_i)
    free(ptr_Rvdw_j)

    return rmaxmin[0], rmaxmin[1], [max_index[0], max_index[1]], [min_index[0], min_index[1]]

def calc_discriminant_okuwaki(np.ndarray[np.float64_t, ndim=2, mode="c"] ref,
                      np.ndarray[np.float64_t, ndim=2, mode="c"] conf,
                      np.ndarray[np.float64_t, ndim=1, mode="c"] n,
                      int Ni,
                      int Nj,
                      np.ndarray[np.float64_t, ndim=2, mode="c"] dij_base):

    cdef double rmaxmin[2]
    cdef int max_index[2]
    cdef int min_index[2]

    cdef int dim1 = ref.shape[1]
    cdef int dim2 = conf.shape[1]
    if dim1 != dim2:
        print ("ERROR. The dimension of points to calculate distances must be equal!!!!")

    c_discriminant_okuwaki(&ref[0,0], &conf[0,0], &n[0], Ni, Nj,
                    dim1, &dij_base[0,0], rmaxmin, max_index, min_index)


    return rmaxmin[0], rmaxmin[1], [max_index[0], max_index[1]], [min_index[0], min_index[1]]


