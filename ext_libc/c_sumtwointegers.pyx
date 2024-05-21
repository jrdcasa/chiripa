
# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# declare the interface to the C code
cdef extern from "calc_sum_two_integers.c":
    cdef bint USED_OPENMP
    int c_sumtwointegers (double *ref, double *other, int n, int m )

OPENMP_ENABLED = True if USED_OPENMP else False

def calc_sum_two_integers(np.ndarray[np.float64_t, ndim=2, mode="c"] ref,
                          np.ndarray[np.float64_t, ndim=2, mode="c"] other,
                          int s1,
                          int s2
                          ):

    cdef int rows = ref.shape[0]
    cdef int cols = other.shape[0]

    r = c_sumtwointegers(&ref[0,0], &other[0,0], s1, s2 )

    return r
