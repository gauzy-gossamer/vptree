# distutils: language = c++

from libcpp cimport bool
from libcpp.vector cimport vector
cimport numpy as np
import numpy as np

cdef extern from "float.h":
    double DBL_MAX

cdef extern from "limits.h":
    int INT_MAX

cdef extern from "vptree.h":
    cdef cppclass DataPoint
    cdef cppclass DataPointSparse:
        DataPointSparse(int D, int ind, double *x)
    cdef cppclass VpTree[T]:
        VpTree()
        void set_sparse()
        void create(vector[T] *items)
        void create(double *vals, int N, int D)
        void search(double *val, int k, int D, double epsilon, vector[int] *indices, vector[double] *distances)
        void view(double *vals, int N)

cdef class vptree(object):
    cdef VpTree[DataPoint]* tree
    cdef VpTree[DataPointSparse]* tree_sparse
    cdef double[:, ::1] Y
    cdef bool sparse
    def __init__(self, sparse = False):
        self.sparse = sparse
        if sparse:
            self.tree_sparse = new VpTree[DataPointSparse]()
            self.tree_sparse.set_sparse()
        else:
            self.tree = new VpTree[DataPoint]()

    def init(self, double[:, ::1] X):
        cdef int N = X.shape[0]
        cdef int D = X.shape[1]
        if self.sparse:
            print 'failed'
        else:
            self.tree.create(&X[0,0], D, N)

    def init_sparse(self, np.ndarray X):
        cdef double [:,:] D
        cdef vector[DataPointSparse] *DP
        if not self.sparse:
            print 'failed'
        else:
            DP = new vector[DataPointSparse]()

            for i in range(X.shape[0]):
                D = np.ascontiguousarray(X[i], dtype=np.double)
                DP.push_back(DataPointSparse(D.shape[0], i, &D[0,0]))
            self.tree_sparse.create(DP)

    def search(self, double[:] V, int k = INT_MAX, double epsilon = DBL_MAX):
        cdef vector[int] res
        cdef vector[double] dists

        if k == INT_MAX and epsilon == DBL_MAX:
            raise ValueError("either k or epsilon should be passed")

        self.tree.search(&V[0], k, V.shape[0], epsilon, &res, &dists)
        return (res, dists)

    def search_sparse(self, double[:,::1] V, int k = INT_MAX, double epsilon = DBL_MAX):
        cdef vector[int] res
        cdef vector[double] dists

        if k == INT_MAX and epsilon == DBL_MAX:
            raise ValueError("either k or epsilon should be passed")

        self.tree_sparse.search(&V[0,0], k, V.shape[0], epsilon, &res, &dists)
        return (res, dists)
