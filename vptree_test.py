import sys
import vptree
from scipy import spatial

import numpy as np


if __name__ == '__main__':

    X = np.random.randn(10, 10)
#    knn = spatial.cKDTree(X, leafsize=10)


    vptree = vptree.vptree(sparse = True)

    X = np.asarray([[[1,2], [2,5], [3,4]],
        [[1,5],[3,1],[5,1]],
        [[1,0.1],[2,5],[10,4]],
        [[2,4],[9,1],[10,1]],
        [[1,3], [3,3], [5,6]],
        [[1,3], [3,3], [5,7]],
        [[1,2], [10, 9], [100,10]],
        [[1,2], [100, 9], [101,10]],
        [[1,2], [100, 9], [121,10]],
        [[1,2], [100, 9], [131,10]],
        [[1,10], [100, 9], [131,10], [140,1]]
        ])

    vptree.init_sparse(X)

    print np.asarray(X).ndim

    for i in range(100):
        Y = np.array([[[1,10], [10, 9], [131,10], [150,12]], [[1,4]]])
        a = np.random.random_integers(100)
        Y[0][0][0] = np.random.random_integers(100)
        Y[0][0][1] = a
        a = np.random.random_integers(100)
        Y[0][1][1] = a
        a = np.random.random_integers(100)
        Y[0][2][1] = a
        a = np.random.random_integers(100)
        Y[1][0][1] = a

        X = np.concatenate((X, Y))
#    Y = np.array([[[1,10], [10, 9], [131,10], [150,12]], [[1,4]]])
#    X = np.concatenate((X, Y))

    val = np.asarray([[1,2], [2,5.], [3,4]])

    print 'search'
    for i in range(10000):
        res = vptree.search_sparse(val, 100, epsilon = 5.)

    print repr(res)
    print X[res[0][1]]

    sys.exit(0)

    for i in range(len(X)):
        res = vptree.search(X[i], 10)
        print repr(res)

    res = vptree.search(X[0][:5], 10)
