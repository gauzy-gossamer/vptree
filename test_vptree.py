from scipy import spatial
import numpy as np

import vptree

if __name__ == '__main__':
    X = np.random.randn(100, 2)

    Y = []
    for i in range(len(X)):
        Y.append([[0,X[i][0]],[1,X[i][1]]])

    Y = np.array(Y)

#    print repr(X)
#    print repr(Y)

    knn = spatial.cKDTree(X, leafsize=100)

    vp = vptree.vptree()
    vp_sparse = vptree.vptree(sparse = True)
    vp.init(X)
    vp_sparse.init_sparse(Y)

    j = 100
    k = 10

    for i in range(len(X)):
        res = knn.query(X[i], k)
        res1 = vp.search(X[i], k)
        res2 = vp_sparse.search_sparse(Y[i], k)

        #print repr(res)
        #print repr(res1)
        #print repr(res2)

        assert(res[1].tolist() == res1[0])
        assert(res1[0] == res2[0])

    print 'test OK'
