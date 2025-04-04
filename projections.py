from sage.modules.vector_rational_dense import Vector_rational_dense
from sage.matrix.constructor import matrix
from sage.rings.rational_field import QQ
from sage.modules.free_module_element import vector
from sage.all import *

# Choices of simple roots sets for G2, F4, E6, E7, E8 are stored as matrices in which:
# - columns are choices of simple roots in a specified dimension
# - in higher dimensional cases (denoted X), orthogonal vectors are appended until the matrix is square
# - in lower dimensional cases (denoted Y), zero vectors are appended until the matrix is square


# G2
simplerootsG2dim3 = matrix(QQ[sqrt(3)],
                           [[ 1,  -1, 1],
                            [-1,   2, 1],
                            [ 0 , -1, 1]])

simplerootsG2dim2 = matrix(QQ[sqrt(3)],
                           [[1,    0,         0],
                            [-3/2, sqrt(3)/2, 0]]
                            ) # the final column is an orthogonal vector

X = simplerootsG2dim3
Y = simplerootsG2dim2
ProjG2 = ( Y * X.transpose() ) * (( X * X.transpose() ).inverse())


# E6
simplerootsE6dim9 = matrix([
                           [0,  0,  0,  0,  0,  1/3, 1, 0, 0],
                           [1,  0,  0,  0,  0, -2/3, 1, 0, 0],
                           [-1, 0,  0,  0,  0,  1/3, 1, 0, 0],
                           [0,  1,  0,  0,  0, -2/3, 0, 1, 0],
                           [0, -1,  1,  0,  0,  1/3, 0, 1, 0],
                           [0,  0, -1,  0,  0,  1/3, 0, 1, 0],
                           [0,  0,  0,  1,  0, -2/3, 0, 0, 1],
                           [0,  0,  0, -1,  1,  1/3, 0, 0, 1],
                           [0,  0,  0,  0, -1,  1/3, 0, 0, 1]
                           ]) # the final three columns are orthogonal vectors

simplerootsE6dim6 = matrix(QQ[sqrt(3)], [
                                   [1,    -1,    0,    0,    0,   0,          0, 0, 0],
                                   [0,     1,   -1,    0,    0,   0,          0, 0, 0],
                                   [0,     0,    1,   -1,    0,   0,          0, 0, 0],
                                   [0,     0,    0,    1,   -1,   0,          0, 0, 0],
                                   [0,     0,    0,    1,    1,   0,          0, 0, 0],
                                   [-1/2, -1/2, -1/2, -1/2, -1/2, sqrt(3)/2,  0, 0, 0]
                                   ])
  
X = simplerootsE6dim9
Y = simplerootsE6dim6
ProjE6 = ( Y * X.transpose() ) * (( X * X.transpose() ).inverse())


# E7
simplerootsE7dim8 = matrix([
                           [0,  0,  0,  0,   0,  0,  1/2, 1],
                           [-1, 0,  0,  0,   0,  0,  1/2, 1],
                           [1, -1,  0,  0,   0,  0,  1/2, 1],
                           [0,  1, -1,  0,   0,  0,  1/2, 1],
                           [0,  0,  1, -1,   0,  0, -1/2, 1],
                           [0,  0,  0,  1,  -1,  0, -1/2, 1],
                           [0,  0,  0,  0,   1, -1, -1/2, 1],
                           [0,  0,  0,  0,   0,  1, -1/2, 1]
                          ]) # the final column is an orthogonal vector

simplerootsE7dim7 = matrix(QQ[sqrt(2)],
                                   [[1,   -1,    0,    0,    0,    0,   0,          0],
                                   [ 0,    1,   -1,    0,    0,    0,   0,          0],
                                   [ 0,    0,    1,   -1,    0,    0,   0,          0],
                                   [ 0,    0,    0,    1,   -1,    0,   0,          0],
                                   [ 0,    0,    0,    0,    1,   -1,   0,          0],
                                   [ 0,    0,    0,    0,    1,    1,   0,          0],
                                   [-1/2, -1/2, -1/2, -1/2, -1/2, -1/2, sqrt(2)/2 , 0]
                                   ])

X = simplerootsE7dim8
Y = simplerootsE7dim7
ProjE7 = ( Y * X.transpose() ) * (( X * X.transpose() ).inverse())
