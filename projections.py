from sage.modules.vector_rational_dense import Vector_rational_dense
from sage.matrix.constructor import matrix
from sage.rings.rational_field import QQ
from sage.modules.free_module_element import vector
from sage.all import *

# Choices of simple roots sets for G2, F4, E6, E7, E8 are stored as matrices in which:
# - columns are choices of simple roots in a specified dimension
# - in higher dimensional cases, an orthogonal vector is appended in the last column
# - in lower dimensional cases, a zero vector is appended in the last column

simplerootsG2dim3 = matrix(QQ[sqrt(3)], 
                           [[1,  -1, 1],
                            [-1,  2, 1],
                            [0 , -1, 1]])

simplerootsG2dim2 = matrix(QQ[sqrt(3)], 
                           [[1, 0, 0],
                            [-3/2, sqrt(3)/2, 0]])


# To-do: input simplerootsE6dim8

simplerootsE6dim6 = matrix(QQ[sqrt(3)],
                                   [[1,-1,0,0,0,0], 
                                   [0,1,-1,0,0,0], 
                                   [0,0,1,-1,0,0], 
                                   [0,0,0,1,-1,0], 
                                   [0,0,0,1,1,0], 
                                   [-1/2, -1/2, -1/2, -1/2, -1/2, sqrt(3)/2]])

# To-do: input simplerootsE7dim8

simplerootsE7dim7 = matrix(QQ[sqrt(2)],
                                   [[1,-1,0,0,0,0,0], 
                                   [0,1,-1,0,0,0,0], 
                                   [0,0,1,-1,0,0,0], 
                                   [0,0,0,1,-1,0,0], 
                                   [0,0,0,0,1,-1,0], 
                                   [0,0,0,0,1,1,0], 
                                   [-1/2, -1/2, -1/2, -1/2, -1/2, -1/2, sqrt(2)/2]])

higherdimroots = [simplerootsG2dim3]  # X's
lowerdimroots =  [simplerootsG2dim2]  # Y's

# solve M * X = Y for M
# currently implemented for G2

for n in range(0,1):
    X = higherdimroots[n]
    Y = lowerdimroots[n]
    M = ( Y * X.transpose() ) * (( X * X.transpose() ).inverse())
    print(M)
    
