from sage.modules.vector_rational_dense import Vector_rational_dense
from sage.matrix.constructor import matrix
from sage.rings.rational_field import QQ
from sage.modules.free_module_element import vector
from sage.all import *

# Input simple roots sets for G2, F4, E6, E7, E8

simplerootsG2dim3 = matrix(QQ[sqrt(3)], 
                           [[1,  -1], 
                            [-1,  2], 
                            [0 , -1]])

simplerootsG2dim2 = matrix(QQ[sqrt(3)], 
                           [[1, 0], 
                            [-3/2, sqrt(3)/2]])

simplerootsF4dim4 = matrix(QQ, 
                           [[0,1,-1,0], 
                            [0,0,1,-1], 
                            [0,0,0,1], 
                            [1/2, -1/2, -1/2, -1/2]])
print(simplerootsF4, '\n')

simplerootsE6dim6 = matrix(QQ[sqrt(3)], 
                                   [[1,-1,0,0,0,0], 
                                   [0,1,-1,0,0,0], 
                                   [0,0,1,-1,0,0], 
                                   [0,0,0,1,-1,0], 
                                   [0,0,0,1,1,0], 
                                   [-1/2, -1/2, -1/2, -1/2, -1/2, sqrt(3)/2]])

simplerootsE7dim7 = matrix(QQ[sqrt(2)], 
                                   [[1,-1,0,0,0,0,0], 
                                   [0,1,-1,0,0,0,0], 
                                   [0,0,1,-1,0,0,0], 
                                   [0,0,0,1,-1,0,0], 
                                   [0,0,0,0,1,-1,0], 
                                   [0,0,0,0,1,1,0], 
                                   [-1/2, -1/2, -1/2, -1/2, -1/2, -1/2, sqrt(2)/2]])

simplerootsE8dim8 = matrix(QQ, 
                          [[1,-1,0,0,0,0,0,0], 
                          [0,1,-1,0,0,0,0,0], 
                          [0,0,1,-1,0,0,0,0], 
                          [0,0,0,1,-1,0,0,0], 
                          [0,0,0,0,1,-1,0,0], 
                          [0,0,0,0,0,1,-1,0], 
                          [0,0,0,0,0,1,1,0], 
                          [-1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2]])

higherdimroots = [simplerootsG2dim3]
lowerdimroots = [simplerootsG2dim2]

for n in range(0,1):
  X = higherdimroots[n]
  Y = lowerdimroots[n]
  M = ( Y * X.transpose() ) * ( X * X.transpose() ).inverse()
  print(M)

