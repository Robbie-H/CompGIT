from sage.modules.vector_rational_dense import Vector_rational_dense
from sage.matrix.constructor import matrix
from sage.rings.rational_field import QQ
from sage.modules.free_module_element import vector
from sage.all import *

#since we are interested in the dual of the lattice where the root system lives, we will have the simple roots as rows and not columns.

# Input choices of simple roots sets for G2, F4, E6, E7, E8
simplerootsG2=matrix(QQ[sqrt(3)], [[1, 0], [-3/2, sqrt(3)/2]])
print(simplerootsG2, '\n')

simplerootsF4=matrix(QQ, [[0,1,-1,0], [0,0,1,-1], [0,0,0,1], [1/2, -1/2, -1/2, -1/2]])
print(simplerootsF4, '\n')

simplerootsE6=matrix(QQ[sqrt(3)], [[1,-1,0,0,0,0], [0,1,-1,0,0,0], [0,0,1,-1,0,0], [0,0,0,1,-1,0], [0,0,0,1,1,0], [-1/2, -1/2, -1/2, -1/2, -1/2, sqrt(3)/2]])
print(simplerootsE6, '\n')

simplerootsE7=matrix(QQ[sqrt(2)], [[1,-1,0,0,0,0,0], [0,1,-1,0,0,0,0], [0,0,1,-1,0,0,0], [0,0,0,1,-1,0,0], [0,0,0,0,1,-1,0], [0,0,0,0,1,1,0], [-1/2, -1/2, -1/2, -1/2, -1/2, -1/2, sqrt(2)/2]])
print(simplerootsE7, '\n')

simplerootsE8=matrix(QQ, [[1,-1,0,0,0,0,0,0], [0,1,-1,0,0,0,0,0], [0,0,1,-1,0,0,0,0], [0,0,0,1,-1,0,0,0], [0,0,0,0,1,-1,0,0], [0,0,0,0,0,1,-1,0], [0,0,0,0,0,1,1,0], [-1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2]])
print(simplerootsE8, '\n')



def find_rays(simpleroots):
    T=set()
    n=simpleroots.rank()
    print(n)
    all_index=tuple([0..n-1])
    for i in [0..n-1]:
        index=list(all_index[0:i]+all_index[i+1:n])
        linearsystem=simpleroots[index, all_index[0:n-1]]
        if (linearsystem.rank()==n-1):
            rhs=-1*simpleroots[index, [n-1]]
            solution=linearsystem.inverse()*rhs
            solution=solution.transpose()
            solution=solution.list()+[1]
        else:
            linearsystem=simpleroots[index, all_index[1:n]]
            rhs=-1*simpleroots[index, [0]]
            solution=linearsystem.inverse()*rhs
            solution=solution.transpose()
            solution=[1]+solution.list()
        solution=tuple(solution)
        T.add(solution)
    return T


#print gamma ray matrices for G2, F4, E6, E7, E8
raysG2=find_rays(simplerootsG2);
print("\nrays G2:")
print(raysG2, "\n")

raysF4=find_rays(simplerootsF4);
print("\nrays F4:")
print(raysF4, "\n")

raysE6=find_rays(simplerootsE6);
print("\nrays E6:")
print(raysE6, "\n")

raysE7=find_rays(simplerootsE7);
print("\nrays E7:")
print(raysE7, "\n")

raysE8=find_rays(simplerootsE8);
print("\nrays E8:")
print(raysE8, "\n")
