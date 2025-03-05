"""
Simple Groups

AUTHORS:

- Patricio Gallardo (2023): initial algorithms
- Jesus Martinez-Garcia (2023-25): initial algorithms and package version
- Han-Bon Moon (2023): initial algorithms
- David Swinarski (2023): initial algorithms
- Robert Hanson (2025): package version


REFERENCES:

- [GMGMS] P. Gallardo, J. Martinez-Garcia, H.-B. Moon, D. Swinarski. "Computation of GIT quotients of semisimple groups". Arxiv pre-print. arXiv:2308.08049
- [HMG] R. Hanson, J. Martinez-Garcia "CompGIT, a Sagemath package for Geometric Invariant Theory". To appear.
"""

from sage.modules.vector_rational_dense import Vector_rational_dense
#from sage.matrix.constructor import Matrix
from sage.rings.rational_field import QQ
from sage.modules.free_module_element import vector
from sage.all import *


#if type_A=True, the type is A and
#the vector has coordinates in basis H
class OneParamSubgroup(Vector_rational_dense):
    """
    Wrapper of Vector_rational_dense to store one-parameter subgroups by storing its weights (in H-coordinates).
    
    It should always be called using one_param_subgroup constructor.
    The reason for the wrapper is to make sure the class Vector_rational_dense is used.
    """
    def __init__(self, parent, value):
        Vector_rational_dense.__init__(self, parent, value)


def upper_triangular_entries(rnk, x, y):
    """
    Constructor for a square matrix with entries equal to 1 if they are in the diagonal or above the diagonal and 0 elsewhere.
    
    The function returns 1 for entry (x, y) if y>=x and 0 otherwise.
    An example of how to combine it with a lambda function is provided.
    
    EXAMPLES::
        
        sage: from SimpleGroup import upper_triangular_entries
        sage: upper_triangular_entries(2, 0, 1)
        1
        sage: upper_triangular_entries(2, 1, 0)
        0
        sage: print(matrix(QQ, 3, 3, lambda x,y: upper_triangular_entries(3,x,y)))
        [1 1 1]
        [0 1 1]
        [0 0 1]
    """
    return 1 if y >= x else 0


def lower_triangular_entries(rnk, x, y):
    """
    Constructor for a square matrix with entries equal to 1 if they are in the diagonal or below the diagonal and 0 elsewhere.
    
    The function returns 1 for entry (x,y) if y<=x and 0 otherwise.
    An example of how to combine it with a lambda function is provided.
    
    EXAMPLES::
        
        sage: from SimpleGroup import lower_triangular_entries
        sage: lower_triangular_entries(2, 0, 1)
        0
        sage: lower_triangular_entries(2, 1, 0)
        1
        sage: print(matrix(QQ, 3, 3, lambda x,y: lower_triangular_entries(3,x,y)))
        [1 0 0]
        [1 1 0]
        [1 1 1]

    """
    return upper_triangular_entries(rnk, y, x)

def one_param_subgroup(data, type_A=False, field=QQ):
    """
    Wrapper for OneParamSubgroup constructor that returns a OneParamSubgroup given in data as a vector.
    
    If type_A = True, it takes data in coordinates given by basis T and returns them in H-coordinates.
    Otherwise, it assumes the data is already in H-coordinates.
    
    EXAMPLES::
        
        sage: from SimpleGroup import one_param_subgroup
        sage: v1 = vector([2, 1, -3])
        sage: v2 = vector([3, 2, 1])
        sage: one_param_subgroup(v1, type_A=True)
        (2, 3)
        sage: one_param_subgroup(v2, type_A=False)
        (3, 2, 1)
    """
    if type_A:
        BasisChange=matrix(QQ, len(data)-1, len(data)-1,
                           lambda x,y: lower_triangular_entries(len(data)-1,x,y))
        v=vector(data[0:len(data)-1])
        v=tuple(BasisChange*v)
    
    else:
        v=tuple(data)
    # return OneParamSubgroup(QQ**len(v), v)
    return vector(field, v)

def A_coord_change_from_T_to_H(rnk,x, y):
    """
    Auxiliary function used to create a change of basis matrix from T coordinates to H coordinates.
    
    It returns what the entry ``x``, ``y`` for this matrix, for a group of rank ``rnk``.
    
    EXAMPLES::
        
        sage: from SimpleGroup import A_coord_change_from_T_to_H
        sage: x = 2
        sage: y = 3
        sage: rnk = 5
        sage: A_coord_change_from_T_to_H(rnk, x ,y)
        0 
    """
    if x==y:
        return 1
    elif x==y+1:
        return -1
    else:
        return 0

def inverse_of_upper_triangular(rnk, x, y):
    """
    Auxiliary function used to create the inverse of the matrix where the diagonal
    and the upper entries are ``1`` and all entries under the diagonal are ``0``.
    
    It returns the entry ``x``, ``y`` for this inverse, for a square matrix of size ``rnk``.
    
    EXAMPLES::
        
        sage: from SimpleGroup import inverse_of_upper_triangular
        sage: x = 2
        sage: y = 3
        sage: rnk = 5
        sage: inverse_of_upper_triangular(rnk, x, y)
        -1
    """
    if x==y:
        return 1
    elif y==x+1:
        return -1
    else:
        return 0

#This commented-out function is a constructor for the matrix of the rays of the
#fundamental chamber (in H basis) of a group of type A.
#It seems such expression in H coordinates is no longer needed, so we comment it out.
#def A_cone_basis_constructor(rnk, x, y):
#    """
#
#    EXAMPLES::
#
#        sage: from SimpleGroup import A_cone_basis_constructor
#        sage: x = 2
#        sage: y = 3
#        sage: rnk = 5
#        sage: A_cone_basis_constructor(rnk, x, y)
#        2
#    """
#    if x<=y:
#        return rnk-y
#    else:
#        return -y - 1


def A_cone_basis_constructor_from_T(rnk, x, y):
    """
    Auxiliary function used to list the rays of the fundamental chamber
    in T coordinates for a group of type A.
    
    It returns what the entry ``x``, ``y`` for this matrix, for a group of rank ``rnk``.
    
    EXAMPLES::
        
        sage: from SimpleGroup import A_cone_basis_constructor_from_T
        sage: x = 2
        sage: y = 3
        sage: rnk = 5
        sage: A_cone_basis_constructor_from_T(rnk, x, y)
        6
    """
    if x<=y:
        return (x+1) * (rnk-y)
    else:
        return (y+1) * (rnk-x)

def A_T_basis_constructor_from_gamma(rnk, x, y):
    """
    Auxiliary function used to created a change of coordinates matrix from T-coordinates
    to gamma-coordinates for a group of type A.
    
    It returns what the entry ``x``, ``y`` for this matrix, for a group of rank ``rnk``.


    EXAMPLES::
        
        sage: from SimpleGroup import A_T_basis_constructor_from_gamma
        sage: x = 2
        sage: y = 3
        sage: rnk = 5
        sage: A_T_basis_constructor_from_gamma(rnk, x, y)
        -1/6
    """
    if x == y:
        return 2/(rnk+1)
    if x == y + 1 or x == y - 1:
        return -1/(rnk+1)

def D_cone_basis_constructor(rnk, x, y):
    """
    Auxiliary function used to list the rays of the fundamental chamber
    in H coordinates for a group of type D.
    
    It returns what the entry ``x``, ``y`` for this matrix, for a group of rank ``rnk``.
    
    EXAMPLES::
        
        sage: from SimpleGroup import D_cone_basis_constructor
        sage: x = 2
        sage: y = 3
        sage: rnk = 5
        sage: D_cone_basis_constructor(rnk, x, y)
        1
    """
    if x<=y: # Index starting from 0?
        return 1
    elif x == rnk - 1 and y == rnk - 2:
        return -1
    else:
        return 0
        

def D_T_basis_constructor_from_gamma(rnk, x, y):
    """
    Auxiliary function used to created a change of coordinates matrix from T-coordinates
    to gamma-coordinates for a group of type D.
    
    It returns what the entry ``x``, ``y`` for this matrix, for a group of rank ``rnk``.

    EXAMPLES::
        
        sage: from SimpleGroup import D_T_basis_constructor_from_gamma
        sage: x = 2
        sage: y = 3
        sage: rnk = 5
        sage: D_T_basis_constructor_from_gamma(rnk, x, y)
        -1
    """
    if x == y:
        if x < rnk - 2:
            return 1
        else:
            return 1/2
    elif x == y-1:
        if x < rnk - 2:
            return -1
        else:
            return -1/2
    elif x == rnk - 1 and y == rnk - 2:
        return 1/2
    else:
        return 0




class SimpleGroup(object):
    """
    Wrapper of WeylGroup that includes additional data about a simple Lie group necessary to solve GIT problems.
    
    It returns an object representing data associated to a
    simple group of a certain Dynkin type and rank that can
    later be used by other methods. In particular, it
    encapsulates a fundamental domain (or fundamental
    (Weyl) chamber) in the Euclidean domain for the root
    system of the group.
    
    
    
    INPUT:
    
    - ``Dynkin_type`` -- The Dynkin type of the simple group, in string format. Currently only types ``"A", "B", "C", "D", "G"`` are implemented.
    - ``rnk`` -- rank of the group.
    
    
    
 
    
    INTERNAL ATTRIBUTES:
    
    Given a fixed maximal torus T in a reductive group G,
    there are a number of lattices and vector spaces that play
    a role studying the group. There are also a number of bases
    that can be useful when studying the group. An important basis
    is the one given by the rays of the fundamental chamber of the
    group (see [GMGMS] for details). Coordinates for this basis
    are referred to as gamma coordinates. We give names
    to the different bases coordinates. These are:
    
    - H-coordinates -- on the hom-spaces Hom(GG_m , T) of one parameter subgroups, 
    basis elements are given by the matrices H_i with only one non-zero element (i, i) of unitary size.  
    - L-coordinates are the dual coordinates to H, on the hom-spaces Hom(T, GG_m) of characters.
    - T-coordinates: basis elements are given by {T_i}_{i = 1, ..., n}, T_i = H_i - H_{i+1} in type A.
    - T-coordinates are equal to H-coordinates in type B, C, D.
    - Note the reduction from n+1 to n dimensions in type A (to account for the fact that weights
    of a one-parameter subgroup add to 0).
    - gamma-coordinates are given by gamma_i = H_1 + ... + H_i for type B, C, D, G.
    - gamma-coordinates can be obtained from a choice of simple roots of the root system in other
    cases. This is provided ad hoc for other groups (see [HMG] for details).
    
    With this in mind, the internal attributes used are:

    - ``Dynkin_type`` -- The Dynkin type of the simple group, in string format. Currently only types ``"A", "B", "C", "D", "G"`` are implemented.
    - ``max_torus_dim`` -- dimension of the maximal torus in the group. This is the same as rnk.
    - ``pairing_matrix`` -- The inner product matrix between characters and one-parameter subgroups. Usually the identity.
    - ``WeylGroup`` -- Object from class ``WeylGroup`` representing the group.
    - ``cone_basis`` -- Matrix containing the rays of a fundamental chamber/domain for in T coordinates.
    - ``T_to_gamma_change`` -- Change of basis matrix from T coordinates to coordinates in the basis given by the rays of the fundamental chamber of G.
    - ``T_to_H_change`` -- Change of basis matrix from T coordinates to H-coordinates.

    
    

    TODO::
    
    Types ``E6, E7, E8`` and ``F4``  should be implemented.
    
    EXAMPLES::
        
        sage: from SimpleGroup import SimpleGroup
        sage: G=SimpleGroup("A", 2)
        sage: G.group_type()
        'A'
        sage: G.rnk()
        2        
        sage: G.fundamental_chamber_generators() # rays of the fundamental chamber in T-coordinates
        [2 1]
        [1 2]
        sage: G.T_to_H_matrix()
        [ 1  0]
        [-1  1]
        [ 0 -1]        
                
        sage: H=SimpleGroup("B", 3)
        sage: H.group_type()
        'B'
        sage: H.rnk()
        3
        sage: H.fundamental_chamber_generators() # in T-coordinates 
        [1 1 1]
        [0 1 1]
        [0 0 1]
        sage: H.T_to_H_matrix() # T = H for type B, C, D
        [1 0 0]
        [0 1 0]
        [0 0 1]
        sage: H.T_to_gamma_matrix()
        [ 1 -1  0]
        [ 0  1 -1]
        [ 0  0  1]            
    """
    def __init__(self, Dynkin_type, rnk):
        self.Dynkin_type = Dynkin_type
        self.max_torus_dim = rnk
        self.pairing_matrix = matrix.identity(QQ, rnk)
        self.WeylGroup = WeylGroup([Dynkin_type, rnk])
        self.field = QQ

        if Dynkin_type == 'A':
#            self.cone_basis_in_H = matrix(QQ, rnk+1, rnk,
#                               lambda x,y: A_cone_basis_constructor(rnk, x, y)) #From gamma-coordinates to H-coordinates. Return rays of F in H-coordinates.
#            It is never used, so we comment it out instead of deleting it in case a future update finds it useful.
            self.cone_basis = matrix(QQ, rnk, rnk,
                                lambda x,y: A_cone_basis_constructor_from_T(rnk, x, y)) #stores rays of the fundamental chamber in T-coordinates. Change of coordinate matrix from gamma-coordinates to T-coordinates.
            self.T_to_gamma_change = matrix(QQ, rnk, rnk,
                                lambda x,y: A_T_basis_constructor_from_gamma(rnk, x, y)) #Change of basis matrix from T-coordinates to gamma-coordinates. Return T-vectors in gamma-coordinates.
            self.T_to_H_change = matrix(QQ, rnk+1, rnk,
                                lambda x,y: A_coord_change_from_T_to_H(rnk, x, y)) #Change of basis matrix from T-coordinates to H-coordinates. Returns T-vectors in H-coordinates.
            self.dual_basis = matrix.identity(QQ,rnk+1)
            for i in range(0,rnk-1):
                self.pairing_matrix[i,i+1]=-1
        
              
        elif Dynkin_type == 'B':
            self.cone_basis = matrix(QQ, rnk, rnk,
                                lambda x,y: upper_triangular_entries(rnk,x,y)) #stores rays of the fundamental chamber in H-coordinates. Change of coordinate matrix from gamma-coordinates to H-coordinates.
            self.T_to_gamma_change = matrix(QQ, rnk, rnk,
                                lambda x,y: inverse_of_upper_triangular(rnk,x,y)) #Change of coordinate matrix from T-coordinates to gamma-coordinates. Returns T-vectors in gamma-coordinates.
            self.T_to_H_change = matrix.identity(QQ, rnk) #Change of coordinate matrix from T-coordinates to H-coordinates. Returns T-vectors in H-coordinates.

        elif Dynkin_type == 'C':
            self.cone_basis = matrix(QQ, rnk, rnk,
                                lambda x,y: upper_triangular_entries(rnk,x,y)) #From gamma-coordinates to H-coordinates. Return rays of F in H-coordinates.
            self.T_to_H_change = matrix.identity(QQ, rnk) #Change of coordinate matrix from T-coordinates to H-coordinates. Returns T-vectors in H-coordinates.
            self.T_to_gamma_change = matrix(QQ, rnk, rnk,
                                lambda x,y: inverse_of_upper_triangular(rnk,x,y)) #Change of coordinate matrix from T-coordinates to gamma-coordinates. Returns T-vectors in gamma-coordinates.
            
        elif Dynkin_type == 'D':
            self.cone_basis = matrix(QQ, rnk, rnk,
                                lambda x,y: D_cone_basis_constructor(rnk,x,y)) #stores rays of the fundamental chamber in H-coordinates. Change of coordinate matrix from gamma-coordinates to H-coordinates.
            self.T_to_gamma_change = matrix(QQ, rnk, rnk,
                                lambda x,y: D_T_basis_constructor_from_gamma(rnk,x,y)) #Change of coordinate matrix from T-coordinates to gamma-coordinates. Returns T-vectors in gamma-coordinates.
            self.T_to_H_change = matrix.identity(QQ, rnk) #Change of coordinate matrix from T-coordinates to H-coordinates. Returns T-vectors in H-coordinates.
            
        elif Dynkin_type == 'E' and rnk == 6:
            self.cone_basis = matrix(QQ[sqrt(3)], [[3,   0,   0,   0,   0,    sqrt(3)],
                                                   [1/5, 1/5, 1/5, 1/5, 1/5,  sqrt(3)],
                                                   [1/3, 1/3, 1/3, 1/3, -1/3, sqrt(3)],
                                                   [1/2, 1/2, 0,   0,   0,    sqrt(3)],
                                                   [1/3, 1/3, 1/3, 0,   0,    sqrt(3)],
                                                   [0,   0,   0,   0,   0,    sqrt(3)]]
                                                   ).transpose()
            self.T_to_gamma_change = self.cone_basis.inverse()
            self.T_to_H_change = matrix.identity(QQ[sqrt(3)], 6)
            
        elif Dynkin_type == 'E' and rnk == 7:
            self.cone_basis = matrix(QQ[sqrt(2)], [[1/4, 1/4, 1/4, 1/4, 1/4, -1/4, sqrt(2)],
                                                   [1/6, 1/6, 1/6, 1/6, 1/6, 1/6,  sqrt(2)],
                                                   [0,   0,   0,   0,   0,   0,    sqrt(2)],
                                                   [2,   0,   0,   0,   0,   0,    sqrt(2)],
                                                   [1/2, 1/2, 0,   0,   0,   0,    sqrt(2)],
                                                   [1/4, 1/4, 1/4, 1/4, 0,   0,    sqrt(2)],
                                                   [1/3, 1/3, 1/3, 0,   0,   0,    sqrt(2)]]
                                                   ).transpose()
            self.T_to_gamma_change = self.cone_basis.inverse()
            self.T_to_H_change = matrix.identity(QQ[sqrt(2)], 7)
            
        elif Dynkin_type == 'E' and rnk == 8:
            self.cone_basis = matrix(QQ, [[-1,    0,    0,    0,    0,    0,    0,   1],
                                          [-1/3, -1/3, -1/3,  0,    0,    0,    0,   1],
                                          [-1/7, -1/7, -1/7, -1/7, -1/7, -1/7, -1/7, 1],
                                          [0,     0,    0,    0,    0,    0,    0,   1],
                                          [-1/5, -1/5, -1/5, -1/5, -1/5, -1/5, 1/5,  1],
                                          [-1/2, -1/2,  0,    0,    0,    0,    0,   1],
                                          [-1/5, -1/5, -1/5, -1/5, -1/5,  0,    0,   1],
                                          [-1/4, -1/4, -1/4, -1/4,  0,    0,    0,   1]]
                                          ).transpose()
            self.T_to_gamma_change = self.cone_basis.inverse()
            self.T_to_H_change = matrix.identity(QQ, 8)
            
        elif Dynkin_type == 'F':
            self.cone_basis = matrix(QQ, [[1, 1, 0, 0],
                                          [1, 0, 0, 0],
                                          [3, 1, 1, 1],
                                          [1, 1/2, 1/2, 0]]
                                          ).transpose()
            self.T_to_gamma_change = self.cone_basis.inverse()
            self.T_to_H_change = matrix.identity(QQ, 4)
            
        elif Dynkin_type == 'G':
            R=QuadraticField(3, 'a')
            self.field = R 
            M = Matrix(R, [[0,   1],[1, 3*sqrt(3)]])
            self.cone_basis = M.transpose() #stores rays of the fundamental chamber in H-coordinates. Change of coordinate matrix from gamma-coordinates to H-coordinates.
            self.T_to_gamma_change = self.cone_basis.inverse() # inverse matrix to self.cone_basis. Change of coordinate matrix from T-coordinates to gamma-coordinates. Returns T-vectors in gamma-coordinates.
            self.T_to_H_change = matrix.identity(R, 2)
        
        else:
            print ('Error: Dynkin type ', Dynkin_type, 'not supported/known')
            return None
    
    
    def Weyl_Group_elements(self):
        """
        Returns Dynkin type and rank within a string. 
        
            EXAMPLES::
            
            sage: from SimpleGroup import SimpleGroup 
            sage: G = SimpleGroup("A", 2)
            sage: G.Weyl_Group_elements()
            Weyl Group of type ['A', 2] (as a matrix group acting on the ambient space)
        """
        return self.WeylGroup

    
    def fundamental_chamber_generators(self):
        """
        Returns a matrix whose columns are the rays that generate the fundamental chamber of the Lie group. 
        
            EXAMPLES::
            
            sage: from SimpleGroup import SimpleGroup 
            sage: G = SimpleGroup("A", 2)
            sage: G.fundamental_chamber_generators()
            [2 1]
            [1 2]
        """
        return self.cone_basis

    
    def pairing(self, OPS, character_tuple):
        """
        
            EXAMPLES::
            
            sage: from SimpleGroup import SimpleGroup 
            sage: G = SimpleGroup("A", 2)
            sage: G.pairing(1, [1,2])
            (-1, 2)
        """
        character=vector(QQ,list(character_tuple))
        return (OPS * self.pairing_matrix) * character

    
    def fetch_pairing_matrix(self):
        """
        Returns the pairing on spaces of one parameter subgroups and characters as a matrix-valued bilinear form. 
        The matrix is defined in H and H-dual coordinates. 
        
            EXAMPLES::
            
            sage: from SimpleGroup import SimpleGroup 
            sage: G = SimpleGroup("A", 2)
            sage: G.fetch_pairing_matrix()
            [ 1 -1]
            [ 0  1]
        """
        return self.pairing_matrix

    
    def in_cone(self, OPS):
        """
        
            EXAMPLES::
            
            sage: from SimpleGroup import SimpleGroup 
            sage: G = SimpleGroup("A", 2)
            sage: G.in_cone(1)
            False
        """
        coordinates = (self.T_to_gamma_change) * OPS
        for i in range(self.max_torus_dim):
            if coordinates[i] < 0:
                return False
        return True

    
    def H_coordinates(self, OPS):
        """
        Returns H coordinates on the hom-spaces Hom(T, GG_m) of characters.
        
            EXAMPLES::
            
            sage: from SimpleGroup import SimpleGroup 
            sage: G = SimpleGroup("A", 2)
            sage: G.H_coordinates(1)
            [ 1  0]
            [-1  1]
            [ 0 -1]        
        """
        return self.T_to_H_change*OPS

    
    def group_type(self):
        """
        Returns the Dynkin type of the group as a string 'X', where X is A, B, C or D. 
        """
        return self.Dynkin_type

    
    def rnk(self):
        """
        Returns the rank of the group.
        """
        return self.max_torus_dim
        
    def T_to_gamma_matrix(self):
        """
        Returns basis T in the coordinates of gamma.
        """
        return self.T_to_gamma_change
    
    def T_to_H_matrix(self):
        """
        Returns basis T in the coordinates of H.
        """
        return self.T_to_H_change
    
    def lattice_field(self):
        """
        Returns the ground field used in the vector space containing the lattice.
        """
        return self.field
    
        
