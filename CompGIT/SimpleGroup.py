from sage.modules.vector_rational_dense import Vector_rational_dense
from sage.matrix.constructor import matrix
from sage.rings.rational_field import QQ
from sage.modules.free_module_element import vector
from sage.all import *


#if type_A=True, the type is A and
#the vector has coordinates in basis H
class OneParamSubgroup(Vector_rational_dense):
    def __init__(self, parent, value):
        Vector_rational_dense.__init__(self, parent, value)


def upper_triangular_entries(dim, x, y):
    """
    
    
    EXAMPLES:
        
        sage: from SimpleGroup import upper_triangular_entries
        sage: upper_triangular_entries(2, [1, 2], [2,1])
        1
    """
    return 1 if y >= x else 0


def lower_triangular_entries(dim, x, y):
    """
    
    
    EXAMPLES:
        
        sage: from SimpleGroup import lower_triangular_entries
        sage: lower_triangular_entries(2, [1, 2], [2,1])
        0
    """
    return upper_triangular_entries(dim, y, x)


def one_param_subgroup(data, type_A=False):
    """
    For type_A=True, take vector in coordinates given by basis T.
    
    EXAMPLES::
        
        sage: from SimpleGroup import one_param_subgroup
        sage: v1 = vector([2,1,-3])
        sage: v2 = vector([3,2,1])
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
    return OneParamSubgroup(QQ**len(v), v)


def A_coord_change_from_T_to_H(dim,x,y):
    """
    
    Examples::
        
        sage: from SimpleGroup import A_coord_change_from_T_to_H
        sage: x=2
        sage: y=3
        sage: dim=5
        sage: A_coord_change_from_T_to_H(dim,x,y)
        0 
    """
    if x==y:
        return 1
    elif x==y+1:
        return -1
    else:
        return 0
    

def inverse_of_upper_triangular(dim,x,y):
    """
    
    Examples::
        
        sage: from SimpleGroup import inverse_of_upper_triangular
        sage: x=2
        sage: y=3
        sage: dim=5
        sage: inverse_of_upper_triangular(dim,x,y)
        -1
    """
    if x==y:
        return 1
    elif y==x+1:
        return -1
    else:
        return 0


def A_cone_basis_constructor(dim, x,y):
    """
    
    Examples::
        
        sage: from SimpleGroup import A_cone_basis_constructor
        sage: x=2
        sage: y=3
        sage: dim=5
        sage: A_cone_basis_constructor(dim,x,y)
        2
    """
    if x<=y:
        return dim-y
    else:
        return -y-1


def A_cone_basis_constructor_from_T(dim,x,y):
    """
    
    Examples::
        
        sage: from SimpleGroup import A_cone_basis_constructor_from_T
        sage: x=2
        sage: y=3
        sage: dim=5
        sage: A_cone_basis_constructor_from_T(dim,x,y)
        6
    """
    if x<=y:
        return (x+1)*(dim-y)
    else:
        return (y+1)*(dim-x)


def A_T_basis_constructor_from_gamma(dim,x,y):
    """
    
    Examples::
        
        sage: from SimpleGroup import A_T_basis_constructor_from_gamma
        sage: x=2
        sage: y=3
        sage: dim=5
        sage: A_T_basis_constructor_from_gamma(dim,x,y)
        -1/6
    """
    if x == y:
        return 2/(dim+1)
    if x == y+1 or x == y-1:
        return -1/(dim+1)


def D_cone_basis_constructor(dim, x,y):
    """
    
    Examples::
        
        sage: from SimpleGroup import D_cone_basis_constructor
        sage: x=2
        sage: y=3
        sage: dim=5
        sage: D_cone_basis_constructor(dim,x,y)
        1
    """
    if x<=y: # Index starting from 0?
        return 1
    elif x == dim-1 and y == dim-2:
        return -1
    else:
        return 0

    
def D_T_basis_constructor_from_gamma(dim,x,y):
    """
    
    Examples::
        
        sage: from SimpleGroup import D_T_basis_constructor_from_gamma
        sage: x=2
        sage: y=3
        sage: dim=5
        sage: D_T_basis_constructor_from_gamma(dim,x,y)
        -1
    """
    if x == y:
        if x < dim-2:
            return 1
        else:
            return 1/2
    elif x == y-1:
        if x < dim-2:
            return -1
        else:
            return -1/2
    elif x == dim-1 and y == dim-2:
        return 1/2
    else:
        return 0


def B_fundamental_weight_constructor(dim,x,y):
    """
    
    Examples::
        
        sage: from SimpleGroup import B_fundamental_weight_constructor
        sage: x=2
        sage: y=3
        sage: dim=5
        sage: B_fundamental_weight_constructor(dim,x,y)
        1
    """
    if y == dim-1:
        return 1/2
    elif x <= y:
        return 1
    else:
        return 0


def D_fundamental_weight_constructor(dim,x,y):
    """
    
    Examples::
        
        sage: from SimpleGroup import D_fundamental_weight_constructor
        sage: x=2
        sage: y=3
        sage: dim=5
        sage: D_fundamental_weight_constructor(dim,x,y)
        0.5
    """
    if x <= y:
        if y < dim-2:
            return 1
        else:
            return 1/2
    elif x == dim-1 and y == dim-2:
        return -1/2
    else:
        return 0


class SimpleGroup(object):
    """
    - Attributes include Dynkin type, rank, weights, characters and Weyl actions.
    - Dynkin types A, B, C and D are each treated seperately.
    - H-coordinates, on the hom-spaces Hom(GG_m , T) of one parameter subgroups, 
    are given by the matrices H_i with only one non-zero element (i, i) of unitary size.  
    - L-coordinates are the dual coordinates to H, on the hom-spaces Hom(T, GG_m) of characters.
    - T-coordinates are equal to H-coordinates in type B, C, D.
    - T-coordinates are given by {T_i}_{i = 1, ..., n}, T_i = H_i - H_{i+1} in in type A.
    - Note the reduction from n+1 to n dimensions in type A (to account for the fact that weights
    of a one-parameter subgroup add to 0).
    - gamma-coordinates are given by gamma_i = H_1 + ... + H_i for type B, C, D.
    
    EXAMPLES::
        
        sage: from SimpleGroup import SimpleGroup
        sage: G=SimpleGroup("A", 2)
        sage: G.Dynkin_type
        'A'
        sage: G.max_torus_dim
        2        
        sage: G.cone_basis # rays of the fundamental chamber (F) in T-coordinates
        [2 1]
        [1 2]
        sage: G.T_to_H_change  
        [ 1  0]
        [-1  1]
        [ 0 -1]        
        sage: G.cone_basis_in_H # rays of F in H-coordinates  
        [ 2  1]
        [-1  1]
        [-1 -2]
        sage: G.fundamental_weights # in L-coordinates
        [1 1]
        [0 1]
                
        sage: H=SimpleGroup("B", 3)
        sage: H.Dynkin_type
        'B'
        sage: H.max_torus_dim
        3
        sage: H.cone_basis # in T-coordinates 
        [1 1 1]
        [0 1 1]
        [0 0 1]
        sage: H.T_to_H_change # T = H for type B, C, D
        [1 0 0]
        [0 1 0]
        [0 0 1]
        sage: H.T_to_gamma_change    
        [ 1 -1  0]
        [ 0  1 -1]
        [ 0  0  1]            
    """
    def __init__(self, Dynkin_type, dim):
        self.Dynkin_type=Dynkin_type
        self.max_torus_dim=dim
        self.pairing_matrix=matrix.identity(QQ,dim)
        self.WeylGroup=WeylGroup([Dynkin_type, dim])

        if Dynkin_type=='A':
            self.lattice_standard_basis=matrix.identity(QQ,dim+1)
            self.cone_basis_in_H=matrix(QQ, dim+1, dim,
                                lambda x,y: A_cone_basis_constructor(dim,x,y)) #From gamma-coordinates to H-coordinates. Return rays of F in H-coordinates.
            self.cone_basis=matrix(QQ, dim, dim,
                                lambda x,y: A_cone_basis_constructor_from_T(dim,x,y)) # Return rays of F in T-coordinates.
            self.T_to_gamma_change=matrix(QQ, dim, dim,
                                lambda x,y: A_T_basis_constructor_from_gamma(dim,x,y)) #From T-coordinates to gamma-coordinates. Return T-vectors in gamma-coordinates.
            self.T_to_H_change=matrix(QQ, dim+1, dim,
                                lambda x,y: A_coord_change_from_T_to_H(dim,x,y)) #From T-coordinates to H-coordinates. Return T-vectors in gamma-coordinates.            self.fundamental_weights=matrix(QQ,
            self.dual_basis=matrix.identity(QQ,dim+1)
            self.fundamental_weights=matrix(QQ, dim, dim, lambda x,y: upper_triangular_entries(dim,x,y)) #Return fundamental weights in L-coordinates.
            for i in range(0,dim-1):
                self.pairing_matrix[i,i+1]=-1
        
              
        elif Dynkin_type=='B':
            self.lattice_standard_basis=matrix.identity(QQ,dim)
            self.cone_basis=matrix(QQ, dim, dim,
                                lambda x,y: upper_triangular_entries(dim,x,y)) #From gamma-coordinates to H-coordinates. Return rays of F in H-coordinates.
            self.T_to_gamma_change=matrix(QQ, dim, dim,
                                lambda x,y: inverse_of_upper_triangular(dim,x,y)) #From T-coordinates to gamma-coordinates. Return T-vectors in H-coordinates.
            self.T_to_H_change=matrix.identity(QQ, dim) #From T-coordinates to H-coordinates. Return T-vectors in gamma-coordinates.            self.fundamental_weights=matrix(QQ, dim, dim, lambda x,y: B_fundamental_weight_constructor(dim,x,y)) #Return fundamental weights in L-coordinates.

        elif Dynkin_type=='C':
            self.lattice_standard_basis=matrix.identity(QQ,dim)
            self.cone_basis=matrix(QQ, dim, dim,
                                lambda x,y: upper_triangular_entries(dim,x,y)) #From gamma-coordinates to H-coordinates. Return rays of F in H-coordinates.
            self.T_to_H_change=matrix.identity(QQ, dim) #From T-coordinates to H-coordinates. Return T-vectors in H-coordinates.
            #self.fundamental_weights=matrix(QQ,
            self.T_to_gamma_change=matrix(QQ, dim, dim,
                                lambda x,y: inverse_of_upper_triangular(dim,x,y)) #From T-coordinates to gamma-coordinates. Return T-vectors in gamma-coordinates.
            self.fundamental_weights=matrix(QQ, dim, dim, lambda x,y: upper_triangular_entries(dim,x,y)) #Return fundamental weights in L-coordinates.
            
        elif Dynkin_type=='D':
            self.lattice_standard_basis=matrix.identity(QQ,dim)
            self.cone_basis=matrix(QQ, dim, dim,
                                lambda x,y: D_cone_basis_constructor(dim,x,y)) #From gamma-coordinates to H-coordinates. Return rays of F in H-coordinates.
            self.T_to_gamma_change=matrix(QQ, dim, dim,
                                lambda x,y: D_T_basis_constructor_from_gamma(dim,x,y)) #Return T-vectors in gamma-coordinates.
            self.T_to_H_change=matrix.identity(QQ, dim) #From T-coordinates to H-coordinates. Return T-vectors in H-coordinates.            self.fundamental_weights=matrix(QQ,
            self.fundamental_weights=matrix(QQ, dim, dim, lambda x,y: D_fundamental_weight_constructor(dim,x,y)) #Return fundamental weights in L-coordinates.
        else:
            print ('Error: Dynkin type ', Dynkin_type, 'not supported/known')
            return None

    
    def Weyl_Group_elements(self):
        """
        
            EXAMPLES::
            
            sage: from SimpleGroup import SimpleGroup 
            sage: G=SimpleGroup("A", 2)
            sage: G.Weyl_Group_elements()
            Weyl Group of type ['A', 2] (as a matrix group acting on the ambient space)
        """
        return self.WeylGroup

    
    def fundamental_chamber_generators(self):
        """
        
            EXAMPLES::
            
            sage: from SimpleGroup import SimpleGroup 
            sage: G=SimpleGroup("A", 2)
            sage: G.fundamental_chamber_generators()
            [2 1]
            [1 2]
        """
        return self.cone_basis

    
    def pairing(self, OPS, character_tuple):
        """
        
            EXAMPLES::
            
            sage: from SimpleGroup import SimpleGroup 
            sage: G=SimpleGroup("A", 2)
            sage: G.pairing(1,[1,2])
            (-1, 2)
        """
        character=vector(QQ,list(character_tuple))
        return (OPS*self.pairing_matrix)*character

    
    def fetch_pairing_matrix(self):
        """
        
            EXAMPLES::
            
            sage: from SimpleGroup import SimpleGroup 
            sage: G=SimpleGroup("A", 2)
            sage: G.fetch_pairing_matrix()
            [ 1 -1]
            [ 0  1]
        """
        return self.pairing_matrix

    
    def in_cone(self, OPS):
        """
        
            EXAMPLES::
            
            sage: from SimpleGroup import SimpleGroup 
            sage: G=SimpleGroup("A", 2)
            sage: G.in_cone(1)
            False
        """
        coordinates=(self.T_to_gamma_change)*OPS
        for i in range(self.max_torus_dim):
            if coordinates[i]<0:
                return False
        return True

    
    def H_coordinates(self, OPS):
        """
        
            EXAMPLES::
            
            sage: from SimpleGroup import SimpleGroup 
            sage: G=SimpleGroup("A", 2)
            sage: G.H_coordinates(1)
            [ 1  0]
            [-1  1]
            [ 0 -1]        
        """
        return self.T_to_H_change*OPS

    
    def Dynkin_type(self):
        return self.Dynkin_type

    
        
