debug=False;
"""
GIT (Geometric Invariant Theory) package

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
from copy import copy
from operator import itemgetter

from sage.combinat.subset import Subsets
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.matrix.constructor import Matrix
from sage.modules.free_module import VectorSpace
from sage.rings.rational_field import QQ
from sage.sets.set import Set
from sage.all import *
from SimpleGroup import *


def proportional(v1, v2):
    """
    Decides if two non-zero vectors are proportional.
    
    The vectors inputed can be of any class that allows iteration and comparison of individual entries. They must be of the same size. Returns ``True`` if they are proportional and ``False`` otherwise.
    
    INPUT:
    - ``v1`` -- Vector v1
    - ``v2`` -- Vector v2
    
    EXAMPLES::
        
        sage: from Git import proportional
        sage: v1 = vector([2,0,6])
        sage: v2 = vector([6,9,1])
        sage: proportional(v1, 2*v1)
        True 
        sage: proportional(v1, v2)
        False
    """
    i = next(i for i, ci in enumerate(v1) if ci != 0)
    m1 = v1[i]
    m2 = v2[i]
    return all(m2 * c1 == m1 * c2 for c1, c2 in zip(v1, v2))


def weights_matrix(weights_set):
    """
    Takes a ``Set`` of weights (vectors) and returns it as a matrix of row vectors (type ``list`` over ``QQ``.
    
    INPUT:
    - ``weights_set`` -- The Set of weights.
    
    EXAMPLES::
                
        sage: from Git import weights_matrix
        sage: weights = ([1,2], [3,4])
        sage: weights_matrix(weights)
        [1 2]
        [3 4]
    """
    return Matrix(QQ, [list(weight) for weight in weights_set])

# computes length vector and picks up minimum entry.


# We will compute the trivial character of the representation using
# the function below, in case we are ever working in coordinates
# for which chi_0 is not the origin


def averageWeight(x):
    """
    Takes a tuple of ``m`` vectors (each of them expressed as a tuple of numbers) ``v^i=(v^i_1, ... v^i_n)``, ``i=1..m`` and returns the average vector ``v``, with j^th entry ``v_j=( v^1_j + ... + v^m_j ) / m``
    
    INPUT:
    ``x`` -- A tuple of vectors, each of them expressed as a tuple of numbers.
    
    EXAMPLES::
        
        sage: from Git import averageWeight
        sage: weights_set = ([1,2], [3,4])
        sage: averageWeight(weights_set)
        (2, 3)
    """
    n=len(x[0])
    N=len(x)
    xbar=[0 for j in range(n)]
    for i in range(N):
        for j in range(n):
            xbar[j] = xbar[j]+x[i][j]
    for j in range(n):
        xbar[j] = xbar[j] / N
    return tuple(xbar)






class GITProblem(object):
    """
    Class to solve GIT problems consisting on a simple connected group acting on projective space.
    
    This class that encapsulates a representation of the action of a simple group, with
    designated weights. Its methods can find the maximal families of non-stable, unstable and strictly polystable loci with respect to a fixed torus, storing these families inside the object.  When an object of the class GITProblem is created, a number of sets of weights (determining any family) is created and stored in the object, to later be used by the class's methods to find the maximal families (determined by maximal sets of weights) of non-stable, unstable and strictly polystable points. When creating the object, an object of the class ``SimpleGroup`` is created. 
    
    INPUT::
    - ``rep`` -- A representation (an object of the class ``sage.combinat.root_system.weyl_characters.WeylCharacterRing_with_category.element_class``, see Examples to see an easy way of creating such objects).
    - ``label`` -- (optional) a string to label the GIT problem for easy tracking, e.g. ``label='plane cubics'``.
    
    
    EXAMPLES::
                 
        # Cubics in P2
        sage: from Git import GITProblem 
        sage: Phi = WeylCharacterRing("A2")
        sage: representation= Phi(3,0,0)
        sage: P=GITProblem(representation,label="Plane cubics")
        sage: P.solve_non_stable(Weyl_optimisation=True)
        {{(1, 2), (2, 1), (0, 0), (-1, 1), (0, 3), (1, -1), (3, 0)}, {(1, 2), (-1, -2), (2, 1), (0, 0), (1, -1), (3, 0)}}
        sage: P.print_solution_nonstable()
        <BLANKLINE>                
        <BLANKLINE>                
        ***************************************
        SOLUTION TO GIT PROBLEM: NONSTABLE LOCI
        ***************************************
        Group: A2
        Representation  A2(3,0,0)
        Set of maximal non-stable states:
        (1) 1-PS = (1, 1, -2) yields a state with 7 characters
        Maximal nonstable state={ (1, 2, 0), (2, 1, 0), (1, 1, 1), (0, 2, 1), (0, 3, 0), (2, 0, 1), (3, 0, 0) }
        (2) 1-PS = (1, -1/2, -1/2) yields a state with 6 characters
        Maximal nonstable state={ (1, 2, 0), (1, 0, 2), (2, 1, 0), (1, 1, 1), (2, 0, 1), (3, 0, 0) }
                  
        sage: P.solve_unstable(Weyl_optimisation=True)
        {{(1, 2), (2, 1), (0, 3), (1, -1), (3, 0)}}
        sage: P.print_solution_unstable()
        <BLANKLINE>       
        <BLANKLINE>                
        **************************************
        SOLUTION TO GIT PROBLEM: UNSTABLE LOCI
        **************************************
        Group: A2
        Representation  A2(3,0,0)
        Set of maximal unstable states:
        (1) 1-PS = (1, 1/4, -5/4) yields a state with 5 characters
        Maximal unstable state={ (1, 2, 0), (2, 1, 0), (0, 3, 0), (2, 0, 1), (3, 0, 0) }
           
        sage: P.solve_strictly_polystable()
        {{(0, 0)}, {(-1, 1), (1, -1), (0, 0)}}
        sage: P.print_solution_strictly_polystable()
        <BLANKLINE>             
        <BLANKLINE>        
        *************************************************************
        SOLUTION TO GIT PROBLEM: STRICTLY POLYSTABLE LOCI
        *************************************************************
        Group: A2
        Representation  A2(3,0,0)
        Set of strictly T-polystable states:
        (1) A state with 1 characters
        Strictly polystable state={ (1, 1, 1) }
        (2) A state with 3 characters
        Strictly polystable state={ (0, 2, 1), (2, 0, 1), (1, 1, 1) }
        """
    def __init__(self, rep,label=''):
        pair=rep.cartan_type()
        self.label=label
        self.Dynkin_type=pair[0]
        self.rank=pair[1]
        self.rep=rep
        self.group=SimpleGroup(self.Dynkin_type, self.rank)
        if self.group is None:
            print ('Group {Dynkin_type}{rank} not yet implemented'.format(Dynkin_type=self.Dynkin_type, rank=self.rank), sep='')
            return None

        #Transforming weights into tuple form
        weights_dict=self.rep.weight_multiplicities()
        weights=tuple(weights_dict.keys())
        if self.Dynkin_type=='A':
            length = self.rank + 1
        else:
            length = self.rank
        weights=tuple([tuple([weight[i] for i in range(length)]) for weight in weights])
        if self.Dynkin_type=='A':
            H_weights=copy(weights)
            self.weights=tuple([tuple([weight[i]-weight[len(weight)-1] for i in range(len(weight)-1)]) for weight in weights])
            conversion_dictionary={}
            for i in range(len(weights)):
                conversion_dictionary[self.weights[i]]=H_weights[i]
            self.L_coord_to_H_dual_conversion=conversion_dictionary
        else:
            self.weights=weights
        self.trivial_character=averageWeight(self.weights)
        self.optimized_weights_non_stable=None
        self.maximal_nonstable_states=None
        self.nonstable_weights_candidates=None
        self.weights_in_all_unstable_states=None
        self.maximal_semistable_states=None
        self.strictly_polystable_states=None
        self.gamma_OPS_nonstable_dictionary={}
        self.gamma_OPS_unstable_dictionary={}
        self.gamma_OPS_strictly_polystable_dictionary={}
        # self.all_states=Set([z for z in powerset(list(self.weights))])   #NOT REALLY NEEDED, HEAVY TO FORM COMPUTATIONALLY
        self.fundamental_chamber_generators=self.group.fundamental_chamber_generators()
        
        #Compute the states destabilised by each generator of the Weyl fundamental chamber
        self.OPS_rays_list=([one_param_subgroup(tuple(self.fundamental_chamber_generators[:,i].transpose())[0], field=self.group.lattice_field()) for i in range(0,self.rank)])
        states_destabilized_by_rays=[Set(self.destabilized_weights(OPS, all_weights_considered=True)) for OPS in self.OPS_rays_list]
        self.states_destabilized_by_rays_strict=[Set(self.destabilized_weights(OPS, all_weights_considered=True, strict_inequality=True)) for OPS in self.OPS_rays_list]
        
            
        #Compute the weights that may be in SOME destabilised state (union) both for later computing the non-stable states and the unstable states
        self.nonstable_weights_destabilized=Set([])
        self.unstable_weights_destabilized=Set([])
        for i in range(self.rank):
            self.nonstable_weights_destabilized=self.nonstable_weights_destabilized.union(states_destabilized_by_rays[i])
            self.unstable_weights_destabilized=self.unstable_weights_destabilized.union(self.states_destabilized_by_rays_strict[i])
        
        #Compute the weights that must be in ALL nonstable states and a superset of the weights that must be in all unstable states (intersection)
        self.weights_in_all_nonstable_states=Set(self.weights)
        for i in range(self.rank):
            self.weights_in_all_nonstable_states=self.weights_in_all_nonstable_states.intersection(self.states_destabilized_by_rays_strict[i])
        
        #Compute the weights that must be in ALL unstable states
        self.compute_weights_in_all_unstable_states()
            
        #Compute the difference of both sets (these are the weights
        #that must be tested) (difference between union and intersection)
        self.nonstable_weights_candidates=self.nonstable_weights_destabilized.difference(self.weights_in_all_nonstable_states)
        self.unstable_weights_candidates=self.unstable_weights_destabilized.difference(self.weights_in_all_unstable_states)


    def Weyl_group(self):
        """
        Returns a set of matrices of the Weyl group of the simple group
        associated to the object in the class ``GITProblem``.
        
        EXAMPLES::
                
            sage: from Git import GITProblem 
            sage: Phi = WeylCharacterRing("A2")
            sage: representation= Phi(3,0,0)
            sage: P=GITProblem(representation,label="Plane cubics")
            sage: print(GITProblem.Weyl_group(P))
            Weyl Group of type ['A', 2] (as a matrix group acting on the ambient space)
        """
        if self.Dynkin_type=='G' or (self.Dynkin_type=='E' and (self.rank==6 or self.rank==7)):
            group_temp=self.group.Weyl_Group_elements();
            return [group_temp.reflection_representation().representation_matrix(g) for g in group_temp]
        return self.group.Weyl_Group_elements()

    def weyl_elt_action_on_state(self,M,state):
        """
        It returns the ``Set`` of weights corresponding to acting on the ``Set`` ``state`` via the element ``M`` in the Weyl group of the problem.
        
        INPUT:
        
        - ``M`` -- a group element in the Weyl group of the group in the problem.
        - ``state`` -- a specific state.
        """
        L=list(state)
        if self.Dynkin_type=='A':
            L=[self.L_coord_to_H_dual_conversion[x] for x in L]
        ML=[M*((Matrix(v)).transpose()) for v in L]
        L=[[x[0] for x in y] for y in ML]
        if self.Dynkin_type=='A':
            L=[ [weight[i]-weight[len(weight)-1] for i in range(len(weight)-1)] for weight in L]
        L=[tuple(x) for x in L]
        return Set(L)

    def intersection_set(self, I_i, ray_i, monomial_0):
        """
        Returns a list of monomials in ``I_i`` whose pairing with ``ray_i`` is larger than the
        pairing with ``monomial_0``.
        """
        returning_list=list()
        for monomial in I_i:
            if self.group.pairing(ray_i, monomial)>self.group.pairing(ray_i, monomial_0):
                returning_list.append(monomial)
        return Set(returning_list)

    def compute_weights_in_all_unstable_states(self):
        """
        Computes the Set of weights that are present in all unstable states and stores it in the object. In most cases this set may be empty.
        
        EXAMPLES::
                
            sage: from Git import GITProblem 
            sage: Phi = WeylCharacterRing("A2")
            sage: representation= Phi(3,0,0)
            sage: P=GITProblem(representation,label="Plane cubics")
            sage: P.compute_weights_in_all_unstable_states()

        """
        weights=list(self.weights)
        #lengths gives a Euclidean norm to the weights, so that we have a way of choosing one.
        lengths=[sum(list(map(lambda i : i * i, weight)))  for weight in weights]
        intersection=self.weights_in_all_nonstable_states

        while weights!=list():
            good_char=False
            #the first loop chooses the weight in the intersection with minimum norm. In the meanwhile, it pops out from the list any weights that are not in the intersection.
            while good_char is False and weights!=list():
                weight_index=min(enumerate(lengths), key=itemgetter(1))[0]
                weight=weights.pop(weight_index)
                lengths.pop(weight_index)
                if weight in intersection:
                    good_char=True
            if good_char:
                J_weight=intersection
                for i in range(self.rank):
                    #Note that since self.intersection_set *only* considers strict inequality, weight itself is not in J_weight_i, thus weight is not in intersection.
                    J_weight_i=self.intersection_set(intersection, self.OPS_rays_list[i], weight)
                    J_weight=J_weight.intersection(J_weight_i)
                intersection=intersection.difference(J_weight)
        self.weights_in_all_unstable_states=self.weights_in_all_nonstable_states.difference(intersection)
        
    def H_dual_coordinates(self, weight):
        """
        Returns ``weight`` in H-dual coordinates, taking in account the type of group of GITProblem.
        
        This method returns H-dual coordinates, on the hom-spaces Hom(T, GG_m) of characters.  These are dual to H-coordinates on the hom-spaces Hom(GG_m , T) of one parameter subgroups, defined by matrices H_i with only one non-zero element (i, i) of unitary size.
        H and T are two bases to express one-parameter subgroups (see ``SimpleGroup``
        documentation for details). For most groups, H-coordinates and T-coordinates are equal,
        but in groups of type ``A`` they differ. As a result, the H-dual coordinates of a weight
        will depend on the group of the representation they live in.   
        
        INPUT:
        
        - ``weight`` -- The weight whose H-coordinates we require, in T-coordinates.
        
        EXAMPLES::
                    
            sage: from Git import GITProblem 
            sage: Phi = WeylCharacterRing("A2")
            sage: representation= Phi(3,0,0)
            sage: P=GITProblem(representation,label="Plane cubics")
            sage: P.H_dual_coordinates((2,1))
            (2, 1, 0)
        """
        if self.Dynkin_type=='A':
            return self.L_coord_to_H_dual_conversion[weight]
        else:
            return weight

    def generate_optimal_weights_non_stable(self):
        """
        This method creates a set of optimal weights to consider to apply Algorithm 3.7 in
        [GMGMS] (Set A_3 in said Algorithm). This method stors that set in
        attribute optimized_weights_non_stable.
        
        EXAMPLES::
        
            sage: from Git import GITProblem
            sage: Phi = WeylCharacterRing("A2")
            sage: representation= Phi(3,0,0)
            sage: P=GITProblem(representation)
            sage: P.generate_optimal_weights_non_stable()

        """
        zero_weight_set=Set(tuple([tuple([0 for i in range(self.rank)])]))#This is set {0} in Algorithm 3.7 in [GMGMS]
        first_optimization=self.nonstable_weights_candidates.difference(zero_weight_set) #This is set A2 in Algorithm 3.7 in [GMGMS] #WARNING: This is a Python set, not a SAGE set
        second_optimization=set([]) #WARNING: This is a Python set, not a SAGE set
        
        #These are lines 5-9 in Algirthm 3.7 in [GMGMS]
        for candidate in first_optimization:
            good = True
            for element in second_optimization:
                if proportional(candidate,element):
                    good = False
                    break
            if good:
                second_optimization.add(candidate)
        #optimized_weights_non_stable is set A3 in Algorithm 3.7 in [GMGMS]
        self.optimized_weights_non_stable = Set(second_optimization)
        
    def destabilized_weights(self, OPS, all_weights_considered=False, strict_inequality=False, nonstable_weights_considered=True):
        """
        Given a one-parameter subgroup (OPS) lambda = t^v determined by a
        vector ``OPS = v = (v_1, ... v_n)``, it returns the weights destabilized by its group action

        INPUT:
        
        - ``OPS`` -- the one-parameter subgroup.
        - ``all_weights_considered`` --
        - ``strict_inequality`` --
        - ``nonstable_weights_considered`` -- If ``True``, it will find all weights in the representation which are destabilised by the one-parameter subgroup. If ``False``, it will consider the parameter ``nonstable_weights_considered`` to determine which weights to find.
        - ``strict_inequality`` -- If ``True`` it will only include those weights whose pairing with ``OPS`` is strictly positive. If ``False`` it will include those weights whose pairing with ``OPS`` is non-negative.
        -- ``nonstable_weights_considered`` -- If ``True``, it will only consider weights that are non-stable with respect to ``OPS``. If ``False`` it will consider weights taht are unstable with respect to ``OPS``.
        
        EXAMPLES::
            
            sage: from Git import GITProblem
            sage: Phi = WeylCharacterRing("A2")
            sage: representation= Phi(3,0,0)
            sage: P=GITProblem(representation)
            sage: P.destabilized_weights(1,2,-3)
            {(1, -1), (-1, -2), (2, 1), (3, 0)}
        """
        if all_weights_considered:
            weights_considered = self.weights
        elif nonstable_weights_considered:
            if self.nonstable_weights_candidates is None:
                print('ERROR: nonstable_weights_candidates not yet computed')
                return None
            weights_considered=self.nonstable_weights_candidates
        else: #unstable weights
            weights_considered=self.unstable_weights_candidates
        nonstable_state = []
        for weight in weights_considered:
            if strict_inequality:
                condition=self.group.pairing(OPS, weight)>0
            else:
                condition=self.group.pairing(OPS, weight)>=0
            if condition:
                nonstable_state.append(weight)
        return Set(nonstable_state)

        

    def solve_non_stable(self, Weyl_optimisation=False):
        """
        Returns (and stores in the object) the non-stable locus of the group action with respect
        to a fixed torus, in the format {{(a_1, ..., a_r), ... }}, where a_1 is the non-stable of weight associated to the first coordinate of weight space, and so on. 

        INPUT:
        
        - ``Weyl_optimisation`` -- If ``True`` it will only consider one-parameter subgroups
        within the fundamental chamber rather than the whole lattice, potentially reducing the final output by eliminating isomorphic unstable/non-stable families. Note that as of
        January 2025, it is unknown whether this really reduces the number of maximal families
        (see Conjecture 7.4 in [GMGMS]) as the output is the same with or without Weyl
        optimisation.
        
        EXAMPLES::
            
            sage: from Git import GITProblem
            sage: Phi = WeylCharacterRing("A2")
            sage: representation= Phi(3,0,0)
            sage: P=GITProblem(representation)
            sage: P.solve_non_stable()
            {{(1, 2), (2, 1), (0, 0), (-1, 1), (0, 3), (1, -1), (3, 0)}, {(1, 2), (-1, -2), (2, 1), (0, 0), (1, -1), (3, 0)}}
        """
        #Reduce the number of sets of weights to be considered (set Q)
        if self.optimized_weights_non_stable is None:
            self.generate_optimal_weights_non_stable()
        candidate_weights_subsets = Set(list(self.optimized_weights_non_stable.subsets(self.rank-1)))
        
        #We find the maximal destabilised states
        maximal_nonstable_candidate_states = set() #WARNING: This is a Python set, not a Sage set. Needed for add/remove
                                                                                                            
        for candidate in candidate_weights_subsets:
            character_matrix=Matrix(QQ, weights_matrix(candidate))
            #Check if they have a unique solution and find it.
            M = character_matrix*self.group.fetch_pairing_matrix().transpose()

            M_kernel = M.right_kernel().basis()
            if len(M_kernel)==1: #The weights have one-dimensional solution
                #Check that the ray perpendicular to the set of weights 'candidate'
                #is in the Weyl fundamental chamber (and choose the right generator)
                gamma_OPS=one_param_subgroup(M_kernel[0], field=self.group.lattice_field())
                if self.group.in_cone(gamma_OPS):
                    destabilizing_OPS = gamma_OPS
                elif self.group.in_cone(-gamma_OPS):
                    destabilizing_OPS = (-1) * gamma_OPS
                else:
                    destabilizing_OPS = None

                if destabilizing_OPS != None:
                    destabilized_state=self.destabilized_weights(destabilizing_OPS)
                    
                    candidate_is_maximal = True
                    for currently_maximal_state in list(maximal_nonstable_candidate_states):
                        if destabilized_state.issubset(currently_maximal_state):
                            candidate_is_maximal = False
                            break
                        elif currently_maximal_state.issubset(destabilized_state):
                            maximal_nonstable_candidate_states.remove(currently_maximal_state)
                            self.gamma_OPS_nonstable_dictionary.pop(currently_maximal_state)
                    if candidate_is_maximal:
                        maximal_nonstable_candidate_states.add(destabilized_state)        # We find the maximal states among all the destabilised states
                        self.gamma_OPS_nonstable_dictionary[destabilized_state]=destabilizing_OPS
        
        
        #Add the weights that are nonstable and in every maximal state back into all maximal states
        enlarged_max_nonstable_states_list=list()
        for reduced_state in maximal_nonstable_candidate_states:
            enlarged_max_nonstable_states_list.append(reduced_state.union(self.weights_in_all_nonstable_states))
            OPS = self.gamma_OPS_nonstable_dictionary[reduced_state]
            self.gamma_OPS_nonstable_dictionary.pop(reduced_state)
            self.gamma_OPS_nonstable_dictionary[reduced_state.union(self.weights_in_all_nonstable_states)]=OPS
        self.unoptimized_maximal_nonstable_states=Set(enlarged_max_nonstable_states_list)
        

        # Perform optimisation step using the Weyl stabilisers
        
        if Weyl_optimisation:
            group_elements = self.Weyl_group()
            maximal_nonstable_final = set() #WARNING: This is a Python set, not a Sage set. Needed for add/remove
            maximal_nonstable_candidate_states_list_copy = list(self.unoptimized_maximal_nonstable_states)
            for candidate in list(self.unoptimized_maximal_nonstable_states):
                is_maximal = True
                lambda_ops=self.gamma_OPS_nonstable_dictionary[candidate]
                for g in group_elements:
                    for state in maximal_nonstable_candidate_states_list_copy:
                        lambda_ops_acted = one_param_subgroup(list(Matrix(self.group.lattice_field(), g.inverse())*(self.group.H_coordinates(lambda_ops))), type_A=self.Dynkin_type=="A", field=self.group.lattice_field())
                        acted_state=self.destabilized_weights(lambda_ops_acted, all_weights_considered=True)
                        if acted_state.issubset(state) and len(acted_state)!=len(state):
                            is_maximal=False
                            break
                    if not is_maximal:
                        break
                if is_maximal:
                    maximal_nonstable_final.add(candidate)
            self.maximal_nonstable_states=Set(list(maximal_nonstable_final))
        else:
            maximal_nonstable_final=set() #WARNING: This is a Python set, not a Sage set. Needed for add/remove
            maximal_nonstable_candidate_states_list_copy=list(self.unoptimized_maximal_nonstable_states)
            for candidate in list(self.unoptimized_maximal_nonstable_states):
                is_maximal = True
                lambda_ops=self.gamma_OPS_nonstable_dictionary[candidate]
                for state in maximal_nonstable_candidate_states_list_copy:
                    # lambda_ops_acted=g.action(self.group.H_coordinates(lambda_ops), type_A=self.Dynkin_type=="A")
                    acted_state=self.destabilized_weights(lambda_ops, all_weights_considered=True)
                    if acted_state.issubset(state) and len(acted_state)!=len(state):
                        is_maximal=False
                        break
                if is_maximal:
                    maximal_nonstable_final.add(candidate)
            self.maximal_nonstable_states=Set(list(maximal_nonstable_final))
        return self.maximal_nonstable_states
            
    def solve_unstable(self, Weyl_optimisation=False):
        """
        Returns (and stores in the object) the unstable locus of the group action with respect to a fixed torus, in the format {{(a_1, ..., a_r), ... }}, where a_1 is the non-stable of weight associated to the first coordinate of weight space, and so on. 
        
        INPUT:
        
        - ``Weyl_optimisation`` -- If ``True`` it will only consider one-parameter subgroups
        within the fundamental chamber rather than the whole lattice, potentially reducing the final output by eliminating isomorphic unstable/non-stable families. Note that as of
        January 2025, it is unknown whether this really reduces the number of maximal families
        (see Conjecture 7.4 in [GMGMS]) as the output is the same with or without Weyl
        optimisation.

        
        EXAMPLES::
            
            sage: from Git import GITProblem
            sage: Phi = WeylCharacterRing("A2")
            sage: representation= Phi(3,0,0)
            sage: P=GITProblem(representation)
            sage: P.solve_unstable()
            {{(1, 2), (2, 1), (0, 3), (1, -1), (3, 0)}}
        """
        candidate_weights_subsets = Set(list(self.unstable_weights_candidates.subsets(self.rank)))
        
        #We find the maximal unstable states
        unstable_states = []
        maximal_unstable_candidate_states = set() #WARNING: This is a Python set, not a Sage set. Needed for add/remove

        for candidate in candidate_weights_subsets:
            
            #Step carried out to make linear system to solve a homogeneous one
            candidate_list = list(candidate)
            substract_weight=candidate_list[len(candidate_list)-1]
            substract_matrix = Matrix(self.group.lattice_field(), [substract_weight for i in range(len(substract_weight)-1)])
            character_matrix = Matrix(self.group.lattice_field(), candidate_list[0:len(candidate_list)-1])
            linear_system_matrix=character_matrix-substract_matrix
            
            #Check if they have a unique solution and find it.
            M = linear_system_matrix*self.group.fetch_pairing_matrix().transpose()
            M_kernel = M.right_kernel().basis()
            if len(M_kernel)==1: #The weights have one-dimensional solution
                #Check that the ray perpendicular to the set of weights 'candidate'
                #is in the Weyl fundamental chamber (and choose the right generator)
                gamma_OPS = one_param_subgroup(M_kernel[0], field=self.group.lattice_field())
                pairing_value = self.group.pairing(gamma_OPS, substract_weight)
                if self.group.in_cone(gamma_OPS):
                    destabilizing_OPS=gamma_OPS
                elif self.group.in_cone(-gamma_OPS):
                    destabilizing_OPS = (-1) * gamma_OPS
                else:
                    destabilizing_OPS = None

                if destabilizing_OPS != None:
                    destabilized_state=self.destabilized_weights(destabilizing_OPS, all_weights_considered=False, strict_inequality=True,nonstable_weights_considered=False)
                    
                    candidate_is_maximal = True
                    for currently_maximal_state in list(maximal_unstable_candidate_states):
                        if destabilized_state.issubset(currently_maximal_state):
                            candidate_is_maximal = False
                            break
                        elif currently_maximal_state.issubset(destabilized_state):
                            maximal_unstable_candidate_states.remove(currently_maximal_state)
                            self.gamma_OPS_unstable_dictionary.pop(currently_maximal_state)
                            if debug:
                                print('currently_maximal_state', currently_maximal_state, 'removed\n'); input('')
                    if candidate_is_maximal:
                        maximal_unstable_candidate_states.add(destabilized_state)        # We find the maximal states among all the destabilised states
                        self.gamma_OPS_unstable_dictionary[destabilized_state]=destabilizing_OPS
                        if debug:
                            print('destabilized_state', destabilized_state, 'added\n'); input('')


        
        #Add the weights that are unnstable and in every maximal state back into all maximal states
        enlarged_max_unstable_states_list=list()
        for reduced_state in maximal_unstable_candidate_states:
            enlarged_max_unstable_states_list.append(reduced_state.union(self.weights_in_all_unstable_states))
            OPS = self.gamma_OPS_unstable_dictionary[reduced_state]
            self.gamma_OPS_unstable_dictionary.pop(reduced_state)
            self.gamma_OPS_unstable_dictionary[reduced_state.union(self.weights_in_all_unstable_states)] = OPS
        self.unoptimized_maximal_unstable_states = Set(enlarged_max_unstable_states_list)

   
        
        if Weyl_optimisation:
            group_elements = self.Weyl_group()
            maximal_unstable_final = set() #WARNING: This is a Python set, not a Sage set. Needed for add/remove
            maximal_unstable_candidate_states_list_copy = list(maximal_unstable_candidate_states)
            if debug:
                print('unoptimized_maximal_unstable_states', self.unoptimized_maximal_unstable_states)
            for candidate in list(self.unoptimized_maximal_unstable_states):
                is_maximal = True
                lambda_ops=self.gamma_OPS_unstable_dictionary[candidate]
                if debug:
                    print('CANDIDATE', candidate); input('')
                for g in group_elements:
                    for state in maximal_unstable_candidate_states_list_copy:
                        lambda_ops_acted = one_param_subgroup(list(g.inverse()*(self.group.H_coordinates(lambda_ops))), type_A=self.Dynkin_type=="A", field=self.group.lattice_field())
                        if debug:
                            print('lambda_ops_acted', lambda_ops_acted); input('');
                        acted_state=self.destabilized_weights(lambda_ops_acted, all_weights_considered=True)
                        if debug:
                            print('acted_state', acted_state, '\nstate', state); input('');
                        if acted_state.issubset(state) and len(acted_state)!=len(state):
                            is_maximal = False
                            if debug:
                                print('is_maximal', is_maximal)
                            if debug:
                                print('acted_state', acted_state, 'state', state)
                            input('')
                            break
                    if not is_maximal:
                        break
                if is_maximal:
                    if debug:
                        print('added', candidate); input('')
                    maximal_unstable_final.add(candidate)
            self.maximal_unstable_states=Set(list(maximal_unstable_final))
        else:
            maximal_unstable_final = set() #WARNING: This is a Python set, not a Sage set. Needed for add/remove
            maximal_unstable_candidate_states_list_copy=list(maximal_unstable_candidate_states)
            for candidate in list(self.unoptimized_maximal_unstable_states):
                is_maximal = True
                lambda_ops=self.gamma_OPS_unstable_dictionary[candidate]
                for state in maximal_unstable_candidate_states_list_copy:
                    # lambda_ops_acted=g.action(self.group.H_coordinates(lambda_ops), type_A=self.Dynkin_type=="A")
                    acted_state = self.destabilized_weights(lambda_ops, all_weights_considered=True)
                    if acted_state.issubset(state) and len(acted_state)!=len(state):
                        is_maximal = False
                        break
                if is_maximal:
                    maximal_unstable_final.add(candidate)
            self.maximal_unstable_states=Set(list(maximal_unstable_final))
            self.maximal_unstable_states=Set(list(self.unoptimized_maximal_unstable_states))
        return self.maximal_unstable_states
        
        
# solve_strictly_polystable is NOT YET FULLY IMPLEMENTED
# Need to add Weyl group optimization

    def solve_strictly_polystable(self):
        """
        Returns the strictly polystable locus of the group action with respect to a fixed torus,
        in the format {{(a_1, ..., a_r), ... }}, where a_1 is the non-stable of weight associated to the first coordinate of weight space, and so on. Note solve_non_stable() and solve_unstable() must have been called first.
        
        EXAMPLES::
            
            sage: from Git import GITProblem
            sage: Phi = WeylCharacterRing("A2")
            sage: representation= Phi(3,0,0)
            sage: P=GITProblem(representation)
            sage: P.solve_non_stable(Weyl_optimisation=True)
            {{(1, 2), (2, 1), (0, 0), (-1, 1), (0, 3), (1, -1), (3, 0)}, {(1, 2), (-1, -2), (2, 1), (0, 0), (1, -1), (3, 0)}}
            sage: P.solve_non_stable()
            {{(1, 2), (2, 1), (0, 0), (-1, 1), (0, 3), (1, -1), (3, 0)}, {(1, 2), (-1, -2), (2, 1), (0, 0), (1, -1), (3, 0)}}
            sage: P.solve_unstable()
            {{(1, 2), (2, 1), (0, 3), (1, -1), (3, 0)}}
            sage: P.solve_strictly_polystable()
            {{(0, 0)}, {(-1, 1), (1, -1), (0, 0)}}  
        """
        # Do lines 2 and 3 of Alg. 3.27
        maximal_states = set()   #This is the set $P_{ps}^F$ in Alg. 3.27
        candidate_states = set() #This is the set $\mathcal{S}_p$ in Alg. 3.27
        # Do lines 4-7 of Alg. 3.27
        for state in self.maximal_nonstable_states:
            OPS = self.gamma_OPS_nonstable_dictionary[state]
            strictly_polystable_state = set() #This is $\Xi_{V,\lambda=0}$ in Alg. 3.27
            for element in state:
                if self.group.pairing(OPS, element)==0:
                    strictly_polystable_state.add(element)
            P=Polyhedron(vertices=list(strictly_polystable_state))
            if self.trivial_character in P.relative_interior():
                strictly_polystable_state = Set(list(strictly_polystable_state))
                maximal_states.add(strictly_polystable_state)
            if self.trivial_character in P:
                strictly_polystable_state = Set(list(strictly_polystable_state))
                candidate_states.add(strictly_polystable_state)
        # Do lines 8-12 of Alg. 3.27
        for state in candidate_states:
            P = Polyhedron(vertices=list(state))
            subsets_of_state=Subsets(state)
            for subset in subsets_of_state:
                Q = Polyhedron(vertices=list(subset))
                if Q.dim() < P.dim():
                    if self.trivial_character in Q.relative_interior():
                        V = VectorSpace(QQ,self.rank)
                        span_of_subset = V.subspace(subset)
                        new_state = Set([v for v in self.weights if V(v) in span_of_subset])
                        maximal_states.add(new_state)
        # Do lines 14-16 of Alg. 3.27
        group_elements = self.Weyl_group()
        maximal_states_list = list(maximal_states)
        Wimages = set()
        for state in maximal_states_list:
            if state not in Wimages:
                for g in group_elements:
                    gstate = self.weyl_elt_action_on_state(g,state)
                    if gstate != state and gstate in maximal_states:
                        #print('When state='+str(state)+' and g='+str(g)+', get gstate='+str(gstate)+'. Removing '+str(gstate))
                        Wimages.add(gstate)
                        maximal_states.remove(gstate)
        self.strictly_polystable_states=Set(list(maximal_states))
        return self.strictly_polystable_states

    def print_solution_nonstable(self):
        """
        It prints the weights of maximal non-stable families with respect to a fixed torus. Note 
        solve_non_stable() must have been called first.
        
        EXAMPLES::
            
            sage: from Git import GITProblem
            sage: Phi = WeylCharacterRing("A2")
            sage: representation= Phi(3,0,0)
            sage: P=GITProblem(representation)
            sage: P.solve_non_stable()
            {{(1, 2), (2, 1), (0, 0), (-1, 1), (0, 3), (1, -1), (3, 0)}, {(1, 2), (-1, -2), (2, 1), (0, 0), (1, -1), (3, 0)}}
            sage: P.print_solution_nonstable()
            <BLANKLINE>           
            <BLANKLINE>            
            ***************************************
            SOLUTION TO GIT PROBLEM: NONSTABLE LOCI
            ***************************************
            Group: A2
            Representation  A2(3,0,0)
            Set of maximal non-stable states:
            (1) 1-PS = (1, 1, -2) yields a state with 7 characters
            Maximal nonstable state={ (1, 2, 0), (2, 1, 0), (1, 1, 1), (0, 2, 1), (0, 3, 0), (2, 0, 1), (3, 0, 0) }
            (2) 1-PS = (1, -1/2, -1/2) yields a state with 6 characters
            Maximal nonstable state={ (1, 2, 0), (1, 0, 2), (2, 1, 0), (1, 1, 1), (2, 0, 1), (3, 0, 0) }
        """
        if self.maximal_nonstable_states is None:
            print('ERROR: The problem is not yet solved. Call solve_non_stable() first and then call print_solution_nonstable()')
            return None
        print('\n\n***************************************\nSOLUTION TO GIT PROBLEM: NONSTABLE LOCI\n***************************************')
        print('Group: {s}{d}'.format(s=self.Dynkin_type, d=self.rank), sep='')
        print('Representation ', self.rep)
        print('Set of maximal non-stable states:')
        i=1
        for state in self.maximal_nonstable_states:
            print('({n}) 1-PS = '.format(n=i), self.group.H_coordinates(self.gamma_OPS_nonstable_dictionary[state]), ' yields a state with ', len(list(state)), ' characters', sep='')
            statelist = [self.H_dual_coordinates(element) for element in state];
            statestr = str(statelist)
            print('Maximal nonstable state={',statestr[1:-1],"}")
            #print('\n')
            i = i+1

        
    def solution_nonstable_str(self):
        """
        It returns a (very long) string describing the weights of maximal non-stable families
        with respect to a fixed torus. This may be useful to save it in a file. Note
        solve_non_stable() must have been called first.
                
        EXAMPLES::
            
            sage: from Git import GITProblem
            sage: Phi = WeylCharacterRing("A2")
            sage: representation= Phi(3,0,0)
            sage: P=GITProblem(representation)
            sage: P.solve_non_stable()
            {{(1, 2), (2, 1), (0, 0), (-1, 1), (0, 3), (1, -1), (3, 0)}, {(1, 2), (-1, -2), (2, 1), (0, 0), (1, -1), (3, 0)}}
            sage: P.solution_nonstable_str()
            '\n\n***************************************\nSOLUTION TO GIT PROBLEM: NONSTABLE LOCI\n***************************************\nGroup: A2 Representation A2(3,0,0)\nSet of maximal non-stable states:\n(1) 1-PS = (1, 1, -2) yields a state with 7 characters\nMaximal nonstable state={(1, 2, 0), (2, 1, 0), (1, 1, 1), (0, 2, 1), (0, 3, 0), (2, 0, 1), (3, 0, 0)}\n(2) 1-PS = (1, -1/2, -1/2) yields a state with 6 characters\nMaximal nonstable state={(1, 2, 0), (1, 0, 2), (2, 1, 0), (1, 1, 1), (2, 0, 1), (3, 0, 0)}\n'
        """
        if self.maximal_nonstable_states is None:
            return 'ERROR: The problem is not yet solved. Call solve_non_stable() first and then call print_solution_nonstable()'
        s = '\n\n***************************************\nSOLUTION TO GIT PROBLEM: NONSTABLE LOCI\n***************************************\n'
        s = s + 'Group: {s}{d}'.format(s=self.Dynkin_type, d=self.rank)
        s = s + ' Representation '+str(self.rep)+'\n'
        s = s + 'Set of maximal non-stable states:\n'
        i = 1
        for state in self.maximal_nonstable_states:
            s = s + '({n}) 1-PS = '.format(n=i)+str(self.group.H_coordinates(self.gamma_OPS_nonstable_dictionary[state]))+' yields a state with '+str(len(list(state)))+' characters\n'
            statelist = [self.H_dual_coordinates(element) for element in state];
            statestr = str(statelist)
            s = s + 'Maximal nonstable state={'+statestr[1:-1]+"}\n"
            i = i + 1
        return s

    
    def print_solution_unstable(self):
        """
        It prints the weights of maximal unstable families with respect to a fixed torus. Note 
        solve_non_unstable() must have been called first.

        EXAMPLES::
            
            sage: from Git import GITProblem
            sage: Phi = WeylCharacterRing("A2")
            sage: representation = Phi(3,0,0)
            sage: P=GITProblem(representation)
            sage: P.solve_unstable()
            {{(1, 2), (2, 1), (0, 3), (1, -1), (3, 0)}}
            sage: P.print_solution_unstable()
            <BLANKLINE>
            <BLANKLINE>            
            **************************************
            SOLUTION TO GIT PROBLEM: UNSTABLE LOCI
            **************************************
            Group: A2
            Representation  A2(3,0,0)
            Set of maximal unstable states:
            (1) 1-PS = (1, 1/4, -5/4) yields a state with 5 characters
            Maximal unstable state={ (1, 2, 0), (2, 1, 0), (0, 3, 0), (2, 0, 1), (3, 0, 0) }
        """
        if self.maximal_unstable_states is None:
            print('ERROR: The problem is not yet solved. Call solve_unstable() first')
            return None
        print('\n\n**************************************\nSOLUTION TO GIT PROBLEM: UNSTABLE LOCI\n**************************************')
        print('Group: {s}{d}'.format(s=self.Dynkin_type, d=self.rank))
        print('Representation ', self.rep)
        print('Set of maximal unstable states:')
        i=1
        for state in self.maximal_unstable_states:
            print ('({d}) 1-PS = '.format(d=i), self.group.H_coordinates(self.gamma_OPS_unstable_dictionary[state]), ' yields a state with ', len(list(state)), ' characters', sep='')
            statelist = [self.H_dual_coordinates(element) for element in state];
            statestr = str(statelist)
            print('Maximal unstable state={',statestr[1:-1],"}")
            #print('\n')
            i = i + 1

    def solution_unstable_str(self):
        """
        It returns a (very long) string describing the weights of maximal unstable families
        with respect to a fixed torus. This may be useful to save it in a file. Note
        solve_unstable() must have been called first.

        
        EXAMPLES::
            
            sage: from Git import GITProblem
            sage: Phi = WeylCharacterRing("A2")
            sage: representation= Phi(3,0,0)
            sage: P=GITProblem(representation)
            sage: P.solve_unstable()
            {{(1, 2), (2, 1), (0, 3), (1, -1), (3, 0)}}
            sage: P.solution_unstable_str()
            '\n\n**************************************\nSOLUTION TO GIT PROBLEM: UNSTABLE LOCI\n**************************************\nGroup: A2 Representation A2(3,0,0)\nSet of maximal unstable states:\n(1) 1-PS = (1, 1/4, -5/4) yields a state with 5 characters\nMaximal unstable state={(1, 2, 0), (2, 1, 0), (0, 3, 0), (2, 0, 1), (3, 0, 0)}\n'
        """
        if self.maximal_unstable_states is None:
            return 'ERROR: The problem is not yet solved. Call solve_unstable() first'
        s = '\n\n**************************************\nSOLUTION TO GIT PROBLEM: UNSTABLE LOCI\n**************************************\n'
        s = s + 'Group: {s}{d}'.format(s=self.Dynkin_type, d=self.rank)
        s = s + ' Representation '+str(self.rep)+'\n'
        s = s + 'Set of maximal unstable states:\n'
        i = 1
        for state in self.maximal_unstable_states:
            s = s + '({n}) 1-PS = '.format(n=i)+str(self.group.H_coordinates(self.gamma_OPS_unstable_dictionary[state]))+' yields a state with '+str(len(list(state)))+' characters\n'
            statelist = [self.H_dual_coordinates(element) for element in state];
            statestr = str(statelist)
            s = s + 'Maximal unstable state={'+statestr[1:-1]+"}\n"
            i = i + 1
        return s

            
    def print_solution_strictly_polystable(self):
        """
        It prints the weights of maximal strictly polystable (polystable but not stable) families
        with respect to a fixed torus. Note solve_non_unstable(), solve_unstable() and
        solve_strictly_polystable() must have been called first.

        
        EXAMPLES::
            
            sage: from Git import GITProblem
            sage: Phi = WeylCharacterRing("A2")
            sage: representation= Phi(3,0,0)
            sage: P=GITProblem(representation)
            sage: P.solve_non_stable(Weyl_optimisation=True)
            {{(1, 2), (2, 1), (0, 0), (-1, 1), (0, 3), (1, -1), (3, 0)}, {(1, 2), (-1, -2), (2, 1), (0, 0), (1, -1), (3, 0)}}
            sage: P.solve_non_stable()
            {{(1, 2), (2, 1), (0, 0), (-1, 1), (0, 3), (1, -1), (3, 0)}, {(1, 2), (-1, -2), (2, 1), (0, 0), (1, -1), (3, 0)}}
            sage: P.solve_unstable()
            {{(1, 2), (2, 1), (0, 3), (1, -1), (3, 0)}}
            sage: P.solve_strictly_polystable()
            {{(0, 0)}, {(-1, 1), (1, -1), (0, 0)}}  
            sage: P.print_solution_strictly_polystable()
            <BLANKLINE>            
            <BLANKLINE>            
            *************************************************************
            SOLUTION TO GIT PROBLEM: STRICTLY POLYSTABLE LOCI
            *************************************************************
            Group: A2
            Representation  A2(3,0,0)
            Set of strictly T-polystable states:
            (1) A state with 1 characters
            Strictly polystable state={ (1, 1, 1) }
            (2) A state with 3 characters
            Strictly polystable state={ (0, 2, 1), (2, 0, 1), (1, 1, 1) }
        """
        if self.strictly_polystable_states is None:
            print ('ERROR: The problem is not yet solved. Call solve_unstable() first. Then, call solve_strictly_polystable() and print_solution_strictly_polystable().')
            return None
        print ('\n\n*************************************************************\nSOLUTION TO GIT PROBLEM: STRICTLY POLYSTABLE LOCI\n*************************************************************')
        print ('Group: {s}{d}'.format(s=self.Dynkin_type, d=self.rank))
        print ('Representation ', self.rep)
        print ('Set of strictly T-polystable states:')
        i = 1
        for state in self.strictly_polystable_states:
            print ('({d}) '.format(d=i),'A state with ', len(list(state)), ' characters', sep='')
            statelist = [self.H_dual_coordinates(element) for element in state];
            statestr = str(statelist)
            print('Strictly polystable state={',statestr[1:-1],"}")
            #print('\n')
            i = i + 1

    def solution_strictly_polystable_str(self):
        """
        It returns a (very long) string describing the weights of maximal strictly polystable
        (polystable but not stable) families with respect to a fixed torus. This may be useful to
        save it in a file. Note solve_non_unstable(), solve_unstable() and
        solve_strictly_polystable() must have been called first.

        EXAMPLES::
            
            sage: from Git import GITProblem
            sage: Phi = WeylCharacterRing("A2")
            sage: representation= Phi(3,0,0)
            sage: P=GITProblem(representation)
            sage: P.solve_non_stable(Weyl_optimisation=True)
            {{(1, 2), (2, 1), (0, 0), (-1, 1), (0, 3), (1, -1), (3, 0)}, {(1, 2), (-1, -2), (2, 1), (0, 0), (1, -1), (3, 0)}}
            sage: P.solve_non_stable()
            {{(1, 2), (2, 1), (0, 0), (-1, 1), (0, 3), (1, -1), (3, 0)}, {(1, 2), (-1, -2), (2, 1), (0, 0), (1, -1), (3, 0)}}
            sage: P.solve_unstable()
            {{(1, 2), (2, 1), (0, 3), (1, -1), (3, 0)}}
            sage: P.solve_strictly_polystable()
            {{(0, 0)}, {(-1, 1), (1, -1), (0, 0)}}    
            sage: P.solution_strictly_polystable_str()
            '\n\n*************************************************************\nSOLUTION TO GIT PROBLEM: STRICTLY POLYSTABLE LOCI\n*************************************************************\nGroup: A2 Representation A2(3,0,0)\nSet of strictly polystable states:\n(1) A state with 1 characters\nStrictly polystable state={(1, 1, 1)}\n(2) A state with 3 characters\nStrictly polystable state={(0, 2, 1), (2, 0, 1), (1, 1, 1)}\n'
        """
        if self.strictly_polystable_states is None:
            return 'ERROR: The problem is not yet solved. Call solve_unstable() first. Then, call solve_strictly_polystable() and solution_strictly_polystable_str().'
        s = '\n\n*************************************************************\nSOLUTION TO GIT PROBLEM: STRICTLY POLYSTABLE LOCI\n*************************************************************\n'
        s = s + 'Group: {s}{d}'.format(s=self.Dynkin_type, d=self.rank)
        s = s + ' Representation '+str(self.rep)+'\n'
        s = s + 'Set of strictly polystable states:\n'
        i = 1
        for state in self.strictly_polystable_states:
            s = s + '({n}) '.format(n=i)+'A state with '+str(len(list(state)))+' characters\n'
            statelist=[self.H_dual_coordinates(element) for element in state];
            statestr=str(statelist)
            s = s + 'Strictly polystable state={'+statestr[1:-1]+"}\n"
            i = i + 1
        return s

            

    def print_solution(self):
        """
        It prints the weights of maximal unstable, non-stable and strictly polystable (polystable but not stable) families with respect to a fixed torus. Note solve_non_unstable(), solve_unstable() and solve_strictly_polystable() must have been called first.

        EXAMPLES::
            
            sage: from Git import GITProblem
            sage: Phi = WeylCharacterRing("A2")
            sage: representation= Phi(3,0,0)
            sage: P=GITProblem(representation)
            sage: P.solve_non_stable()
            {{(1, 2), (2, 1), (0, 0), (-1, 1), (0, 3), (1, -1), (3, 0)}, {(1, 2), (-1, -2), (2, 1), (0, 0), (1, -1), (3, 0)}}
            sage: P.solve_unstable()
            {{(1, 2), (2, 1), (0, 3), (1, -1), (3, 0)}}
            sage: P.solve_strictly_polystable()
            {{(0, 0)}, {(-1, 1), (1, -1), (0, 0)}}            
            sage: P.print_solution()
            <BLANKLINE>
            <BLANKLINE>
            ***************************************
            SOLUTION TO GIT PROBLEM: NONSTABLE LOCI
            ***************************************
            Group: A2
            Representation  A2(3,0,0)
            Set of maximal non-stable states:
            (1) 1-PS = (1, 1, -2) yields a state with 7 characters
            Maximal nonstable state={ (1, 2, 0), (2, 1, 0), (1, 1, 1), (0, 2, 1), (0, 3, 0), (2, 0, 1), (3, 0, 0) }
            (2) 1-PS = (1, -1/2, -1/2) yields a state with 6 characters
            Maximal nonstable state={ (1, 2, 0), (1, 0, 2), (2, 1, 0), (1, 1, 1), (2, 0, 1), (3, 0, 0) }
            <BLANKLINE>
            <BLANKLINE>
            **************************************
            SOLUTION TO GIT PROBLEM: UNSTABLE LOCI
            **************************************
            Group: A2
            Representation  A2(3,0,0)
            Set of maximal unstable states:
            (1) 1-PS = (1, 1/4, -5/4) yields a state with 5 characters
            Maximal unstable state={ (1, 2, 0), (2, 1, 0), (0, 3, 0), (2, 0, 1), (3, 0, 0) }
            <BLANKLINE>
            <BLANKLINE>
            *************************************************************
            SOLUTION TO GIT PROBLEM: STRICTLY POLYSTABLE LOCI
            *************************************************************
            Group: A2
            Representation  A2(3,0,0)
            Set of strictly T-polystable states:
            (1) A state with 1 characters
            Strictly polystable state={ (1, 1, 1) }
            (2) A state with 3 characters
            Strictly polystable state={ (0, 2, 1), (2, 0, 1), (1, 1, 1) }
        """
        self.print_solution_nonstable()
        self.print_solution_unstable()
        self.print_solution_strictly_polystable()

    def solve_all(self):
        """
        It stores in the object the non-stable and unstable locus of the group action with respect to a fixed torus. It does not use Weyl Optimisation.
        """
        self.solve_nonstable()
        self.solve_unstable()
