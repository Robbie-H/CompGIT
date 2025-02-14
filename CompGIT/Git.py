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
    Decide if two non-zero vectors are proportional.
    
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
    Take a set of weights and return it as a matrix.
    
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


def timedRunProblem(representation,label='', separateOutputs=False):
    P=GITProblem(representation,label=label)
    t0=time()
    P.solve_non_stable(Weyl_optimisation=True)
    t1=time()
    if separateOutputs:
        f = open(P.label+' output.txt', 'a')
        s1=P.solution_nonstable_str()
        f.write(s1)
        f.close()
    t2=time()
    P.solve_unstable(Weyl_optimisation=True)
    t3=time()
    if separateOutputs:
        f = open(P.label+' output.txt', 'a')
        s2=P.solution_unstable_str()
        f.write(s2)
        f.close()
    t4=time()
    P.solve_strictly_polystable()
    t5=time()
    s0=str([P.label,P.rep,round(t1-t0,3),round(t3-t2,3),round(t5-t4,3),len(P.weights),len(P.optimized_weights_non_stable),len(P.unstable_weights_candidates),len(P.maximal_nonstable_states),len(P.maximal_unstable_states),len(P.strictly_polystable_states)])
    if separateOutputs:
        f = open(P.label+' output.txt', 'a')
        s3=P.solution_strictly_polystable_str()
        f.write(s3)
        f.write("\n")
        f.write(s0)
        s0=str([P.label,P.rep,round(t1-t0,3),round(t3-t2,3),round(t5-t4,3),len(P.weights),len(P.optimized_weights_non_stable),len(P.unstable_weights_candidates),len(P.maximal_nonstable_states),len(P.maximal_unstable_states),len(P.strictly_polystable_states)])
        f.close()
    if not separateOutputs:
        f = open(P.label+' output.txt', 'w')
        s1=P.solution_nonstable_str()
        f.write(s1)
        s2=P.solution_unstable_str()
        f.write(s2)
        s3=P.solution_strictly_polystable_str()
        f.write(s3)
        f.write("\n")
        s0=str([P.label,P.rep,round(t1-t0,3),round(t3-t2,3),round(t5-t4,3),len(P.weights),len(P.optimized_weights_non_stable),len(P.unstable_weights_candidates),len(P.maximal_nonstable_states),len(P.maximal_unstable_states),len(P.strictly_polystable_states)])
        f.write(s0)
        f.close()
    print('Timed run complete for '+str(representation))





class GITProblem(object):
    """
    
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
        self.OPS_rays_list=([one_param_subgroup(tuple(self.fundamental_chamber_generators[:,i].transpose())[0]) for i in range(0,self.rank)])
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
            
        EXAMPLES::
                
            sage: from Git import GITProblem 
            sage: Phi = WeylCharacterRing("A2")
            sage: representation= Phi(3,0,0)
            sage: P=GITProblem(representation,label="Plane cubics")
            sage: print(GITProblem.Weyl_group(P))
            Weyl Group of type ['A', 2] (as a matrix group acting on the ambient space)
        """
        return self.group.Weyl_Group_elements()

    def weyl_elt_action_on_state(self,M,state):
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
        returning_list=list()
        for monomial in I_i:
            if self.group.pairing(ray_i, monomial)>self.group.pairing(ray_i, monomial_0):
                returning_list.append(monomial)
        return Set(returning_list)

    def compute_weights_in_all_unstable_states(self):
        """
            
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
        
        EXAMPLES::
        
            sage: from Git import GITProblem
            sage: Phi = WeylCharacterRing("A2")
            sage: representation= Phi(3,0,0)
            sage: P=GITProblem(representation)
            sage: P.generate_optimal_weights_non_stable()

        """
        zero_weight_set=Set(tuple([tuple([0 for i in range(self.rank)])]))
        first_optimization=self.nonstable_weights_candidates.difference(zero_weight_set) #WARNING: This is a Python set, not a SAGE set
        second_optimization=set([]) #WARNING: This is a Python set, not a SAGE set
        for candidate in first_optimization:
            good = True
            for element in second_optimization:
                if proportional(candidate,element):
                    good = False
                    break
            if good:
                second_optimization.add(candidate)
        self.optimized_weights_non_stable=Set(second_optimization)
        
    def destabilized_weights(self, OPS, all_weights_considered=False, strict_inequality=False, nonstable_weights_considered=True):
        """
        
        EXAMPLES::
            
            sage: from Git import GITProblem
            sage: Phi = WeylCharacterRing("A2")
            sage: representation= Phi(3,0,0)
            sage: P=GITProblem(representation)
            sage: P.destabilized_weights(1)
            {(1, -1), (-1, -2), (0, 0)}
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
                gamma_OPS=one_param_subgroup(M_kernel[0])
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
                        lambda_ops_acted = one_param_subgroup(list(g.inverse()*(self.group.H_coordinates(lambda_ops))), type_A=self.Dynkin_type=="A")
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
            substract_matrix = Matrix(QQ, [substract_weight for i in range(len(substract_weight)-1)])
            character_matrix = Matrix(QQ, candidate_list[0:len(candidate_list)-1])
            linear_system_matrix=character_matrix-substract_matrix
            
            #Check if they have a unique solution and find it.
            M = linear_system_matrix*self.group.fetch_pairing_matrix().transpose()
            M_kernel = M.right_kernel().basis()
            if len(M_kernel)==1: #The weights have one-dimensional solution
                #Check that the ray perpendicular to the set of weights 'candidate'
                #is in the Weyl fundamental chamber (and choose the right generator)
                gamma_OPS = one_param_subgroup(M_kernel[0])
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
                    if candidate_is_maximal:
                        maximal_unstable_candidate_states.add(destabilized_state)        # We find the maximal states among all the destabilised states
                        self.gamma_OPS_unstable_dictionary[destabilized_state]=destabilizing_OPS

        
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
            for candidate in list(self.unoptimized_maximal_unstable_states):
                is_maximal = True
                lambda_ops=self.gamma_OPS_unstable_dictionary[candidate]
                for g in group_elements:
                    for state in maximal_unstable_candidate_states_list_copy:
                        lambda_ops_acted = one_param_subgroup(list(g.inverse()*(self.group.H_coordinates(lambda_ops))), type_A=self.Dynkin_type=="A")
                        acted_state=self.destabilized_weights(lambda_ops_acted, all_weights_considered=True)
                        if acted_state.issubset(state) and len(acted_state)!=len(state):
                            is_maximal = False
                            break
                    if not is_maximal:
                        break
                if is_maximal:
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
        self.solve_nonstable()
        self.solve_unstable()
