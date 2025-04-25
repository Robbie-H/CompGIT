from CompGIT import *

# # Enter choice of simple group by specifying Dynkin type and rank 
# G=SimpleGroup("A", 2)

# # List Lie theoretic facts 
# G.group_type() 
# G.rnk() 
# G.fundamental_chamber_generators() # rays of the fundamental chamber
# G.fetch_pairing_matrix() # pairing between one parameter subgroups and characters 

# Solve a GIT problem for plane cubics
Phi = WeylCharacterRing("A2")
representation= Phi(3,0,0) # SL_2 is acting on a 3 dimensional weight system with highest weight 3omega_1
P=GITProblem(representation,label="Plane cubics")
P.solve_non_stable(Weyl_optimisation=True)
P.print_solution_nonstable() 
P.solve_unstable(Weyl_optimisation=True)
P.print_solution_unstable()
P.solve_strictly_polystable()
P.print_solution_strictly_polystable()

# # Solve a GIT problem for cubic surfaces
# Phi = WeylCharacterRing("A3")
# representation= Phi(3,0,0, 0) # SL_3 is acting on a 3 dimensional weight system with highest weight 3omega_1
# P=GITProblem(representation,label="Plane cubics")
# P.solve_non_stable(Weyl_optimisation=True)
# P.print_solution_nonstable() 
# P.solve_unstable(Weyl_optimisation=True)
# P.print_solution_unstable()
# P.solve_strictly_polystable()
# P.print_solution_strictly_polystable()
# P.solve_non_stable(Weyl_optimisation=False)
# P.print_solution_nonstable() 
# P.solve_unstable(Weyl_optimisation=False)
# P.print_solution_unstable()
# P.solve_strictly_polystable()
# P.print_solution_strictly_polystable()

# # Basic operations on weights  
# weights=P.weights 
# averageWeight(weights) # yeilds zero for type A  
# weights_matrix(weights) # take a set of weights and return it as a matrix

# # compute the trivial character of the representation
# P.trivial_character 

# # Compute weights they may be in some destabilised state 
# P.nonstable_weights_destabilized
# P.unstable_weights_destabilized

# # Compute maximal non-stable states 
# P.maximal_nonstable_states

# Solve a GIT problem for G2
Phi2 = WeylCharacterRing("G2")
representation= Phi2(3,0,0) # G_2 is acting on a 7 dimensional weight system with highest weight 3omega_1
P2=GITProblem(representation,label="Weird cubics")
P2.solve_non_stable(Weyl_optimisation=True)
P2.print_solution_nonstable() 

# Solve a GIT problem for G2
Phi2 = WeylCharacterRing("G2")
representation= Phi2(3,0,0) # G_2 is acting on a 7 dimensional weight system with highest weight 3omega_1
P2=GITProblem(representation,label="Weird cubics")
P2.solve_non_stable(Weyl_optimisation=False)
P2.print_solution_nonstable() 

# # Solve a GIT problem for G2
# Phi2 = WeylCharacterRing("G2")
# representation= Phi2(3,0,0) # G_2 is acting on a 7 dimensional weight system with highest weight 3omega_1
# P2.solve_unstable(Weyl_optimisation=True)
# P2.print_solution_unstable()
# P2.solve_strictly_polystable()
# P2.print_solution_strictly_polystable()

# # Solve a GIT problem for G2
# Phi2 = WeylCharacterRing("G2")
# representation= Phi2(3,0,0) # G_2 is acting on a 7 dimensional weight system with highest weight 3omega_1
# P2.solve_unstable(Weyl_optimisation=False)
# P2.print_solution_unstable()
# P2.solve_strictly_polystable()
# # P2.print_solution_strictly_polystable()