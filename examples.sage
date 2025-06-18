from CompGIT import *

# Enter choice of simple group by specifying Dynkin type and rank 
G=SimpleGroup("A", 2)

# List Lie theoretic facts 
G.group_type() 
G.rnk() 
G.fundamental_chamber_generators() # rays of the fundamental chamber
G.fetch_pairing_matrix() # pairing between one parameter subgroups and characters 

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

print( '\n************\n\nTESTING CREATION OF A REPRESENTATION V\'=Sp4(2,0)\n')
Sp4 = WeylCharacterRing("C2")
representation3 = Sp4(2,0);
#####Sp_4 acting on C^2 with maximum weight (2,0).

print( representation3)
print( '\n************\n\nTESTING CREATION OF PROBLEM FOR V\'\n')
P3=GITProblem(representation3)
P3.solve_non_stable()
P3.print_solution_nonstable()
P3.solve_unstable()
P3.print_solution_unstable()
P3.solve_strictly_polystable()
P3.print_solution_strictly_polystable()

#########################################################

print( '\n************\n\nTESTING CREATION OF A REPRESENTATION V\'=SL4(5,0,0)\n')
SL4 = WeylCharacterRing("A3")
representation4 = SL4(5,0,0);
P4=GITProblem(representation4)
P4.solve_non_stable()
P4.print_solution_nonstable()
P4.solve_unstable()
P4.print_solution_unstable()

SO5 = WeylCharacterRing("B2")
representation5= SO5(2,1); print (representation)
print( representation5)
P5=GITProblem(representation5)
P5.solve_non_stable()
P5.print_solution_nonstable()
P5.solve_unstable()
P5.print_solution_unstable()


Sp4 = WeylCharacterRing("C2")
representation6= Sp4(3,1); print(representation)
P6=GITProblem(representation6)

print( '\n************\n\nTESTING CREATION OF PROBLEM FOR V\'\n')
P6.solve_non_stable()
P6.print_solution_nonstable()
P6.solve_unstable()
P6.print_solution_unstable()

print ('\n************\n\nTESTING CREATION OF A REPRESENTATION\n')
G = WeylCharacterRing("A2")
V = G(3,0,0)
representation7 = V.exterior_power(2)
print (representation7)
P7=GITProblem(representation7)
P7.solve_non_stable()
P7.print_solution_nonstable()
P7.solve_unstable()
P7.print_solution_unstable()
############################################################


