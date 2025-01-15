# Computational GIT

The CompGIT package is a tool for computing GIT quotients in algebraic geometry. Given an action of a classical group on a projective variety, CompGIT computes the certain orbits (called unstable/non-stable/strictly polystable) that need to be removed to form a well-behaved quotient. 

For instance, let the non-zero complex numbers C^* act on C^2 by λ · (x, y) = (x/λ, λy). The orbits are
conics parametrised by xy = c , the axes x = 0, y = 0, and the origin (x, y) = (0, 0). The output of
CompGIT will instruct the user to remove the unstable orbits given by the axes - these are precisely the
orbits that are not closed under the C^*-action. One then obtains a GIT quotient space (C^2−{axes})/C^* that is topologically equal to C.

The contents of our code performs several simplifications that reduce the problem to analysing torus actions on projective space. This procedure is based on the paper _Computing GIT Quotients of Semisimple Groups_ by Patricio Gallardo, Jesus Martinez Garcia, Han-Bom Moon, and David Swinarski. For an introduction to GIT theory we recommend notes of Richard Thomas available on the arXiv at https://arxiv.org/abs/math/0512411. 

# Dependencies 

 - Python
 - SageMath

# Outputs 

Each problem output will state the Dynkin type (root system) of a group 𝐺 as Xn where X is A,B,C,D,E,F,G and n is a positive integer. For instance, A3 corresponds to SL3. The output will also state the representation the group acts on by listing the highest weight(s) of the representation. For instance, A3(3,0,0,0) means the group is SL3 acting on a 4-dimensional weight system whose highest weight is 3𝜔_1.

The output also presents a list of non-stable, unstable and strictly polystable loci. Essentially, this is a list of one-parameter subgroups in 𝐺. Note these are presented just by their weights up to multiplication by scalar. Note that groups of type 𝐴 and groups of type 𝐵−𝐸 follow different basis and thus the weights of a one-parameter subgroup of type 𝐴 add up to 0 while the rest do not. We hope this does not cause confusion.

For each one-parameter subgroup, a state is listed. This state will contain all the weights of the representation which are non-stable, unstable or strictly polystable with respect to the one-parameter subgroup. The program (and the results in the paper) guarantee that any 𝑇-non-stable (or unstable, strictly polystable, respectively) point in 𝑋 must belong to one of these states. 

An example for cubic surfaces is worked out in the paper. In it, one can read how to interpret this output to find all stable and strictly polystable points. More examples can be found at https://jesusmartinezgarcia.net/git/ or https://faculty.fordham.edu/dswinarski/ComputationalGIT/ 
