# Computational GIT

This software is a tool for computing GIT quotients. It is based on the paper “Computing GIT Quotients of Semisimple Groups” by Patricio Gallardo, Jesus Martinez Garcia, Han-Bom Moon, and David Swinarski. It can be used to calculate the stable, semistable and polystable locus of a polarised hypersurface equipped with the action of a classical algebraic group. For more information visit https://jesusmartinezgarcia.net/git/ or https://faculty.fordham.edu/dswinarski/ComputationalGIT/. 

# Dependencies 

 - python
 - sage
 - sage.geometry
 - sage.modules

# Outputs 

Each problem output will state the Dynkin type (root system) of a group 𝐺 as Xn where X is A,B,C,D,E,F,G and n is a positive integer. For instance, A3 corresponds to SL3. The output will also state the representation the group acts on by listing the highest weight(s) of the representation. For instance, A3(3,0,0,0) means the group is SL3 acting on a 4-dimensional weight system whose highest weight is 3𝜔_1.

The output also presents a list of non-stable, unstable and strictly polystable loci. Essentially, this is a list of one-parameter subgroups in 𝐺. Note these are presented just by their weights up to multiplication by scalar. Note that groups of type 𝐴 and groups of type 𝐵−𝐸 follow different basis and thus the weights of a one-parameter subgroup of type 𝐴 add up to 0 while the rest do not. We hope this does not cause confusion.

For each one-parameter subgroup, a state is listed. This state will contain all the weights of the representation which are non-stable, unstable or strictly polystable with respect to the one-parameter subgroup. The program (and the results in the paper) guarantee that any 𝑇-non-stable (or unstable, strictly polystable, respectively) point in 𝑋 must belong to one of these states. An example for cubic surfaces (whose output is listed below) is worked out in the paper. In it, one can read how to interpret this output to find all stable and strictly polystable points. 

More examples can be found at https://jesusmartinezgarcia.net/git/ or https://faculty.fordham.edu/dswinarski/ComputationalGIT/ 
