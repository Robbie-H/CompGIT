# Computational GIT

The CompGIT package is a tool for computing Geometric Invariant Theory (GIT) quotients in algebraic geometry. In a nutshell, GIT is a theory to model orbit spaces of algebraic varieties. Given an action of a simple complex reductive group $G$ on a projective space $\mathbb{P}^n$, CompGIT gives a description of the $G$-orbits of $\mathbb{P}^n$, called unstable/non-stable/strictly polystable orbits, that need to be removed/treated specially to form a well-behaved quotient. 

## Paradigmatic example of GIT quotient
As an example of the kind of problem being considered, let the 1-dimensional algebraic group of non-zero complex numbers $G=\mathbb C^*$, acting on two dimensional complex space $\mathbb C^2$ by $\lambda · (x, y) = (x\lambda^{-1}, \lambda y)$. This affine example is not considered by our code, but described here for pedagogical purposes. The orbits of the action are:
* Conics parametrised by $xy = c$, $c\in \mathbb C^*$. These are stable orbits (they are closed and their stabilisers are finite).
* the origin $(x, y) = (0, 0)$. This is a strictly polystable orbit: it is semistable and it is closed in the set of semistable orbits, but it is not stable as its stabiliser is not finite (it is $\mathbb C^* $).
* The punctured axes $\{x = 0, y\neq 0\}$  $\{y = 0, x\neq 0\}$. These are strictly semistable orbits which are neither polystable nor stable (they are not closed).

We have a quotient map $\mathbb C^2 \rightarrow \mathbb C^2//\mathbb C^* \cong \mathbb C$ to the GIT quotient, where the fibre at any $c\neq 0$ is the smooth conic $xy=c$. At $0\in \mathbb C$, the fibre is the union of the two punctured axes and the origin. In this example, the quotient parametrises polystable orbits, whereas projective GIT quotients also have unstable orbits that need to be removed when taking GIT quotients (unstable orbits on $V$ are the orbits that contain $0\in V$ in their closure).

## References
Several simplifications and observations can be carried out to reduce the problem to analysing torus actions on projective space. For a description of these, and the algorithms used, we refer the user to the paper [Computing GIT Quotients of Semisimple Groups](https://arxiv.org/abs/2308.08049) by Patricio Gallardo, Jesus Martinez Garcia, Han-Bom Moon, and David Swinarski, to which this code is a companion of and where several examples are also discussed. Further modifications to the code, its packaging and documentation have been led by Robert Hanson, with help by Martinez-Garcia and the Sagemath community, most notably Fréderic Chapoton. 

For a quick introduction to GIT we recommend [notes by Richard Thomas](https://arxiv.org/abs/math/0512411).

The aforementioned paper by Gallardo et.al. considers the more general case of projective varieties by a reductive group. However, in most applications the group is simple or semisimple. This has the advantage that we can use existing libraries in Sagemath for simple groups.

# Dependencies 

* SageMath version 9.0 or later and the SageMath standard library, 
* Python version 3.9 or later

The Sage installation guide is available [here](https://doc.sagemath.org/html/en/installation/index.html)

Commands are run within a console, such as the SageMath Shell on linux or terminal on MacOS. Once a sage session is launched, we include the prompt ‘sage:’ in our code.  
 
# Installation with pip 

To add CompGIT to your SageMath installation, first download CompGIT as a ```.zip``` file from the [CompGIT GitHub homepage](https://github.com/Robbie-H/CompGIT) and then run the console command 

```
sage -pip install (location of CompGIT)
```

# Example 

To use CompGIT, start a sage session and type 

```
sage: from CompGIT import *
```

One could begin by choosing a simple group by specifying the Dynkin type and rank: 

```
sage: G=SimpleGroup("A", 2)
```

The following computations print some basic properties of the group: 

```
sage: G.Dynkin_type 
'A'
sage: G.max_torus_dim 
2
sage: G.cone_basis # rays of the fundamental chamber
[2 1]
[1 2]
```

The main output of CompGIT is the solution to a GIT problem. We fix a representation. In the case below, $SL_3$ acts on a 3 dimensional weight system with weights $(3,0,0)$.

```
sage: Phi = WeylCharacterRing("A2")
sage: representation= Phi(3,0,0)
sage: P=GITProblem(representation,label="Plane cubics")
sage: P.solve_non_stable(Weyl_optimisation=True)
{{(1, 2), (2, 1), (0, 0), (-1, 1), (0, 3), (1, -1), (3, 0)}, {(1, 2), (-1, -2), (2, 1), (0, 0), (1, -1), (3, 0)}}
sage: P.print_solution_nonstable()


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
```

# Outputs 

The output states the Dynkin type by the notation Xn, where X is either A, B, C, D, E, F or G, and n is a positive integer. In the example above, A2 corresponds to the special linear group $SL_3$. The output will also state the highest weight(s) of the group action. In the example above, $SL_3$ acts on a 3-dimensional weight system, with highest weight 3𝜔_1.

Outputs can also state a list of unstable and strictly polystable loci. Essentially, this is a list of one-parameter subgroups in 𝐺, presented by their weights up to multiplication by scalar. For each one-parameter subgroup, a state is listed. This state will contain all the weights of the representation which are non-stable, unstable or strictly polystable with respect to the one-parameter subgroup. The program (and the results in the paper) guarantee that any 𝑇-non-stable (or unstable, strictly polystable, respectively) point in 𝑋 must belong to one of these states. 

An example for cubic surfaces is worked out in the paper. In it, one can read how to interpret this output to find all stable and strictly polystable points. More examples can be found at https://jesusmartinezgarcia.net/git/ or https://faculty.fordham.edu/dswinarski/ComputationalGIT/ 

# Running doctests

To run doctests for CompGIT, extract the files from the ```.zip``` into some folder. Make sure that you have a system variable called ```PYTHONPATH``` and include in it the path to the subfolder with ```CompGIT``` in it (the one with the file ```GIT.py``` in it). To check if the system variable exist you can run 

```printenv PYTHONPATH```

If something appears but not the CompGIT directory, write

```export PYTHONPATH="${PYTHONPATH}:directory"```

where ```directory``` is replaced by the full CompGIT source path.

If when calling ```printenv``` above no paths are produced, then write

```export PYTHONPATH="directory"```

where ```directory``` is replaced by the full CompGIT source path.

If when calling ```printenv``` you do see the CompGIT source path, you do not have to do anything, continue to the next step.

Then, from the ```CompGIT``` subfolder (where ```GIT.py``` is stored) run

```
sage -t .
```

See also the sagemath documentation on doctests available [here](https://doc.sagemath.org/html/en/developer/doctesting.html). If doctests pass, you will receive the message 

```
----------------------------------------------------------------------
All tests passed!
----------------------------------------------------------------------
```

# License
CompGIT is distributed under the terms of the GNU General Public License v3.0. 
