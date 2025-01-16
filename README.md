# Computational GIT

The CompGIT package is a tool for computing Geometric Invariant Theory (GIT) quotients in algebraic geometry. In a nutshell, GIT is a theory to model orbit spaces of algebraic varieties into algebraic varieties. Given an action of a complex reductive connected simple Lie group $G$ on a projective space $\mathbb{P}^n$, CompGIT gives a description (up to $G$-equivalence) of orbits of $\mathbb{P}^n$ (called unstable/non-stable/strictly polystable orbits) that need to be removed/treated specially to form a well-behaved quotient of $\mathbb{P}^n$ by $G$. Because any $G$-linearised projective variety $X$ can be embedded $G$-equivariantly into some $\mathbb{P}^n$, this package is (potentially) sufficient to study any GIT quotient of a projective variety by a group $G$.

## Paradigmatic example of GIT quotient
As an example of the kind of problem being considered, let the 1-dimensional algebraic group of non-zero complex numbers $G=\mathbb C^*$ act on $\mathbb C^2$ by $\lambda · (x, y) = (x\lambda^{-1}, \lambda y)$. Note this is not the case contemplated by our code. Indeed, this group is not simple and $\mathbb C$ is not a projective variety. However, we describe this example for pedagogical purposes as the paradigm is the same, and this example is easier to understand. The orbits are:
* Conics parametrised by $xy = c$, $c\in \mathbb C^*$. These are stable orbits (they are closed and their stabilisers are finite).
* the origin $(x, y) = (0, 0)$. This is a strictly polystable orbit: it is semistable and it is closed in the set of semistable orbits, but it is not stable as its stabiliser is not finite (it is $\mathbb C^* $).
* The punctured axes $\{x = 0, y\neq 0\}$  $\{y = 0, x\neq 0\}$. These are strictly semistable orbits which are not polystable or stable (they are not closed).

We thus have an equivariant quotient map $\mathbb C^2 \rightarrow \mathbb C^2//\mathbb C^* \cong \mathbb C$, where the fibre at any $c\neq 0$ is the smooth conic $xy=c$. At $0\in \mathbb C$, the fibre is the union of the two punctured axes and the origin. Note the quotient does parametrise polystable orbits.

On a projective GIT quotient $\mathbb P(V)$, we would also have unstable orbits. These are orbits that need to be removed for any sense of quotient to make sense. Unstable orbits are the orbits that contain $0\in V$ in their closure.

## References
Several simplifications and observations can be carried out to reduce the problem to analysing torus actions on projective space. For a description of these, and the algorithms used, we refer the user to the paper [Computing GIT Quotients of Semisimple Groups](https://arxiv.org/abs/2308.08049) by Patricio Gallardo, Jesus Martinez Garcia, Han-Bom Moon, and David Swinarski, to which this code is a companion of and where several examples are also discussed. Further modifications to the code, its packaging and documentation have been led by Robert Hanson, with help by Martinez-Garcia and the Sagemath community, most notably Fréderic Chapoton. 

For a quick introduction to GIT we recommend [the notes by Richard Thomas](https://arxiv.org/abs/math/0512411).

The aforementioned paper by Gallardo et.al. considers the more general case of projective varieties by a reductive group. However, in most applications the group is always simple or semisimple. This has the advantage that we can use existing libraries in Sagemath for simple groups.

# Dependencies 

 - Python
 - SageMath

# Outputs 

Each problem output will state the Dynkin type (root system) of a group 𝐺 as Xn where X is A,B,C,D,E,F,G and n is a positive integer. For instance, A3 corresponds to SL3. The output will also state the representation the group acts on by listing the highest weight(s) of the representation. For instance, A3(3,0,0,0) means the group is SL3 acting on a 4-dimensional weight system whose highest weight is 3𝜔_1.

The output also presents a list of non-stable, unstable and strictly polystable loci. Essentially, this is a list of one-parameter subgroups in 𝐺. Note these are presented just by their weights up to multiplication by scalar. Note that groups of type 𝐴 and groups of type 𝐵−𝐸 follow different basis and thus the weights of a one-parameter subgroup of type 𝐴 add up to 0 while the rest do not. We hope this does not cause confusion.

For each one-parameter subgroup, a state is listed. This state will contain all the weights of the representation which are non-stable, unstable or strictly polystable with respect to the one-parameter subgroup. The program (and the results in the paper) guarantee that any 𝑇-non-stable (or unstable, strictly polystable, respectively) point in 𝑋 must belong to one of these states. 

An example for cubic surfaces is worked out in the paper. In it, one can read how to interpret this output to find all stable and strictly polystable points. More examples can be found at https://jesusmartinezgarcia.net/git/ or https://faculty.fordham.edu/dswinarski/ComputationalGIT/ 
