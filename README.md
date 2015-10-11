#Conjugate Gradients

###Simple Implementation of the Iterative Solver

See [this paper](https://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf) for a great explanation of the subject.

It's suitable for Symmetric Linear (Sparse) Problems.

The system matrix is modeled as a **List Of Lists**.

It features two simple preconditioners:

* [Jacobi](https://en.wikipedia.org/wiki/Preconditioner#Jacobi_.28or_diagonal.29_preconditioner)
* [SSOR](https://en.wikipedia.org/wiki/Successive_over-relaxation)

I used it as a solver for simple 2D [FEM](https://en.wikipedia.org/wiki/Finite_element_method) [Poisson Problems](https://en.wikipedia.org/wiki/Poisson's_equation).
