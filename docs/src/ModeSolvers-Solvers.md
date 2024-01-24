# OpticalFibers.ModeSolvers - Solvers

```@meta
CurrentModule = OpticalFibers.ModeSolvers
```

## Quasi-analytical multi-step index fibers solver
This solver assumes that the fiber consists of several concentric layers of uniform refractive index. The cladding is assumed to be infinite. In each layer, the analytical solutions are the Bessel functions. The solver uses the interface conditions to predict the effective index and the profile of the modes. This method is described in the book written by J. Bures [Bures2009](@cite).  

This solver only returns guided modes: the effective index of the mode is real and cannot be lower than the refractive index of the outer cladding. 
```@docs
    multi_step_fiber_modes
```

## Finite difference mode solvers
In this package, the FD solvers only return guided modes. The method used in the vectorial case is described in the paper of Zhu [Zhu2002](@cite).

The computation of the effective index amounts to an eigenvalue problem. Three solutions are available to solve this eigenvalue problem: the use of the package `Arpack.jl` and the use of the package `ArnoldiMethod.jl` combined with `LinearAlgebra.jl` (LU decomposition) or `MUMPS.jl`.

For the moment, it is only possible to use a PML in the 1D scalar case.

```@docs
    FD
```

## Finite element mode solvers
The FEM solvers are based on `Gridap.jl` and can compute modes of isotropic and anisotropic fibers (useful when using a PML).

The computation of the effective index amounts to an eigenvalue problem. Two solutions are available to solve this eigenvalue problem: the use of the package `ArnoldiMethod.jl` combined with `LinearAlgebra.jl` (LU decomposition) or with `MUMPS.jl`.

In the case of anisotropic fibers, the dimension of the matrix is twice as large as in the case of isotropic fibers.

```@docs
    FEM
```
