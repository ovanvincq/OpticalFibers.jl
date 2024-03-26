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

## Finite element mode solvers
The FEM solvers are based on `Gridap.jl` and can compute modes of isotropic (functions `FEM1D` and `FEM2D`) and anisotropic fibers (function `FEM2D_general`).
If a PML is not used, the functions `FEM1D` and `FEM2D` compute guided modes only. To compute the leaky modes, the user must add a PML by setting the value of dPML in these functions. 

The computation of the effective index amounts to an eigenvalue problem which is solved by using the package `ArnoldiMethod.jl` combined with `LinearAlgebra.jl` (LU decomposition) or `MUMPS.jl`.
Note that in the case of anisotropic fibers, the dimension of the matrix is twice as large as in the case of isotropic fibers.

```@docs
    FEM1D
    FEM2D
    FEM2D_anisotropic
```
