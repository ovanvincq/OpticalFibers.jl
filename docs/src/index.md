# OpticalFibers.jl

[![Stable](https://img.shields.io/badge/docs-dev-blue.svg)](https://ovanvincq.github.io/OpticalFibers.jl) [![DOI](https://zenodo.org/badge/696327811.svg)](https://zenodo.org/badge/latestdoi/696327811)

OpticalFibers is a package that allows to compute modes of optical fibers. Different methods are implemented to find scalar or vector modes:
- A semi-analytical solver (based on Bessel functions) for multi-step index fibers.
- Finite element method (using `Gridap.jl`) for any kind of isotropic or anisotropic fiber (useful to find leaky modes using a PML for example)

## Installation
OpticalFibers requires at least julia 1.9 and can be installed with:

```julia
using Pkg
Pkg.add("OpticalFibers")
```

## Quickstart
Computation of the scalar fundamental mode (l=0) of a step index fiber with a core-radius of 2 µm, a refractive index of 1.47 for core and 1.45 for cladding at a wavelength of 1 µm:
```julia
julia> using OpticalFibers
julia> using OpticalFibers.ModeSolvers
julia> ms=multi_step_fiber_modes(1,0,2,[1.47,1.45])
1-element Vector{Mode}:
 [LP 0,1,1.463179347605715,1,Nothing]
```
Computation of the fundamental vector mode (l=1) of the same fiber:
```julia
julia> mv=multi_step_fiber_modes(1,1,2,[1.47,1.45],type=:Vector)
1-element Vector{Mode}:
 [HE 1,1,1.4631371608572663,1,Nothing]
```

Computation of the scalar modes of a parabolic-index fiber with a core-radius of 4 µm, a refractive index of 1.48 for core center and 1.45 for cladding at a wavelength of 1 µm by using the finite element method with 1000 nodes between r=0 and r=20 µm:
```julia
julia> using OpticalFibers
julia> using OpticalFibers.ModeSolvers
julia> using Gridap
julia> m=FEM1D(1,0,x->(1.45+0.03*(1-x[1]^2/16)*(x[1]<=4))^2,CartesianDiscreteModel((0,20),1000),field=true,neigs=5)
2-element Vector{Mode{ScalarFieldFEM1D}}:
 [Mode LP n°1,1.471980656971672,1,ScalarFieldFEM1D]
 [Mode LP n°2,1.4561502566053002,1,ScalarFieldFEM1D]

julia> using Plots
julia> r=0:0.01:10
julia> plot(r,abs.(computeField(m[1],r)),label="LP01",xlabel="r (µm)",ylabel="|E|²")
```
![Fundamental mode example](docs/src/assets/fig1.png)

## Credits
OpticalFibers.jl is maintained by Olivier Vanvincq ([University of Lille](https://www.univ-lille.fr/), [PhLAM laboratory](https://phlam.univ-lille.fr/)).