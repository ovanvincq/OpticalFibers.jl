# OpticalFibers.ModeSolvers

```@meta
CurrentModule = OpticalFibers.ModeSolvers
```

This module contains different mode solvers for optical fibers. They are based on the interface conditions of analytical solutions between the layers of the fiber, finite difference method or finite element method.

## Fields
This module defines two types of field: the scalar field (valid in the case of weakly guiding fibers) and the vector field.
```@docs
Field
ScalarField
VectorField
ScalarFieldFEM
VectorFieldFEM
```
Two fields of the same type can be add or substract if their coordinates x and y are similar. It is also possible to multiply or divide a field by a scalar.

Two convenient function to interpolate the fields on another grid (changing x and y) are implemented:
```@docs
Interpolation(::ScalarField,::Vector{<:Real},::Vector{<:Real})
Interpolation(::VectorField,::Vector{<:Real},::Vector{<:Real})
```

## Modes
```@docs
Mode
ScalarMode1D
ScalarMode2D
VectorMode
ScalarModeFEM
VectorModeFEM
```
A ScalarMode1D can be convert into a ScalarMode2D or a VectorMode:
```@docs
ScalarMode2D(::ScalarMode1D;::char)
VectorMode(::ScalarMode1D;::Char,::Char)
```
A ScalarMode2D can be convert into a VectorMode:
```@docs
VectorMode(::ScalarMode2D;::Char)
```

A ScalarModeFEM can be convert into a VectorModeFEM:
```@docs
VectorModeFEM(::ScalarModeFEM;::Char)
```

## Propagation of modes
The field due to the propagation of a mode after a propagation distance z can be compute be using these functions:
```@docs
ScalarField(::ScalarMode1D,::Real;::Char)
ScalarField(::ScalarMode2D,::Real)
VectorField(::VectorMode,::Real)
VectorField(::ScalarMode2D,::Real;::Char)
VectorField(::ScalarMode1D,::Real;::Char,::Char)
ScalarFieldFEM(::ScalarModeFEM,::Real)
VectorFieldFEM(::VectorModeFEM,::Real)
```

## Poynting vector
Convenient functions are impleted to compute the poynting vector of modes or fields:
```@docs
PoyntingVector(::VectorField)
PoyntingVector(::VectorMode)
PoyntingVector(::ScalarMode2D)
PoyntingVector(::ScalarMode1D;::Char)
PoyntingVector(::ScalarModeFEM)
PoyntingVector(::VectorModeFEM)
```

## Mode normalization
The modes given by the mode solvers of this package are not normalized. You can normalize them by using the function `normalize!`.
The vector modes are normalized so as $\vert \iint \frac{1}{2}\left(\vec{E}\wedge\vec{H}^*\right) dS \vert=1$ and the scalar mode are normalized like vector modes if the optional parameter `unitIntegral` is `false` or so as $ \iint \vert E \vert^2 dS=1 $ if `unitIntegral` is `true`.
```@docs
normalize!(::ScalarMode1D;::Bool)
normalize!(::ScalarMode2D;::Bool)
normalize!(::VectorMode)
normalize!(::ScalarModeFEM;::Bool)
normalize!(::VectorModeFEM)
```

## Overlap integral
The overlap integral between two fields f$_1$ and f$_2$ is $\langle f_1 \vert f_2 \rangle = \iint E_1 E_2^* dS$ for scalar fields and $\langle f_1 \vert f_2 \rangle = \iint \frac{1}{2}\left(\vec{E}_1 \wedge \vec{H}_2^*\right)_z dS$ for vector fields.  

It is also possible to compute overlap integrals between two modes or between a mode and a field. However, when modes are involved, they are normalized before the calculation.

```@docs
overlap(::ScalarField,::ScalarField)
overlap(::VectorField,::VectorField)
overlap(::ScalarMode1D,::ScalarMode1D)
overlap(::ScalarMode2D,::ScalarMode2D)
overlap(::VectorMode,::VectorMode)
overlap(::ScalarField,::ScalarMode2D)
overlap(::ScalarMode2D,::ScalarField)
overlap(::VectorField,::VectorMode)
overlap(::VectorMode,::VectorField)
```

## Effective area
The effective area is defined by:
- for scalar modes: $\frac{\left( \iint \vert E \vert^2 dS \right)^2}{\iint \vert E \vert^4 dS}$
- for vector modes, there are two different values [Laegsgaard2012](@cite) : $\frac{\mu_0}{\varepsilon_0}\frac{\left(Re\left(\iint(\vec{E}\wedge\vec{H}^*)_z dS\right)\right)^2}{\iint n_0^2 \vert \vec{E}.\vec{E}^* \vert^2 dS}$ and $\frac{\mu_0}{\varepsilon_0}\frac{\left(Re\left(\iint(\vec{E}\wedge\vec{H}^*)_z dS\right)\right)^2}{\iint n_0^2 \left(\vec{E}.\vec{E}\right)\left(\vec{E}^*.\vec{E}^*\right) dS}$ where $n_0$ is the refractive index involved in the relation between the nonlinear index $n_2$ and the third order susceptibility $\chi^{(3)}$.

```@docs
    Aeff(::ScalarMode1D)
    Aeff(::ScalarMode2D)
    Aeff(::VectorMode,::Union{Real,Function})
    Aeff(::ScalarModeFEM)
    Aeff(::VectorModeFEM,::Union{Real,Function})
```

## Non-linear coefficient
If all the lentghs are in microns and $n_2$ is in SI (m²/W), you have to multiply the result by $10^{18}$ to obtain the result in W$^{-1}$.m$^{-1}$
```@docs
    nonLinearCoefficient(::ScalarMode1D,::Union{Real,Function})
    nonLinearCoefficient(::ScalarMode2D,::Union{Real,Function})
    nonLinearCoefficient(::VectorMode,::Union{Real,Function},::Union{Real,Function})
    nonLinearCoefficient(::ScalarModeFEM,::Union{Real,Function})
    nonLinearCoefficient(::VectorModeFEM,::Union{Real,Function},::Union{Real,Function})
```

## Mode Field Diameter (MFD)
The computation MFD is only valid for Gaussian-like beam (maximum at the center of the fiber and electric field with a constant sign). The MFD is calculated by finding the positions where $E=\frac{\max{\left(E\right)}}{\exp(1)}$ for scalar modes and $P_z=\frac{\max{\left(P_z\right)}}{\exp(2)}$ for vector modes.
```@docs
    MFD(::ScalarMode1D)
    MFD(::ScalarMode2D,::Real)
    MFD(::VectorMode,::Real)
    MFD(::ScalarModeFEM,::Real)
    MFD(::VectorModeFEM,::Real)
```
