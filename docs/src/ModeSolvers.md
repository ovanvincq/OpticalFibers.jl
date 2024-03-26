# OpticalFibers.ModeSolvers - Modes and Fields

```@meta
CurrentModule = OpticalFibers.ModeSolvers
```

## Fields
This module defines two types of fields: the scalar field (valid in the case of weakly guiding fibers) and the vector field. If the fiber has a cylindrical symmetry, it is sufficient to describe the variations of a scalar field along a radius (the 2D scalar field can be then obtained by multiplying the 1D field by $\cos(\ell \theta)$ or $\sin(\ell \theta)$ where $\ell$ is the azimuthal number).

Four abstract structures are defined to describe an electromagnetic field:
```@docs
Field
ScalarField1D
ScalarField2D
VectorField2D
```

Each kind of field can be described by a function, a matrix (a vector in the 1D case) or a `CellField` when using a FEM solver.

Note that functions are always function of a tuple `(r,)` in the case of a 1D field or `(x,y)` in case of a 2D field. This ensures a direct compatibility with the packages `Gridap`(FEM method) and `HCubature` (for integration).
```@docs
ScalarFieldFunction1D
ScalarFieldFunction2D
VectorFieldFunction2D
ScalarFieldMatrix1D
ScalarFieldMatrix2D
VectorFieldMatrix2D
ScalarFieldFEM1D
ScalarFieldFEM2D
VectorFieldFEM2D
```

The fields can be multiplied by a scalar and conjugated. It is also possible to take the real part or the imaginary part of a field.
Two fields of the same type can also be added or subtracted.

A function is implemented to compute the value of a field at a given position (r,) in 1D or (x,y) in 2D:
```@docs
computeField
```

## Modes
A mode is a structure that contains a name, an effective index, the wavelength and, optionally, a field.

This package assumes that the field is proportionnal to $\exp\left(i\beta z-\omega t\right)$.
```@docs
Mode
```

## Check if the field/mode is valid
A function was implemented to check if a field or a mode is valid by checking the number of arguments of the functions, the dimension of the matrices or the triangulaiton of the FEM fields. 
```@docs
isvalidField
isvalidMode
```

## Modes/Fields conversion
It is possible to convert a 1D function or matrix mode/field to 2d. The argument `angle` corresponds to the orientation in degrees (0 correspond to $\cos(\ell \theta)$ and -90 to $\sin(\ell \theta)$). This function is not available for FEM mode/field since it requires the creation of a 2D mesh.
```@docs
convertTo2D
```

It is also possible to convert a function or FEM field to a matrix mode/field.
```@docs
convertToMatrix
```

A mode that contains a ScalarField2D field can be converted to a mode with a VectorField2D with the approximation $\vec{k}=\beta \vec{e}_z$ and $E_z=H_z=0$.
```@docs
convertToVector
```

## Propagation of modes
The field due to the propagation of a mode after a propagation distance z can be calculated by using this function:
```@docs
getField
```

## Losses
The mode losses are computed from the imaginary part of the effective index.

The function returns the losses of a mode in dB/km if the wavelength is in meters. If the wavelength is in microns, you have to multiply the result by 1E6 to obtain the losses in dB/km. 
```@docs
losses(::Mode)
```

## Poynting vector
Convenient functions that compute the Poynting vector of 2D modes and fields. All functions return a tuple (Px,Py,Pz) that contains the three components of the Poynting vector in the same type as the field.
```@docs
PoyntingVector
```

## Mode normalization
The modes given by the mode solvers of this package are not normalized. You can normalize them by using the function `normalize!`.
The vector modes are normalized so as $\vert \iint \frac{1}{2}\left(\vec{E}\wedge\vec{H}^*\right) dS \vert=1$ and the scalar mode are normalized like vector modes if the optional parameter `unitIntegral` is `false` or so as $ \iint \vert E \vert^2 dS=1 $ if `unitIntegral` is `true`.
```@docs
normalize!
```

## Overlap integral
The overlap integral between two fields or modes f$_1$ and f$_2$ is $\langle f_1 \vert f_2 \rangle = \iint E_1 E_2^* dS$ for scalar fields/modes and $\langle f_1 \vert f_2 \rangle = \iint \frac{1}{2}\left(\vec{E}_1 \wedge \vec{H}_2^*\right)_z dS$ for vector fields/modes.  

```@docs
overlap
```

## Effective area
The effective area is defined by:
- for scalar modes: $\frac{\left( \iint \vert E \vert^2 dS \right)^2}{\iint \vert E \vert^4 dS}$
- for vector modes, there are two different values [Laegsgaard2012](@cite) : $\frac{\mu_0}{\varepsilon_0}\frac{\left(Re\left(\iint(\vec{E}\wedge\vec{H}^*)_z dS\right)\right)^2}{\iint n_0^2 \vert \vec{E}.\vec{E}^* \vert^2 dS}$ and $\frac{\mu_0}{\varepsilon_0}\frac{\left(Re\left(\iint(\vec{E}\wedge\vec{H}^*)_z dS\right)\right)^2}{\iint n_0^2 \left(\vec{E}.\vec{E}\right)\left(\vec{E}^*.\vec{E}^*\right) dS}$ where $n_0$ is the refractive index involved in the relation between the nonlinear index $n_2$ and the third order susceptibility $\chi^{(3)}$. If n0=0, the function assumes that n0≃neff (correct for weakly-guiding fibers).

```@docs
    Aeff
```

## Non-linear coefficient
If all the lengths are in microns and $n_2$ is in SI (m²/W), you have to multiply the result by $10^{18}$ to obtain the result in W$^{-1}$.m$^{-1}$.
In this function, n2 is the nonlinear index and n0 is the refractive index involved in the relation between $n_2$ and the third order susceptibility $\chi^{(3)}$. They must be constants or functions of the tuple (r,) in the 1D case or (x,y) in the 2D case. If n0=0, the function assumes that n0≃neff (correct for weakly-guiding fibers).
In the case of modes describes by functions, it is possible to specify the relative and the absolute tolerances for the integral calculations.
```@docs
    nonLinearCoefficient
```

## Mode Field Diameter (MFD)
The computation MFD is only valid for Gaussian-like beam (maximum at the center of the fiber and electric field with a constant sign). The MFD is calculated by finding the positions where $\vert E \vert=\frac{\max{\left(\vert E\vert \right)}}{\exp(1)}$ for scalar modes and $P_z=\frac{\max{\left(P_z\right)}}{\exp(2)}$ for vector modes. The optional parameters `theta` is the angle (in degrees) between the direction along which the MFD is computed and the x-axis.
If the field is not quasi-Gaussian (in the case of a higher order mode for example), the value given by this function has no physical meaning.
```@docs
    MFD
```

## Functions specific to FEM fields/modes
The function writevtk of the package Gridap is overloaded to save a FEM field or mode to a file that can be opened with ParaView.
```@docs
    writevtk
```
The function get_triangulation of the package Gridap is overloaded to get the triangulation of a FEM field (which is needed to plot a field with GridapMakie for example).
```@docs
    get_triangulation
```

To compute the value of a `CellField` at a position (x,y), the function getValue is implemented to overcome the bug encountered in Gridap.jl when using an unregular triangular mesh (see [issue 981](https://github.com/gridap/Gridap.jl/issues/981)). This function returns the value of the `CellField` at the positions given by the matrix of `Gridap.Point` p.
```@docs
    getValue
``` 



