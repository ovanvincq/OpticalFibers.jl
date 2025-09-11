# OpticalFibers.ModeSolvers - Modes and Fields

```@meta
CurrentModule = OpticalFibers.ModeSolvers
```

## Electromagnetic Fields
This module defines two types of fields: the scalar field (valid in the case of weakly guiding fibers) and the vector field. If the fiber has a cylindrical symmetry, it is sufficient to describe the variations of a scalar field along a radius (the 2D scalar field can be then obtained by multiplying the 1D field by $\cos(\ell \theta)$ or $\sin(\ell \theta)$ where $\ell$ is the azimuthal number).

```@docs
FiberEMField
ScalarFiberEMField
ScalarFiberEMField1D
ScalarFiberEMField2D
VectorFiberEMField
```

The fields can be multiplied by a scalar and conjugated. It is also possible to take the real part or the imaginary part of a field.
Two fields of the same type can also be added or subtracted.


## Modes
A mode is a structure that contains a name, an effective index, the wavelength and, optionally, an electromagnetic field.

This package assumes that the field is proportionnal to $\exp\left(i\beta z-\omega t\right)$.
```@docs
Mode
```

## Modes/Fields conversion
It is possible to convert a scalar 1D mode or field to a scalar 2d one. The argument `orientation_angle` corresponds to the orientation in degrees (0 correspond to $\cos(\ell \theta)$ and 90 to $\sin(\ell \theta)$).
```@docs
convertTo2D
```

A mode that contains a scalar field can be converted to a vector one by using the approximation: $\vec{k}=\beta \vec{e}_z$ and $E_z=H_z=0$. The argument `polarization_angle` corresponds to the direction of polarization of the vector mode (0 correspond to $x$ and 90 to $y$)
```@docs
convertToVector
```

## Propagation of modes
The field due to the propagation of a mode after a propagation distance z can be calculated. This function simply returns the EMField multiplied by $\exp(i\beta z)$.
```@docs
getEMField
```

## Mode normalization
The modes given by the mode solvers of this package are not normalized. You can normalize them by using the function `normalize!`.
The vector modes are normalized so as $\vert \iint \frac{1}{2}\left(\vec{E}\wedge\vec{H}^*\right) dS \vert=1$ and the scalar mode are normalized like vector modes if the optional parameter `unitIntegral` is `false` or so as $ \iint \vert E \vert^2 dS=1 $ if `unitIntegral` is `true`.
```@docs
normalize
```

## Overlap integral
The overlap integral between two modes or between an EM field and a mode f$_1$ and f$_2$ is $\langle f_1 \vert f_2 \rangle = \iint E_1 E_2^* dS$ for scalar fields/modes and $\langle f_1 \vert f_2 \rangle = \iint \frac{1}{2}\left(\vec{E}_1 \wedge \vec{H}_2^*\right)_z dS$ for vector fields/modes.  

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
In this function, n2 is the nonlinear index and n0 is the refractive index involved in the relation between $n_2$ and the third order susceptibility $\chi^{(3)}$. They must be constants or functions of the tuple (r,) in the 1D case or (x,y) in the 2D case. If n0=0, this function assumes that n0≃neff (which is correct for weakly-guiding fibers).
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
If the field or mode is computed with a FEM solver, the function writevtk of the package Gridap is overloaded to save the EM field or the mode to a file that can be opened with ParaView.
```@docs
    writevtk
```





