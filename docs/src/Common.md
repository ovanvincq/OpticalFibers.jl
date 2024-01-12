# OpticalFibers - Common structs and functions

```@meta
CurrentModule = OpticalFibers
```

## Convenient functions
```@docs
    piecewiseIndex
    meshgrid
    derivative
    ring
    approx_nFSM_PCF
    approx_neff_PCF
    add_cylindrical_PML
```

## Tensor
```@docs
    tensor3
    tensor3(::Union{Function,Number},::Union{Function,Number},::Union{Function,Number},::Union{Function,Number},::Union{Function,Number},::Union{Function,Number},::Union{Function,Number},::Union{Function,Number},::Union{Function,Number})
    tensor3(::Union{Function,Number},::Union{Function,Number},::Union{Function,Number})
    tensor3(::Union{Function,Number})
    inverse
```


