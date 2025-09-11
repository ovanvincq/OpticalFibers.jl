# OpticalFibers - Using Unitful and Gridap

```@meta
CurrentModule = OpticalFibers
```

`OpticalFibers.jl` is based on two main packages:
- `Unitful.jl` that allows to add a unit to quantities
- `Gridap.jl` which is a finite element library

## Unitful
For convenience `OpticalFibers.jl` exports some functions of `Unitful.jl`: `uconvert`, `uconvertp`, `uconvertrp`, `ustrip`, `@u_str`, `unit`, `Units`, `NoUnits`, `dimension`.

Basically, you can create a quantity with a unit and do conversions:
```julia
julia> using OpticalFibers
julia> lambda=1.0u"µm"
1.0 μm
julia> lambda^2
1.0 μm^2
julia> uconvert(u"m",lambda)
1.0e-6 m
julia> ustrip(lambda)
1.0
julia> unit(lambda)
μm
julia> dimension(lambda)
𝐋
```

If you want to use others functionnalities, please install and load `Unitful.jl`.

## Gridap
For convenience `OpticalFibers.jl` exports some functions of `Gridap.jl` and `GridapGmsh.jl`: `integrate`, `Point`, `VectorValue`, `TensorValue`, `diagonal_tensor`, `CartesianDiscreteModel`, `GmshDiscreteModel`, `writevtk`, `CellField`,`cross`, `simplexify`, `Triangulation`, `Measure`.

With this package, you can create Point, Vectors, Tensor:
```julia
julia> point2D=Point(1.0,2.0)
VectorValue{2, Float64}(1.0, 2.0)
julia> vector=VectorValue(5.0,6.0,7.0)
VectorValue{3, Float64}(5.0, 6.0, 7.0)
julia> tensor=TensorValue(1.0,3.0,im,6.0)
TensorValue{2, 2, ComplexF64, 4}(1.0 + 0.0im, 3.0 + 0.0im, 0.0 + 1.0im, 6.0 + 0.0im)
julia> tensor2=diagonal_tensor(VectorValue(1.0,2.0,3.0))
TensorValue{3, 3, Float64, 9}(1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0)
julia> tensor2⋅vector
VectorValue{3, Float64}(5.0, 12.0, 21.0)
julia> cross(VectorValue(2.5,2.5,0),vector)
VectorValue{3, Float64}(17.5, -17.5, 2.5)
julia> cross(VectorValue(2.5,2.5,0),vector)
VectorValue{3, Float64}(17.5, -17.5, 2.5)
```

You can also create or load meshes:
```julia
julia> model1D = CartesianDiscreteModel((0,15),1500)
CartesianDiscreteModel()
julia> model2D = GmshDiscreteModel("./models/example1.msh")
Info    : Reading './models/example1.msh'...
Info    : 7 entities
Info    : 4733 nodes
Info    : 9564 elements
Info    : Done reading './models/example1.msh'
```

If you want to use others functionnalities, please install and load `Gridap.jl` and `GridapGmsh.jl`.

## UnitfulModel

In `OpticalFibers.jl`, a length unit can be added to a model:

```@docs
    UnitfulModel
```

For example:
```julia
julia> unitful_model1D=UnitfulModel(model1D,u"µm")
UnitfulModel(CartesianDiscreteModel(), μm)
julia> unitful_model2D=model2D*u"nm"
UnitfulModel(UnstructuredDiscreteModel(), nm)
```

## UnitfulField
`UnitfulField` is an abstract structure that describes a constant field that depends on space coordinates.

```@docs
    UnitfulField
    ScalarUnitfulField
    num_space_dims
    num_field_dims
```

A `UnitfulField` can be created with a matrix, a function or a Gridap field:
```@docs
    ArrayField
    FunctionField
    FEMField
```

Example of creation of a `FunctionField`:
```julia
julia>  f(x) = VectorValue(x[1]^2/1u"m^2", x[2]^2/1u"m^2")u"V"
f (generic function with 1 method)
julia> ff = FunctionField(2, f)
(::FunctionField{2, (2,), 𝐋^2 𝐌 𝐈^-1 𝐓^-3})     (generic function with 3 methods)
julia> f2(x) = VectorValue(x[1]/1u"m",x[2]/1u"m")u"A"
f2 (generic function with 1 method)
julia> ff2 = FunctionField(2, f2)
(::FunctionField{2, (2,), 𝐈})  (generic function with 3 methods)
```

`UnitfulField` is a `Function` so as the field can be evaluated directly:
```julia
julia> ff(Point(2.0,3.0)u"m")
VectorValue{2, Unitful.Quantity{Float64, 𝐋 ^2 𝐌  𝐈 ^-1 𝐓 ^-3, Unitful.FreeUnits{(V,), 𝐋 ^2 𝐌  𝐈 ^-1 𝐓 ^-3, nothing}}}(4.0 V, 9.0 V)
julia> ff([Point(2.0,3.0)u"m",Point(1.0,1.0)u"m"])
2-element Vector{VectorValue{2, Unitful.Quantity{Float64, 𝐋 ^2 𝐌  𝐈 ^-1 𝐓 ^-3, Unitful.FreeUnits{(V,), 𝐋 ^2 𝐌  𝐈 ^-1 𝐓 ^-3, nothing}}}}:
 VectorValue{2, Unitful.Quantity{Float64, 𝐋^2 𝐌 𝐈^-1 𝐓^-3, Unitful.FreeUnits{(V,), 𝐋^2 𝐌 𝐈^-1        𝐓^-3, nothing}}}(4.0 V, 9.0 V)
 VectorValue{2, Unitful.Quantity{Float64, 𝐋^2 𝐌 𝐈^-1 𝐓^-3, Unitful.FreeUnits{(V,), 𝐋^2 𝐌 𝐈^-1        𝐓^-3, nothing}}}(1.0 V, 1.0 V)
julia> ff([2.0u"m",5.0u"m"],[3.0u"m",4.0u"m"])
2×2 Matrix{VectorValue{2, Unitful.Quantity{Float64, 𝐋 ^2 𝐌  𝐈 ^-1 𝐓 ^-3, Unitful.FreeUnits{(V,), 𝐋 ^2 𝐌  𝐈 ^-1 𝐓 ^-3, nothing}}}}:
  VectorValue{2, Quantity{Float64, 𝐋^2 𝐌 𝐈^-1 𝐓^-3, FreeUnits{(V,), 𝐋^2 𝐌 𝐈^-1 𝐓^-3, nothing}}        }(4.0 V, 9.0 V)   VectorValue{2, Quantity{Float64, 𝐋^2 𝐌 𝐈^-1 𝐓^-3, FreeUnits{(V,), 𝐋^2 𝐌 𝐈^-1 𝐓^-3, nothing}}        }(4.0 V, 16.0 V)
 VectorValue{2, Quantity{Float64, 𝐋^2 𝐌 𝐈^-1 𝐓^-3, FreeUnits{(V,), 𝐋^2 𝐌 𝐈^-1 𝐓^-3, nothing}}        }(25.0 V, 9.0 V)  VectorValue{2, Quantity{Float64, 𝐋^2 𝐌 𝐈^-1 𝐓^-3, FreeUnits{(V,), 𝐋^2 𝐌 𝐈^-1 𝐓^-3, nothing}}        }(25.0 V, 16.0 V)
```

A lot of operator can be used with `UnitfulField` (+, *, cross, norm, real...):
```julia
julia> ffdiv=ff/2u"W"
(::FunctionField{2, (2,), 𝐈^-1})  (generic function with 3 methods)
julia> dot_ff = dot(ff, VectorValue(1.0, 1.0)u"V")
(::FunctionField{2, (), 𝐋^4 𝐌^2 𝐈^-2 𝐓^-6})     (generic function with 3 methods)
julia> dot_ff(Point(2.0,3.0)u"m")
13.0 V^2
julia> product_ff = ff⋅ff2
(::FunctionField{2, (), 𝐋^2 𝐌 𝐓^-3})    (generic function with 3 methods)
julia> uconvert(u"W",product_ff(Point(2.0,3.0)u"m"))
35.0 W
```

Integration is also implemented for `UnitfulField`
```julia
julia> integrate(ff,[-1,1]u"m",[-2,2]u"m")
VectorValue{2, Unitful.Quantity{Float64, 𝐋 ^4 𝐌  𝐈 ^-1 𝐓 ^-3, Unitful.FreeUnits{(m^2, V), 𝐋 ^4 𝐌  𝐈 ^-1 𝐓 ^-3, nothing}}}(-2.33333333333317 m^2 V, -2.3333333333328277 m^2 V)
```