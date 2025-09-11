using Gridap
using GridapGmsh
using Unitful
using LinearAlgebra

import Base: *, /, -, +, ^, real, imag, conj, sin, cos, tan, sqrt, abs, abs2, exp, log, log10
import LinearAlgebra: norm

export integrate, Point, VectorValue, TensorValue, diagonal_tensor, CartesianDiscreteModel, GmshDiscreteModel, writevtk, CellField, simplexify, Triangulation, Measure
export norm

@generated function _bc(f,a::NTuple{N},b::Union{Unitful.Units,Quantity}) where N
    s = "("
    for i in 1:N
        s *= "f(a[$i],b), "
    end
    s *= ")"
    Meta.parse(s)
end

  
@generated function _bc(f,b::Union{Unitful.Units,Quantity},a::NTuple{N}) where N
    s = "("
    for i in 1:N
      s *= "f(b,a[$i]), "
    end
    s *= ")"
    Meta.parse(s)
end

for op in (:*,:/)
    @eval begin
      function ($op)(a::Gridap.TensorValues.MultiValue,b::Unitful.Units)
        r = _bc($op,a.data,b)
        T = Gridap.TensorValues._eltype($op,r,a,b)
        M  = Gridap.TensorValues.change_eltype(a,T)
        M(r)
      end
    end
end

for op in (:*,:/)
    @eval begin
      function ($op)(a::Gridap.TensorValues.MultiValue,b::Quantity)
        r = _bc($op,a.data,b)
        T = Gridap.TensorValues._eltype($op,r,a,b)
        M  = Gridap.TensorValues.change_eltype(a,T)
        M(r)
      end
    end
end

for op in (:*,)
    @eval begin
      function ($op)(b::Unitful.Units,a::Gridap.TensorValues.MultiValue)
        r = _bc($op,b,a.data)
        T = Gridap.TensorValues._eltype($op,r,b,a)
        M  = Gridap.TensorValues.change_eltype(a,T)
        M(r)
      end
    end
end

for op in (:*,)
    @eval begin
      function ($op)(b::Quantity,a::Gridap.TensorValues.MultiValue)
        r = _bc($op,b,a.data)
        T = Gridap.TensorValues._eltype($op,r,b,a)
        M  = Gridap.TensorValues.change_eltype(a,T)
        M(r)
      end
    end
end

function LinearAlgebra.:norm(a::CellField)
  return Operation(norm)(a)
end

function Unitful.dimension(a::Gridap.TensorValues.MultiValue)
    return dimension(eltype(a))
end

function Unitful.unit(a::Gridap.TensorValues.MultiValue)
  return unit(eltype(a))
end

function Unitful.ustrip(u::Unitful.Units,a::Gridap.TensorValues.MultiValue)
  r=ustrip.(u,a.data)
  T = Gridap.TensorValues._eltype(ustrip,r,u,a)
  M  = Gridap.TensorValues.change_eltype(a,T)
  M(r)
end