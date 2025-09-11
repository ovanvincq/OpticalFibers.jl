using Gridap
using Gridap.ReferenceFEs
using Gridap.CellData
using Gridap.Arrays
using Unitful
using StaticArrays
using NearestNeighbors
using LinearAlgebra
using Interpolations
using Trapz

export uconvert, uconvertp, uconvertrp, ustrip, @u_str, unit, Units, NoUnits, dimension, cross,×

import LinearAlgebra: dot,⋅
import Gridap: cross,×
import Base: zero

export dot,⋅

for op in (:sin,:cos,:tan,:exp,:log,:sqrt,:log10)
    @eval begin
      ($op)(a::CellField) = Operation($op)(a)
    end
end

export UnitfulField, ScalarUnitfulField, ArrayField, FunctionField, FEMField, num_space_dims, num_field_dims, UnitfulModel, get_triangulation

"""
    struct UnitfulModel

- model :: `DiscreteModel`
- unit :: `Unitful.Units{<:Any,Unitful.dimension(u"m")}`
"""
struct UnitfulModel
    model::DiscreteModel
    unit::Unitful.Units{<:Any,Unitful.dimension(u"m")}
end

Base.:*(m::DiscreteModel, u::Unitful.Units{<:Any,Unitful.dimension(u"m")}) = UnitfulModel(m,u)
Base.:*(u::Unitful.Units{<:Any,Unitful.dimension(u"m")},m::DiscreteModel) = UnitfulModel(m,u)

function Unitful.unit(m::UnitfulModel)
    return m.unit
end

"""
    abstract type UnitfulField{Dp,Df,dim} <:Function where dim<:Unitful.Dimensions end 

- Dp: dimension of the space: 1, 2, 3...
- Df: dimension of the field: () for a scalar field, (2,) for 2D vector or (3,3) for a 3x3 tensor
- dim: unit dimension of the field
"""
abstract type UnitfulField{Dp,Df,dim} <:Function where dim<:Unitful.Dimensions end 

"""
    const ScalarUnitfulField{Dp,dim}=UnitfulField{Dp,(),dim} 
"""
const ScalarUnitfulField{Dp,dim}=UnitfulField{Dp,(),dim}

"""
    function num_space_dims(::UnitfulField{Dp,Df,dim}) where {Dp,Df,dim}

returns Dp
"""
function num_space_dims(::UnitfulField{Dp,Df,dim}) where {Dp,Df,dim}
    return Dp 
end

"""
    function num_field_dims(::UnitfulField{Dp,Df,dim}) where {Dp,Df,dim}

returns Df
"""
function num_field_dims(::UnitfulField{Dp,Df,dim}) where {Dp,Df,dim}
    return Df 
end

function Unitful.dimension(::UnitfulField{Dp,Df,dim}) where {Dp,Df,dim}
    return dim
end

#Base.:transpose(f::UnitfulField) = f
Broadcast.:broadcastable(f::UnitfulField)=Ref(f)

(f::UnitfulField{Dp,<:Any,<:Any})(x::Vararg{Union{realLength,AbstractVector{<:realLength}}}) where Dp =
begin
    if length(x)!=Dp
        throw(ArgumentError("The number of arguments must be equal to the number of space dimensions"))
    end
    t=Vector{Array}(undef,0)
    scalar=all(isa.(x,realLength))
    if (!scalar)
        for k=1:length(x)
            if isa(x[k],realLength)
                d=ones(Int64,length(x))
                push!(t,reshape([x[k]],Tuple(d)...))
            else
                d=ones(Int64,length(x))
                d[k]=length(x[k])
                push!(t,reshape(x[k],Tuple(d)...))
            end
        end
        return f(Point.(t...))
    else
        return f(Point([k[1][1] for k in x]))
    end
end 

######################   ArrayField   ######################
#Field defined by an array on a rectangular grid
"""
    struct ArrayField{Dp,Df,dim} <:UnitfulField{Dp,Df,dim}

Field defined by an array on a rectangular grid
- position: `AbstractVector{AbstractVector{<:realLength}}`
- value: `AbstractArray{<:Number}`
"""
struct ArrayField{Dp,Df,dim} <:UnitfulField{Dp,Df,dim}
    position::AbstractVector{AbstractVector{<:realLength}}
    value::AbstractArray{<:Number}
    function ArrayField(position::Union{AbstractVector{<:AbstractVector{<:realLength}},AbstractVector{<:realLength}},value::AbstractArray{<:Number})
        if isa(position,AbstractVector{<:realLength})
            position=[position]
        end
        Dp=length(position)
        x=Tuple(length(position[j]) for j=1:Dp)
        if size(value)!=x
            throw(ArgumentError("Dimensions do not match"))
        end
        if length(unique(dimension.(value)))!=1
            throw(ArgumentError("All values must have the same dimension"))
        end
        Df=size(value[1])
        dim=dimension(value[1])
        new{Dp,Df,dim}(deepcopy(position),deepcopy(value))
    end
end

(f::ArrayField{Dp,<:Any,<:Any})(r::VectorValue{Dp,<:realLength}) where Dp = 
begin 
    interp=linear_interpolation(Tuple(f.position),f.value,extrapolation_bc=Tuple(0*f.value[1]))
    result=interp(r...)
    if isa(result,Tuple)
        return zero(f.value[1])
    else
        return result
    end
end
(f::ArrayField{Dp,<:Any,<:Any})(r::AbstractArray{<:VectorValue{Dp,<:realLength}}) where Dp = 
begin
    interp=linear_interpolation(Tuple(f.position),f.value,extrapolation_bc=Tuple(0*f.value[1]))
    result=map(x->interp(x...),r)
    result[isa.(result,Tuple)].=zero(f.value[1])
    return result
end

function Base.zero(f::ArrayField)
    return f*0
end

function concatenate_positions(f1::ArrayField{Dp,<:Any,<:Any},f2::ArrayField{Dp,<:Any,<:Any}) where {Dp}
    pos=sort.(unique.(vcat.(f1.position,f2.position)))
    n=length.(pos)
    t=Vector{Array}(undef,0)
    for k in axes(pos,1)
        d=ones(Int64,length(n))
        d[k]=n[k]
        push!(t,reshape(pos[k],Tuple(d)...))
    end
    v=VectorValue.(t...)
    return pos,v
end

Base.:*(f::ArrayField,k::Number) = ArrayField(deepcopy(f.position),k.*f.value);
Base.:*(k::Number,f::ArrayField) = ArrayField(deepcopy(f.position),k.*f.value);
Base.:/(f::ArrayField,k::Number) = ArrayField(deepcopy(f.position),f.value./k);
Base.:/(k::Number,f::ArrayField) = ArrayField(deepcopy(f.position),k./f.value);
Base.:+(f::ArrayField{<:Any,Df,<:Any},k::Number) where Df = (size(k)==Df) ? ArrayField(deepcopy(f.position),f.value.+k) : throw(ArgumentError("The arguments must have the same size"));
Base.:+(k::Number,f::ArrayField{<:Any,Df,<:Any}) where Df = (size(k)==Df) ? ArrayField(deepcopy(f.position),k.+f.value) : throw(ArgumentError("The arguments must have the same size"));
Base.:-(f::ArrayField{<:Any,Df,<:Any},k::Number) where Df = (size(k)==Df) ? ArrayField(deepcopy(f.position),f.value.-k) : throw(ArgumentError("The arguments must have the same size"));
Base.:-(k::Number,f::ArrayField{<:Any,Df,<:Any}) where Df = (size(k)==Df) ? ArrayField(deepcopy(f.position),k.-f.value) : throw(ArgumentError("The arguments must have the same size"));

Base.:^(f::ArrayField{<:Any,(),<:Any},k::Number) = ArrayField(deepcopy(f.position),(f.value).^k);
Base.:inv(f::ArrayField{<:Any,<:Any,<:Any}) = ArrayField(deepcopy(f.position),inv.(f.value));

LinearAlgebra.:dot(f::ArrayField,k::Number) = ArrayField(deepcopy(f.position),dot.(f.value,Ref(k)));
LinearAlgebra.:dot(k::Number,f::ArrayField) = ArrayField(deepcopy(f.position),dot.(Ref(k),f.value));
LinearAlgebra.:cross(f::ArrayField,k::Number) = ArrayField(deepcopy(f.position),cross.(f.value,Ref(k)));
LinearAlgebra.:cross(k::Number,f::ArrayField) = ArrayField(deepcopy(f.position),cross.(Ref(k),f.value));

for op in (:+,:-,:*,:/,:cross,:dot)
    @eval begin
        function ($op)(f1::ArrayField,f2::ArrayField)
            if (f1.position==f2.position)
                ArrayField(deepcopy(f1.position),($op).(f1.value,f2.value));
            else
                pos,v=concatenate_positions(f1,f2)
                ArrayField(pos,($op).(f1(v),f2(v)))
            end
        end
    end
    @eval begin
        function ($op)(f1::ArrayField,f2::Function)
            t=Vector{Array}(undef,0)
            for k=1:length(f1.position)
                d=ones(Int64,length(f1.position))
                d[k]=length(f1.position[k])
                push!(t,reshape(f1.position[k],Tuple(d)...))
            end
            ArrayField(f1.position,($op).(f1.value,f2.(Point.(t...))))
        end
    end
    @eval begin
        function ($op)(f2::Function,f1::ArrayField)
            t=Vector{Array}(undef,0)
            for k=1:length(f1.position)
                d=ones(Int64,length(f1.position))
                d[k]=length(f1.position[k])
                push!(t,reshape(f1.position[k],Tuple(d)...))
            end
            ArrayField(f1.position,($op).(f2.(Point.(t...)),f1.value))
        end
    end
end

LinearAlgebra.:norm(f::ArrayField) = ArrayField(deepcopy(f.position),norm.(f.value));
for op in (:+,:-,:sqrt,:real,:imag,:conj,:sin,:cos,:tan,:abs,:abs2,:exp,:log,:log10)
    @eval begin
        function ($op)(f::ArrayField)
            ArrayField(deepcopy(f.position),($op).(f.value));
        end
    end 
end

function Gridap.integrate(f::ArrayField{<:Any,Df,<:Any};kwargs...) where Df
    if Df==()
        return trapz(Tuple(f.position),f.value);
    else
        N=length(f.value[1].data)
        result=Vector(undef,0)
        for i=1:N
            push!(result,trapz(Tuple(f.position),getindex.(getfield.(f.value,:data),i)));
        end
        T = eltype(result[1])
        M  = Gridap.TensorValues.change_eltype(f.value[1],T)
        return M(Tuple(result))
    end
end

function Gridap.integrate(f::ArrayField{Dp,Df,<:Any},a::AbstractVector{<:realLength},b::AbstractVector{<:realLength};kwargs...) where {Dp,Df}
    if ((length(a)!=length(b)) || (length(a)!=Dp))
        throw(ArgumentError("Lengths of a and b must be equal to Dp"))
    end
    pos=Vector{Vector{Int64}}(undef,length(a))
    signe=1
    for k=1:length(a)
        if a[k]<b[k]
            pos[k]=findall(f.position[k].>=a[k] .&& f.position[k].<=b[k])
        else
            pos[k]=findall(f.position[k].<=a[k] .&& f.position[k].>=b[k])
            signe=-1*signe
        end
    end
    position2=[f.position[k][pos[k]] for k=1:length(f.position)]
    f2=ArrayField(position2,f.value[pos...])
    return signe*integrate(f2)
end


######################   FunctionField   ######################

#the function must return a Number or a VectorValue (Gridap)
"""
    struct FunctionField{Dp,Df,dim} <:UnitfulField{Dp,Df,dim}

- value: `Function`

The constructor of a FunctionField requires the value of Dp : `FunctionField(Dp::Int,value::Function)`
"""
struct FunctionField{Dp,Df,dim} <:UnitfulField{Dp,Df,dim}
    value::Function
    function FunctionField(Dp::Int,value::Function)
        x=VectorValue(Tuple(0.0u"m" for j=1:Dp))
        try
            value(x)
        catch e
            #println(x)
            throw(ArgumentError("The function is not valid"))
        end
        Df=size(value(x))
        dim=dimension(value(x))
        new{Dp,Df,dim}(x->deepcopy(value)(x))
    end
    function FunctionField(f::ArrayField{Dp,Df,dim}) where {Dp,Df,dim}
        new{Dp,Df,dim}(f)
    end
end

function ArrayField(position::Union{AbstractVector{<:AbstractVector{<:realLength}},AbstractVector{<:realLength}},f::FunctionField)
    n=length.(position)
    t=Vector{Array}(undef,0)
    for k in axes(position,1)
        d=ones(Int64,length(n))
        d[k]=n[k]
        push!(t,reshape(position[k],Tuple(d)...))
    end
    v=f(Point.(t...))
    return ArrayField(position,v)
end



(f::FunctionField{Dp,<:Any,<:Any})(r::VectorValue{Dp,<:realLength}) where Dp = f.value(r)
(f::FunctionField{Dp,<:Any,<:Any})(r::AbstractArray{<:VectorValue{Dp,<:realLength}}) where Dp=f.value.(r)

function Base.zero(f::FunctionField{Dp,<:Any,<:Any}) where{Dp}
    x=VectorValue(Tuple(0u"m" for j=1:Dp))
    zz=zero(f(x))
    return FunctionField(Dp,x->zz)
end

Base.:*(f::FunctionField{Dp,<:Any,<:Any},k::Number) where Dp = FunctionField(Dp,x->k*deepcopy(f.value)(x));
Base.:*(k::Number,f::FunctionField{Dp,<:Any,<:Any}) where Dp = FunctionField(Dp,x->k*deepcopy(f.value)(x));
Base.:/(f::FunctionField{Dp,<:Any,<:Any},k::Number) where Dp = FunctionField(Dp,x->deepcopy(f.value)(x)/k);
Base.:/(k::Number,f::FunctionField{Dp,<:Any,<:Any}) where Dp = FunctionField(Dp,x->k/deepcopy(f.value)(x));
Base.:+(f::FunctionField{Dp,Df,<:Any},k::Number) where {Dp,Df}  = (size(k)==Df) ? FunctionField(Dp,x->deepcopy(f.value)(x)+k) : throw(ArgumentError("The arguments must have the same size"));
Base.:+(k::Number,f::FunctionField{Dp,Df,<:Any}) where {Dp,Df}  = (size(k)==Df) ? FunctionField(Dp,x->k+deepcopy(f.value)(x)) : throw(ArgumentError("The arguments must have the same size"));
Base.:-(f::FunctionField{Dp,Df,<:Any},k::Number) where {Dp,Df}  = (size(k)==Df) ? FunctionField(Dp,x->deepcopy(f.value)(x)-k) : throw(ArgumentError("The arguments must have the same size"));
Base.:-(k::Number,f::FunctionField{Dp,Df,<:Any}) where {Dp,Df}  = (size(k)==Df) ? FunctionField(Dp,x->k-deepcopy(f.value)(x)) : throw(ArgumentError("The arguments must have the same size"));

Base.:^(f::FunctionField{Dp,(),<:Any},k::Number) where Dp = FunctionField(Dp,x->(deepcopy(f.value)(x))^k);
Base.:inv(f::FunctionField{Dp,<:Any,<:Any}) where Dp = FunctionField(Dp,x->inv(deepcopy(f.value)(x)));

LinearAlgebra.:dot(f1::FunctionField{Dp,<:Any,<:Any},k::Number) where Dp = FunctionField(Dp,x->dot(deepcopy(f1.value)(x),k));
LinearAlgebra.:dot(k::Number,f1::FunctionField{Dp,<:Any,<:Any}) where Dp = FunctionField(Dp,x->dot(k,deepcopy(f1.value)(x)));
Gridap.:cross(f1::FunctionField{Dp,<:Any,<:Any},k::Number) where Dp = FunctionField(Dp,x->cross(deepcopy(f1.value)(x),k));
Gridap.:cross(k::Number,f1::FunctionField{Dp,<:Any,<:Any}) where Dp = FunctionField(Dp,x->cross(k,deepcopy(f1.value)(x)));

LinearAlgebra.:norm(f::FunctionField{Dp,<:Any,<:Any}) where Dp = FunctionField(Dp,x->LinearAlgebra.norm(deepcopy(f.value)(x)));

for op in (:+,:-,:*,:/,:cross,:dot)
    @eval begin
        function ($op)(f1::FunctionField{Dp,<:Any,<:Any},f2::FunctionField{Dp,<:Any,<:Any}) where Dp
            FunctionField(Dp,x->($op)(deepcopy(f1.value)(x),deepcopy(f2.value)(x)));
        end
    end
    @eval begin
        function ($op)(f1::FunctionField{Dp,<:Any,<:Any},f2::ArrayField{Dp,<:Any,<:Any}) where Dp
            ($op)(f1,FunctionField(Dp,f2))
        end
    end
    @eval begin
        function ($op)(f1::ArrayField{Dp,<:Any,<:Any},f2::FunctionField{Dp,<:Any,<:Any}) where Dp
            ($op)(FunctionField(Dp,f1),f2)
        end
    end
    @eval begin
        function ($op)(f1::FunctionField{Dp,<:Any,<:Any},f2::Function) where Dp
            FunctionField(Dp,x->($op)(deepcopy(f1.value)(x),deepcopy(f2)(x)));
        end
    end
    @eval begin
        function ($op)(f1::Function,f2::FunctionField{Dp,<:Any,<:Any}) where Dp
            FunctionField(Dp,x->($op)(deepcopy(f1)(x),deepcopy(f2.value)(x)));
        end
    end
end

for op in (:+,:-,:sqrt,:real,:imag,:conj,:sin,:cos,:tan,:abs,:abs2,:exp,:log,:log10)
    @eval begin
        function ($op)(f::FunctionField{Dp,<:Any,<:Any}) where Dp
            FunctionField(Dp,x->($op)(deepcopy(f.value)(x)));
        end
    end 
end

function Gridap.integrate(f::FunctionField{Dp,Df,<:Any};characteristic_length::AbstractVector{<:realLength}=ones(Dp)u"m",kwargs...) where {Dp,Df}
    return function_integrate_unitful(f.value,-Inf*ones(Dp)*u"m",Inf*ones(Dp)*u"m";characteristic_length=characteristic_length,kwargs...)
end

function Gridap.integrate(f::FunctionField{Dp,Df,<:Any},a::AbstractVector{<:realLength},b::AbstractVector{<:realLength};characteristic_length::AbstractVector{<:realLength}=ones(Dp)u"m",kwargs...) where {Dp,Df}
    if ((length(a)!=length(b)) || (length(a)!=Dp))
        throw(ArgumentError("Lengths of a and b must be equal to Dp"))
    end
    return function_integrate_unitful(f.value,a,b;characteristic_length=characteristic_length,kwargs...)
end

######################   FEMField   ######################
#FEMField is a field defined on a triangulation with a measure
"""
    struct FEMField{Dp,Df,dim} <:UnitfulField{Dp,Df,dim}

- value: `CellField`
- dΩ: `Gridap.CellData.Measure`
- PositionUnit: `Unitful.Units{<:Any,Unitful.dimension(u"m")}`
- FieldUnit: `Unitful.Units`
"""
struct FEMField{Dp,Df,dim} <:UnitfulField{Dp,Df,dim}
    value::CellField
    dΩ::Gridap.CellData.Measure #For integration
    PositionUnit::Unitful.Units{<:Any,Unitful.dimension(u"m")}
    FieldUnit::Unitful.Units
    function FEMField(value::CellField,dΩ::Gridap.CellData.Measure,pu,fu)
        trian=get_triangulation(value)
        trian2=get_triangulation(dΩ.quad)
        if trian!=trian2
            throw(ArgumentError("The Field and the measure must have the same Triangulation"))
        end
        Dp=num_point_dims(trian)
        Df=size(getValue(value,get_node_coordinates(trian)[1]))
        new{Dp,Df,dimension(fu)}(value,dΩ,pu,fu)
    end
    function FEMField(f::FunctionField{Dp,Df,dim},dΩ::Gridap.CellData.Measure,pu,fu) where {Dp,Df,dim}
        trian=get_triangulation(dΩ.quad)
        if dimension(fu)!=dim
            throw(ArgumentError("Dimension error"))
        end
        cf=CellField(x->ustrip(fu,f.value(x*pu)),trian)
        new{Dp,Df,dim}(cf,dΩ,pu,fu)
    end
    function FEMField(f::ArrayField{Dp,Df,dim},dΩ::Gridap.CellData.Measure,pu,fu) where {Dp,Df,dim}
        trian=get_triangulation(dΩ.quad)
        if dimension(fu)!=dim
            throw(ArgumentError("Dimension error"))
        end
        cf=CellField(x->ustrip(fu,f(x*pu)),trian)
        new{Dp,Df,dim}(cf,dΩ,pu,fu)
    end
end

for op in (:+,:-,:*,:/,:cross,:dot)
    @eval begin
        function ($op)(f1::FunctionField{Dp,<:Any,<:Any},f2::FEMField{Dp,<:Any,<:Any}) where Dp
            ($op)(f1,FunctionField(Dp,f2))
        end
    end
    @eval begin
        function ($op)(f1::FEMField{Dp,<:Any,<:Any},f2::FunctionField{Dp,<:Any,<:Any}) where Dp
            ($op)(FunctionField(Dp,f1),f2)
        end
    end
end

Gridap.get_triangulation(f::FEMField) = get_triangulation(f.value)

function FunctionField(f::FEMField{Dp,<:Any,<:Any}) where Dp
    return FunctionField(Dp,f)
end
function ArrayField(position::Union{AbstractVector{<:AbstractVector{<:realLength}},AbstractVector{<:realLength}},f::FEMField)
    n=length.(position)
    t=Vector{Array}(undef,0)
    for k in axes(position,1)
        d=ones(Int64,length(n))
        d[k]=n[k]
        push!(t,reshape(position[k],Tuple(d)...))
    end
    v=f(Point.(t...))
    return ArrayField(position,v)
end

(f::FEMField{Dp,<:Any,<:Any})(r::VectorValue{Dp,<:realLength}) where Dp=getValue(f.value,ustrip(f.PositionUnit,r))*f.FieldUnit
(f::FEMField{Dp,<:Any,<:Any})(r::AbstractArray{<:VectorValue{Dp,<:realLength}}) where Dp=getValue(f.value,ustrip.(f.PositionUnit,r)).*f.FieldUnit

function Base.zero(f::FEMField{<:Any,Df,<:Any}) where{Df}
    trian=get_triangulation(f.value)
    zz=zero(getValue(f.value,get_node_coordinates(trian)[1]))
    return FEMField(CellField(x->zz,trian),f.dΩ,f.PositionUnit,f.FieldUnit)    
end

Base.:*(f::FEMField,k::Number) = FEMField(ustrip(k)*f.value,f.dΩ,f.PositionUnit,f.FieldUnit*unit(k));
Base.:*(k::Number,f::FEMField) = FEMField(ustrip(k)*f.value,f.dΩ,f.PositionUnit,f.FieldUnit*unit(k));
Base.:/(f::FEMField,k::Number) = FEMField(f.value/ustrip(k),f.dΩ,f.PositionUnit,f.FieldUnit/unit(k));
Base.:/(k::Number,f::FEMField) = FEMField(ustrip(k)/f.value,f.dΩ,f.PositionUnit,f.FieldUnit/unit(k));
Base.:+(f::FEMField{<:Any,Df,<:Any},k::Number) where Df  = (size(k)==Df) ? FEMField(f.value+ustrip(f.FieldUnit,k),f.dΩ,f.PositionUnit,f.FieldUnit) : throw(ArgumentError("The arguments must have the same size"));
Base.:+(k::Number,f::FEMField{<:Any,Df,<:Any}) where Df  = (size(k)==Df) ? FEMField(ustrip(f.FieldUnit,k)+f.value,f.dΩ,f.PositionUnit,f.FieldUnit) : throw(ArgumentError("The arguments must have the same size"));
Base.:-(f::FEMField{<:Any,Df,<:Any},k::Number) where Df  = (size(k)==Df) ? FEMField(f.value-ustrip(f.FieldUnit,k),f.dΩ,f.PositionUnit,f.FieldUnit) : throw(ArgumentError("The arguments must have the same size"));
Base.:-(k::Number,f::FEMField{<:Any,Df,<:Any}) where Df  = (size(k)==Df) ? FEMField(ustrip(f.FieldUnit,k)-f.value,f.dΩ,f.PositionUnit,f.FieldUnit) : throw(ArgumentError("The arguments must have the same size"));

Base.:^(f::FEMField{<:Any,(),<:Any},k::Int64) = 
begin
    if (k==0)
        trian=get_triangulation(f.value)
        return FEMField(CellField(x->1.0,trian),f.dΩ,f.PositionUnit,NoUnits)
    elseif (k<0)
        val=1.0/f.value
        for j=2:(-k)
            val=val/f.value
        end
        return FEMField(val,f.dΩ,f.PositionUnit,f.FieldUnit^(-k));
    else
        val=f.value
        for j=2:k
            val=val*f.value
        end
        return FEMField(val,f.dΩ,f.PositionUnit,f.FieldUnit^k);
    end
end
Base.:inv(f::FEMField{<:Any,<:Any,<:Any}) =
begin
    FEMField(inv(f.value),f.dΩ,f.PositionUnit,f.FieldUnit^(-1));
end

LinearAlgebra.:dot(f1::FEMField,k::Number) = FEMField(dot(f1.value,ustrip(k)),f1.dΩ,f1.PositionUnit,f1.FieldUnit*unit(k));
LinearAlgebra.:dot(k::Number,f1::FEMField) = FEMField(dot(ustrip(k),f1.value),f1.dΩ,f1.PositionUnit,f1.FieldUnit*unit(k));
LinearAlgebra.:cross(f1::FEMField,k::Number) = FEMField(cross(f1.value,ustrip(k)),f1.dΩ,f1.PositionUnit,f1.FieldUnit*unit(k));
LinearAlgebra.:cross(k::Number,f1::FEMField) = FEMField(cross(ustrip(k),f1.value),f1.dΩ,f1.PositionUnit,f1.FieldUnit*unit(k));

Base.:-(f1::FEMField{Dp,Df,dim},f2::FEMField{Dp,Df,dim}) where{Dp,Df,dim} = ((f1.dΩ==f2.dΩ) && (f1.PositionUnit==f2.PositionUnit)) ? FEMField(f1.value-f2.value*ustrip(NoUnits,f2.FieldUnit/(1*f1.FieldUnit)),f1.dΩ,f1.PositionUnit,f1.FieldUnit) : throw(ArgumentError("dΩ and PositionUnit must be the same for both arguments"));
Base.:+(f1::FEMField{Dp,Df,dim},f2::FEMField{Dp,Df,dim}) where{Dp,Df,dim} = ((f1.dΩ==f2.dΩ) && (f1.PositionUnit==f2.PositionUnit)) ? FEMField(f1.value+f2.value*ustrip(NoUnits,f2.FieldUnit/(1*f1.FieldUnit)),f1.dΩ,f1.PositionUnit,f1.FieldUnit) : throw(ArgumentError("dΩ and PositionUnit must be the same for both arguments"));
Base.:/(f1::FEMField,f2::FEMField) = ((f1.dΩ==f2.dΩ) && (f1.PositionUnit==f2.PositionUnit)) ? FEMField(f1.value/f2.value,f1.dΩ,f1.PositionUnit,f1.FieldUnit/f2.FieldUnit) : throw(ArgumentError("dΩ and PositionUnit must be the same for both arguments"));
Base.:*(f1::FEMField,f2::FEMField) = ((f1.dΩ==f2.dΩ) && (f1.PositionUnit==f2.PositionUnit)) ? FEMField(f1.value*f2.value,f1.dΩ,f1.PositionUnit,f1.FieldUnit*f2.FieldUnit) : throw(ArgumentError("dΩ and PositionUnit must be the same for both arguments"));
LinearAlgebra.:dot(f1::FEMField,f2::FEMField) = ((f1.dΩ==f2.dΩ) && (f1.PositionUnit==f2.PositionUnit)) ? FEMField(dot(f1.value,f2.value),f1.dΩ,f1.PositionUnit,f1.FieldUnit*f2.FieldUnit) : throw(ArgumentError("dΩ and PositionUnit must be the same for both arguments"));
Gridap.:cross(f1::FEMField,f2::FEMField) = ((f1.dΩ==f2.dΩ) && (f1.PositionUnit==f2.PositionUnit)) ? FEMField(cross(f1.value,f2.value),f1.dΩ,f1.PositionUnit,f1.FieldUnit*f2.FieldUnit) : throw(ArgumentError("dΩ and PositionUnit must be the same for both arguments"));
LinearAlgebra.:norm(f::FEMField) = FEMField(norm(f.value),f.dΩ,f.PositionUnit,f.FieldUnit)

for op in (:*,:cross,:dot)
    @eval begin
        function ($op)(f::FEMField{Dp,<:Any,<:Any},fun::Function) where Dp
            x0=VectorValue(Tuple(0*f.PositionUnit for j=1:Dp))
            fun_unit=unit(fun(x0));
            fun2(x)=ustrip(fun_unit,fun(x*f.PositionUnit))
            return FEMField(($op)(f.value,fun2),f.dΩ,f.PositionUnit,f.FieldUnit*fun_unit)
        end
    end
    @eval begin
        function ($op)(fun::Function,f::FEMField{Dp,<:Any,<:Any}) where Dp
            x0=VectorValue(Tuple(0*f.PositionUnit for j=1:Dp))
            fun_unit=unit(fun(x0));
            fun2(x)=ustrip(fun_unit,fun(x*f.PositionUnit))
            return FEMField(($op)(fun2,f.value),f.dΩ,f.PositionUnit,f.FieldUnit*fun_unit)
        end
    end
end

for op in (:+,:-)
    @eval begin
        function ($op)(f::FEMField{Dp,<:Any,<:Any},fun::Function) where Dp
            x0=VectorValue(Tuple(0*f.PositionUnit for j=1:Dp))
            fun_unit=unit(fun(x0));
            if dimension(fun_unit)!=dimension(f.FieldUnit)
                throw(ArgumentError("The function must return a value with the same dimension as the FieldUnit"))
            end
            fun2=x->ustrip(f.FieldUnit,fun(x*f.PositionUnit))
            return FEMField(($op)(f.value,fun2),f.dΩ,f.PositionUnit,f.FieldUnit)
        end
    end
    @eval begin
        function ($op)(fun::Function,f::FEMField{Dp,<:Any,<:Any}) where Dp
            x0=VectorValue(Tuple(0*f.PositionUnit for j=1:Dp))
            fun_unit=unit(fun(x0));
            if dimension(fun_unit)!=dimension(f.FieldUnit)
                throw(ArgumentError("The function must return a value with the same dimension as the FieldUnit"))
            end
            fun2=x->ustrip(f.FieldUnit,fun(x*f.PositionUnit))
            return FEMField(($op)(fun2,f.value),f.dΩ,f.PositionUnit,f.FieldUnit)
        end
    end
end

Base.:/(f::FEMField,fun::Function)=
begin
    Ω=get_triangulation(f.value)
    Dp=num_point_dims(Ω)
    x0=VectorValue(Tuple(0*f.PositionUnit for j=1:Dp))
    fun_unit=unit(fun(x0));
    fun2=x->ustrip(fun_unit,fun(x*f.PositionUnit))
    return FEMField(f.value/fun2,f.dΩ,f.PositionUnit,f.FieldUnit/fun_unit)
end
Base.:/(fun::Function,f::FEMField)=
begin
    Ω=get_triangulation(f.value)
    Dp=num_point_dims(Ω)
    x0=VectorValue(Tuple(0*f.PositionUnit for j=1:Dp))
    fun_unit=unit(fun(x0));
    fun2=x->ustrip(fun_unit,fun(x*f.PositionUnit))
    return FEMField(fun2/f.value,f.dΩ,f.PositionUnit,fun_unit/f.FieldUnit)
end


for op in (:+,:-,:sqrt,:real,:imag,:conj,:sin,:cos,:tan,:abs,:abs2,:exp,:log,:log10)
    @eval begin
        function ($op)(f::FEMField)
            FEMField(($op)(f.value),f.dΩ,f.PositionUnit,unit(($op)(1*f.FieldUnit)));
        end
    end 
end

function Gridap.integrate(f::FEMField{Dp,<:Any,<:Any};kwargs...) where {Dp}
    return sum(integrate(f.value,f.dΩ))*f.FieldUnit*f.PositionUnit^Dp;
end

function Gridap.integrate(f::FEMField{Dp,<:Any,<:Any},a::AbstractVector{<:realLength},b::AbstractVector{<:realLength};kwargs...) where {Dp}
    if ((length(a)!=length(b)) || (length(a)!=Dp))
        throw(ArgumentError("Lengths of a and b must be equal to Dp"))
    end
    aa=ustrip.(f.PositionUnit,a)
    bb=ustrip.(f.PositionUnit,b)
    Ω=get_triangulation(f.value)
    ff=f.value
    for k=1:length(a)
        if (aa[k]<bb[k])
            ff=ff*CellField(x->(x[k]<=bb[k])*(x[k]>=aa[k])*1,Ω)
        else
            ff=-ff*CellField(x->(x[k]<=aa[k])*(x[k]>=bb[k])*1,Ω)
        end
    end
    return sum(integrate(ff,f.dΩ))*f.FieldUnit*f.PositionUnit^Dp;
end

############################# getValue for FEM Fields #############################
#From cellfields.jl
function distance(polytope::ExtrusionPolytope, inv_cmap::Gridap.Fields.Field, x::Point)
    extrusion = polytope.extrusion
    isempty(extrusion) && return zero(eltype(x))
    p = inv_cmap(x)
    if all(e == HEX_AXIS for e in extrusion)
        # Boundaries are at `a=0` and `a=1` in each direction
        return maximum(max(0 - a, a - 1) for a in p)
    else all(e == TET_AXIS for e in extrusion)
        # Calculate barycentric coordinates
        λ = Point(p..., 1 - sum(p))
        return maximum(-λ)
    end
end

function get_cell_vector(cache1,x::Point,inv_cmap2,polytope2)
    searchmethod, kdtree, vertex_to_cells, cell_to_ctype, ctype_to_polytope, cell_map, table_cache = cache1
    cell_final=zero(eltype(vertex_to_cells[1]))
    dist_final=Inf
    for (id,dist) in zip(knn(kdtree, SVector(Tuple(x)), searchmethod.num_nearest_vertices, true)...)
        cells = getindex!(table_cache,vertex_to_cells,id)
        for jcell in cells
            jdist=distance(polytope2[jcell], inv_cmap2[jcell], x)
            if (jdist < -1000eps(Float64))
                return jcell
            elseif (jdist ≤ 1000eps(Float64)) && (jdist < dist_final)
                cell_final = jcell
                dist_final = jdist
            end
        end
    end
    return cell_final
end

function get_cell(cache1,x::Point)
    searchmethod, kdtree, vertex_to_cells, cell_to_ctype, ctype_to_polytope, cell_map, table_cache = cache1
    cell_final=zero(eltype(vertex_to_cells[1]))
    dist_final=Inf
    for (id,dist) in zip(knn(kdtree, SVector(Tuple(x)), searchmethod.num_nearest_vertices, true)...)
        cells = getindex!(table_cache,vertex_to_cells,id)
        function cell_distance(cell::Integer)
            ctype = cell_to_ctype[cell]
            polytope = ctype_to_polytope[ctype]
            cmap = cell_map[cell]
            inv_cmap = inverse_map(cmap)
            return distance(polytope, inv_cmap, x)
        end
        for jcell in cells
            jdist = cell_distance(jcell)
            if (jdist < -1000eps(Float64))
                return jcell
            elseif (jdist ≤ 1000eps(Float64)) && (jdist < dist_final)
                cell_final = jcell
                dist_final = jdist
            end
        end
    end
    return cell_final
end

"""
    getValue(f::Gridap.CellField,p::AbstractArray{<:Point})
"""
function getValue(f::Gridap.CellField,p::AbstractArray{<:Point})
    #result=zeros(ComplexF64,size(p))
    Dp=num_dims(get_triangulation(f))
    x=VectorValue(Tuple(0 for j=1:Dp))
    cache1,cache2=Gridap.Fields.return_cache(f,x);
    searchmethod, kdtree, vertex_to_cells, cell_to_ctype, ctype_to_polytope, cell_map, table_cache = cache1
    cell_f_cache, f_cache, cell_f, f₀ = cache2
    #Pour avoir le type de sortie
    cf=getindex!(cell_f_cache, cell_f, 1)
    a=evaluate!(f_cache,cf,p[1])
    result=fill(zero(a),size(p))
    ###
    nt=Threads.nthreads()
    cache1t=Vector{typeof(cache1)}(undef,nt)
    cell_f_cachet=Vector{typeof(cell_f_cache)}(undef,nt)
    f_cachet=Vector{typeof(f_cache)}(undef,nt)
    inv_cmap = inverse_map.(cell_map)
    polytope = ctype_to_polytope[cell_to_ctype]
    for i=1:Threads.nthreads()
        cache1t[i]=(KDTreeSearch(num_nearest_vertices=Dp), kdtree, vertex_to_cells, cell_to_ctype, ctype_to_polytope, cell_map, deepcopy(table_cache))
        cell_f_cachet[i]=deepcopy(cell_f_cache)
        f_cachet[i]=deepcopy(f_cache)
    end
    Threads.@threads for i in eachindex(p)
        cell=get_cell_vector(cache1t[Threads.threadid()],p[i],inv_cmap,polytope)
        if (cell!=0)
            cf = getindex!(cell_f_cachet[Threads.threadid()], cell_f, cell)
            result[i]=evaluate!(f_cachet[Threads.threadid()],cf,p[i])
        end
    end
    return result
end

"""
    getValue(f::Gridap.CellField,p::Point)
"""
function getValue(f::Gridap.CellField,p::Point)
    Dp=num_dims(get_triangulation(f))
    x=VectorValue(Tuple(0.0 for j=1:Dp))
    cache1,cache2=Gridap.Fields.return_cache(f,x);
    cell_f_cache, f_cache, cell_f, f₀ = cache2
    cell=get_cell(cache1,p)
    #Pour avoir le type de sortie
    cf=getindex!(cell_f_cache, cell_f, 1)
    ###
    if (cell!=0)
        cf = getindex!(cell_f_cache, cell_f, cell)
        return evaluate!(f_cache,cf,p)
    else
        a=evaluate!(f_cache,cf,p)
        return zero(a)
    end
end



