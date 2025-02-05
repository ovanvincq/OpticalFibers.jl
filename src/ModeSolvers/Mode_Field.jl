"""
    abstract type Field end
"""
abstract type Field end

Base.:transpose(f::Field) = f
Broadcast.:broadcastable(f::Field)=Ref(f)

"""
    abstract type ScalarField1D <: Field end
"""
abstract type ScalarField1D <: Field end

"""
    abstract type ScalarField2D <: Field end
"""
abstract type ScalarField2D <: Field end

"""
    abstract type VectorField2D <: Field end
"""
abstract type VectorField2D <: Field end

"""
    mutable struct ScalarFieldFunction1D <: ScalarField1D
- nu :: `Int` - Azimuthal number
- E :: `Function`
"""
mutable struct ScalarFieldFunction1D <: ScalarField1D
    nu::Int
    E::Function
    ScalarFieldFunction1D(nu,E)=new(nu,deepcopy(E))
end

Base.:+(f1::ScalarFieldFunction1D, f2::ScalarFieldFunction1D) = ((f1.nu==f2.nu)) ?  ScalarFieldFunction1D(f1.nu,x->f1.E(x)+f2.E(x)) : throw(ArgumentError("Both fields must have the same value of nu"))
Base.:-(f1::ScalarFieldFunction1D, f2::ScalarFieldFunction1D) = ((f1.nu==f2.nu)) ?  ScalarFieldFunction1D(f1.nu,x->f1.E(x)-f2.E(x)) : throw(ArgumentError("Both fields must have the same value of nu"))
Base.:*(k::Number, f::ScalarFieldFunction1D) = ScalarFieldFunction1D(f.nu,x->k*f.E(x));
Base.:*(f::ScalarFieldFunction1D,k::Number) = ScalarFieldFunction1D(f.nu,x->k*f.E(x));
Base.:/(f::ScalarFieldFunction1D,k::Number) = ScalarFieldFunction1D(f.nu,x->f.E(x)/k);
Base.:+(f::ScalarFieldFunction1D) = ScalarFieldFunction1D(f.nu,x->f.E(x));
Base.:-(f::ScalarFieldFunction1D) = ScalarFieldFunction1D(f.nu,x->-f.E(x));
Base.:real(f::ScalarFieldFunction1D) = ScalarFieldFunction1D(f.nu,x->real(f.E(x)));
Base.:imag(f::ScalarFieldFunction1D) = ScalarFieldFunction1D(f.nu,x->imag(f.E(x)));
Base.:conj(f::ScalarFieldFunction1D) = ScalarFieldFunction1D(f.nu,x->conj(f.E(x)));
Base.:abs(f::ScalarFieldFunction1D) = ScalarFieldFunction1D(f.nu,x->abs(f.E(x)));
Base.:abs2(f::ScalarFieldFunction1D) = ScalarFieldFunction1D(f.nu,x->abs2(f.E(x)));

"""
    mutable struct ScalarFieldMatrix1D <: ScalarField1D
- nu :: `Int` - Azymuthal number
- r :: `AbstractVector{<:Real}`
- E :: `AbstractVector{<:Number}`
"""
mutable struct ScalarFieldMatrix1D <: ScalarField1D
    nu::Int
    r::AbstractVector{<:Real}
    E::AbstractVector{<:Number}
    ScalarFieldMatrix1D(nu,r,E)=new(nu,copy(r),copy(E))
end

Base.:+(f1::ScalarFieldMatrix1D, f2::ScalarFieldMatrix1D) = ((f1.nu==f2.nu) && (f1.r==f2.r)) ?  ScalarFieldMatrix1D(f1.nu,f1.r,f1.E+f2.E) : throw(ArgumentError("Both fields must have the same value of nu and the same vector r"))
Base.:-(f1::ScalarFieldMatrix1D, f2::ScalarFieldMatrix1D) = ((f1.nu==f2.nu) && (f1.r==f2.r)) ?  ScalarFieldMatrix1D(f1.nu,f1.r,f1.E-f2.E) : throw(ArgumentError("Both fields must have the same value of nu and the same vector r"))
Base.:*(k::Number, f::ScalarFieldMatrix1D) = ScalarFieldMatrix1D(f.nu,f.r,k*f.E);
Base.:*(f::ScalarFieldMatrix1D,k::Number) = ScalarFieldMatrix1D(f.nu,f.r,k*f.E);
Base.:/(f::ScalarFieldMatrix1D,k::Number) = ScalarFieldMatrix1D(f.nu,f.r,f.E/k);
Base.:+(f::ScalarFieldMatrix1D) = ScalarFieldMatrix1D(f.nu,f.r,f.E);
Base.:-(f::ScalarFieldMatrix1D) = ScalarFieldMatrix1D(f.nu,f.r,-f.E);
Base.:real(f::ScalarFieldMatrix1D) = ScalarFieldMatrix1D(f.nu,f.r,real(f.E));
Base.:imag(f::ScalarFieldMatrix1D) = ScalarFieldMatrix1D(f.nu,f.r,imag(f.E));
Base.:conj(f::ScalarFieldMatrix1D) = ScalarFieldMatrix1D(f.nu,f.r,conj(f.E));
Base.:abs(f::ScalarFieldMatrix1D) = ScalarFieldMatrix1D(f.nu,f.r,abs.(f.E));
Base.:abs2(f::ScalarFieldMatrix1D) = ScalarFieldMatrix1D(f.nu,f.r,abs2.(f.E));

"""
    mutable struct ScalarFieldFEM1D <: ScalarField1D
- nu :: `Int` - Azymuthal number
- dΩ :: `Gridap.CellData.Measure`
- E :: `Gridap.CellField`
"""
mutable struct ScalarFieldFEM1D <: ScalarField1D
    nu::Int
    dΩ::Gridap.CellData.Measure
    E::Gridap.CellField
end

Base.:+(f1::ScalarFieldFEM1D, f2::ScalarFieldFEM1D) = (f1.nu==f2.nu) ?  ScalarFieldFEM1D(f1.nu,f1.dΩ,f1.E+f2.E) : throw(ArgumentError("Both fields must have the same value of nu and the same vector r"))
Base.:-(f1::ScalarFieldFEM1D, f2::ScalarFieldFEM1D) = (f1.nu==f2.nu) ?  ScalarFieldFEM1D(f1.nu,f1.dΩ,f1.E-f2.E) : throw(ArgumentError("Both fields must have the same value of nu and the same vector r"))
Base.:*(k::Number, f::ScalarFieldFEM1D) = ScalarFieldFEM1D(f.nu,f.dΩ,k*f.E);
Base.:*(f::ScalarFieldFEM1D,k::Number) = ScalarFieldFEM1D(f.nu,f.dΩ,k*f.E);
Base.:/(f::ScalarFieldFEM1D,k::Number) = ScalarFieldFEM1D(f.nu,f.dΩ,f.E/k);
Base.:+(f::ScalarFieldFEM1D) = ScalarFieldFEM1D(f.nu,f.dΩ,f.E);
Base.:-(f::ScalarFieldFEM1D) = ScalarFieldFEM1D(f.nu,f.dΩ,-f.E);
Base.:real(f::ScalarFieldFEM1D) = ScalarFieldFEM1D(f.nu,f.dΩ,real(f.E));
Base.:imag(f::ScalarFieldFEM1D) = ScalarFieldFEM1D(f.nu,f.dΩ,imag(f.E));
Base.:conj(f::ScalarFieldFEM1D) = ScalarFieldFEM1D(f.nu,f.dΩ,conj(f.E));
Base.:abs(f::ScalarFieldFEM1D) = ScalarFieldFEM1D(f.nu,f.dΩ,abs(f.E));
Base.:abs2(f::ScalarFieldFEM1D) = ScalarFieldFEM1D(f.nu,f.dΩ,abs2(f.E));

"""
    mutable struct ScalarFieldFunction2D <: ScalarField2D
- E :: `Function`
"""
mutable struct ScalarFieldFunction2D <: ScalarField2D
    E::Function
    ScalarFieldFunction2D(E)=new(deepcopy(E))
end

Base.:+(f1::ScalarFieldFunction2D, f2::ScalarFieldFunction2D) = ScalarFieldFunction2D(x->f1.E(x)+f2.E(x))
Base.:-(f1::ScalarFieldFunction2D, f2::ScalarFieldFunction2D) = ScalarFieldFunction2D(x->f1.E(x)-f2.E(x))
Base.:*(k::Number, f::ScalarFieldFunction2D) = ScalarFieldFunction2D(x->k*f.E(x));
Base.:*(f::ScalarFieldFunction2D,k::Number) = ScalarFieldFunction2D(x->k*f.E(x));
Base.:/(f::ScalarFieldFunction2D,k::Number) = ScalarFieldFunction2D(x->f.E(x)/k);
Base.:+(f::ScalarFieldFunction2D) = ScalarFieldFunction2D(x->f.E(x));
Base.:-(f::ScalarFieldFunction2D) = ScalarFieldFunction2D(x->-f.E(x));
Base.:real(f::ScalarFieldFunction2D) = ScalarFieldFunction2D(x->real(f.E(x)));
Base.:imag(f::ScalarFieldFunction2D) = ScalarFieldFunction2D(x->imag(f.E(x)));
Base.:conj(f::ScalarFieldFunction2D) = ScalarFieldFunction2D(x->conj(f.E(x)));
Base.:abs(f::ScalarFieldFunction2D) = ScalarFieldFunction2D(x->abs(f.E(x)));
Base.:abs2(f::ScalarFieldFunction2D) = ScalarFieldFunction2D(x->abs2(f.E(x)));

"""
    mutable struct ScalarFieldMatrix2D <: ScalarField2D
- x :: `AbstractVector{<:Real}`
- y :: `AbstractVector{<:Real}`
- E :: `AbstractMatrix{<:Number}`
"""
mutable struct ScalarFieldMatrix2D <: ScalarField2D
    x::AbstractVector{<:Real}
    y::AbstractVector{<:Real}
    E::AbstractMatrix{<:Number}
    ScalarFieldMatrix2D(x,y,E)=new(copy(x),copy(y),copy(E))
end

Base.:+(f1::ScalarFieldMatrix2D, f2::ScalarFieldMatrix2D) = ((f1.x==f2.x) && (f1.y==f2.y)) ?  ScalarFieldMatrix2D(f1.x,f1.y,f1.E+f2.E) : throw(ArgumentError("Both fields must have the same vectors x and y"))
Base.:-(f1::ScalarFieldMatrix2D, f2::ScalarFieldMatrix2D) = ((f1.x==f2.x) && (f1.y==f2.y)) ?  ScalarFieldMatrix2D(f1.x,f1.y,f1.E-f2.E) : throw(ArgumentError("Both fields must have the same vectors x and y"))
Base.:*(k::Number, f::ScalarFieldMatrix2D) = ScalarFieldMatrix2D(f.x,f.y,k*f.E);
Base.:*(f::ScalarFieldMatrix2D,k::Number) = ScalarFieldMatrix2D(f.x,f.y,k*f.E);
Base.:/(f::ScalarFieldMatrix2D,k::Number) = ScalarFieldMatrix2D(f.x,f.y,f.E/k);
Base.:+(f::ScalarFieldMatrix2D) = ScalarFieldMatrix2D(f.x,f.y,f.E);
Base.:-(f::ScalarFieldMatrix2D) = ScalarFieldMatrix2D(f.x,f.y,-f.E);
Base.:real(f::ScalarFieldMatrix2D) = ScalarFieldMatrix2D(f.x,f.y,real(f.E));
Base.:imag(f::ScalarFieldMatrix2D) = ScalarFieldMatrix2D(f.x,f.y,imag(f.E));
Base.:conj(f::ScalarFieldMatrix2D) = ScalarFieldMatrix2D(f.x,f.y,conj(f.E));
Base.:abs(f::ScalarFieldMatrix2D) = ScalarFieldMatrix2D(f.x,f.y,abs.(f.E));
Base.:abs2(f::ScalarFieldMatrix2D) = ScalarFieldMatrix2D(f.x,f.y,abs2.(f.E));

"""
    mutable struct ScalarFieldFEM2D <: ScalarField2D
- dΩ :: `Gridap.CellData.Measure`
- E :: `Gridap.CellField`
"""
mutable struct ScalarFieldFEM2D <: ScalarField2D
    dΩ::Gridap.CellData.Measure
    E::Gridap.CellField
end

Base.:+(f1::ScalarFieldFEM2D, f2::ScalarFieldFEM2D) = ScalarFieldFEM2D(f1.dΩ,f1.E+f2.E)
Base.:-(f1::ScalarFieldFEM2D, f2::ScalarFieldFEM2D) = ScalarFieldFEM2D(f1.dΩ,f1.E-f2.E) 
Base.:*(k::Number, f::ScalarFieldFEM2D) = ScalarFieldFEM2D(f.dΩ,k*f.E);
Base.:*(f::ScalarFieldFEM2D,k::Number) = ScalarFieldFEM2D(f.dΩ,k*f.E);
Base.:/(f::ScalarFieldFEM2D,k::Number) = ScalarFieldFEM2D(f.dΩ,f.E/k);
Base.:+(f::ScalarFieldFEM2D) = ScalarFieldFEM2D(f.dΩ,f.E);
Base.:-(f::ScalarFieldFEM2D) = ScalarFieldFEM2D(f.dΩ,-f.E);
Base.:real(f::ScalarFieldFEM2D) = ScalarFieldFEM2D(f.dΩ,real(f.E));
Base.:imag(f::ScalarFieldFEM2D) = ScalarFieldFEM2D(f.dΩ,imag(f.E));
Base.:conj(f::ScalarFieldFEM2D) = ScalarFieldFEM2D(f.dΩ,conj(f.E));
Base.:abs(f::ScalarFieldFEM2D) = ScalarFieldFEM2D(f.dΩ,abs(f.E));
Base.:abs2(f::ScalarFieldFEM2D) = ScalarFieldFEM2D(f.dΩ,abs2(f.E));
#Base.:sum(f::AbstractVector{ScalarFieldFEM2D}) = ScalarFieldFEM2D(f[1].dΩ,sum(getproperty.(f,:E)))

"""
    mutable struct VectorFieldFunction2D <: VectorField2D
- Ex :: `Function`
- Ey :: `Function`
- Ez :: `Function`
- Hx :: `Function`
- Hy :: `Function`
- Hz :: `Function`
"""
mutable struct VectorFieldFunction2D <: VectorField2D
    Ex::Function
    Ey::Function
    Ez::Function
    Hx::Function
    Hy::Function
    Hz::Function
    VectorFieldFunction2D(Ex,Ey,Ez,Hx,Hy,Hz)=new(deepcopy(Ex),deepcopy(Ey),deepcopy(Ez),deepcopy(Hx),deepcopy(Hy),deepcopy(Hz))
end

Base.:+(f1::VectorFieldFunction2D, f2::VectorFieldFunction2D) = VectorFieldFunction2D(x->f1.Ex(x)+f2.Ex(x),x->f1.Ey(x)+f2.Ey(x),x->f1.Ez(x)+f2.Ez(x),x->f1.Hx(x)+f2.Hx(x),x->f1.Hy(x)+f2.Hy(x),x->f1.Hz(x)+f2.Hz(x))
Base.:-(f1::VectorFieldFunction2D, f2::VectorFieldFunction2D) = VectorFieldFunction2D(x->f1.Ex(x)-f2.Ex(x),x->f1.Ey(x)-f2.Ey(x),x->f1.Ez(x)-f2.Ez(x),x->f1.Hx(x)-f2.Hx(x),x->f1.Hy(x)-f2.Hy(x),x->f1.Hz(x)-f2.Hz(x))
Base.:*(k::Number, f::VectorFieldFunction2D) = VectorFieldFunction2D(x->k*f.Ex(x),x->k*f.Ey(x),x->k*f.Ez(x),x->k*f.Hx(x),x->k*f.Hy(x),x->k*f.Hz(x));
Base.:*(f::VectorFieldFunction2D,k::Number) = VectorFieldFunction2D(x->k*f.Ex(x),x->k*f.Ey(x),x->k*f.Ez(x),x->k*f.Hx(x),x->k*f.Hy(x),x->k*f.Hz(x));
Base.:/(f::VectorFieldFunction2D,k::Number) = VectorFieldFunction2D(x->f.Ex(x)/k,x->f.Ey(x)/k,x->f.Ez(x)/k,x->f.Hx(x)/k,x->f.Hy(x)/k,x->f.Hz(x)/k);
Base.:+(f::VectorFieldFunction2D) = VectorFieldFunction2D(x->f.Ex(x),x->f.Ey(x),x->f.Ez(x),x->f.Hx(x),x->f.Hy(x),x->f.Hz(x));
Base.:-(f::VectorFieldFunction2D) = VectorFieldFunction2D(x->-f.Ex(x),x->-f.Ey(x),x->-f.Ez(x),x->-f.Hx(x),x->-f.Hy(x),x->-f.Hz(x));
Base.:real(f::VectorFieldFunction2D) = VectorFieldFunction2D(x->real(f.Ex(x)),x->real(f.Ey(x)),x->real(f.Ez(x)),x->real(f.Hx(x)),x->real(f.Hy(x)),x->real(f.Hz(x)));
Base.:imag(f::VectorFieldFunction2D) = VectorFieldFunction2D(x->imag(f.Ex(x)),x->imag(f.Ey(x)),x->imag(f.Ez(x)),x->imag(f.Hx(x)),x->imag(f.Hy(x)),x->imag(f.Hz(x)));
Base.:conj(f::VectorFieldFunction2D) = VectorFieldFunction2D(x->conj(f.Ex(x)),x->conj(f.Ey(x)),x->conj(f.E(x)),x->conj(f.Hx(x)),x->conj(f.Hy(x)),x->conj(f.Hz(x)));

"""
    mutable struct VectorFieldMatrix2D <: VectorField2D
- x :: `AbstractVector{<:Real}`
- y :: `AbstractVector{<:Real}`
- Ex :: `AbstractMatrix{<:Number}`
- Ey :: `AbstractMatrix{<:Number}`
- Ez :: `AbstractMatrix{<:Number}`
- Hx :: `AbstractMatrix{<:Number}`
- Hy :: `AbstractMatrix{<:Number}`
- Hz :: `AbstractMatrix{<:Number}`
"""
mutable struct VectorFieldMatrix2D <: VectorField2D
    x::AbstractVector{<:Real}
    y::AbstractVector{<:Real}
    Ex::AbstractMatrix{<:Number}
    Ey::AbstractMatrix{<:Number}
    Ez::AbstractMatrix{<:Number}
    Hx::AbstractMatrix{<:Number}
    Hy::AbstractMatrix{<:Number}
    Hz::AbstractMatrix{<:Number}
    VectorFieldMatrix2D(x,y,Ex,Ey,Ez,Hx,Hy,Hz)=new(copy(x),copy(y),copy(Ex),copy(Ey),copy(Ez),copy(Hx),copy(Hy),copy(Hz))
end

Base.:+(f1::VectorFieldMatrix2D, f2::VectorFieldMatrix2D) = ((f1.x==f2.x) && (f1.y==f2.y)) ?  VectorFieldMatrix2D(f1.x,f1.y,f1.Ex+f2.Ex,f1.Ey+f2.Ey,f1.Ez+f2.Ez,f1.Hx+f2.Hx,f1.Hy+f2.Hy,f1.Hz+f2.Hz) : throw(ArgumentError("Both fields must have the same vectors x and y"))
Base.:-(f1::VectorFieldMatrix2D, f2::VectorFieldMatrix2D) = ((f1.x==f2.x) && (f1.y==f2.y)) ?  VectorFieldMatrix2D(f1.x,f1.y,f1.Ex-f2.Ex,f1.Ey-f2.Ey,f1.Ez-f2.Ez,f1.Hx-f2.Hx,f1.Hy-f2.Hy,f1.Hz-f2.Hz) : throw(ArgumentError("Both fields must have the same vectors x and y"))
Base.:*(k::Number, f::VectorFieldMatrix2D) = VectorFieldMatrix2D(f.x,f.y,k*f.Ex,k*f.Ey,k*f.Ez,k*f.Hx,k*f.Hy,k*f.Hz);
Base.:*(f::VectorFieldMatrix2D,k::Number) = VectorFieldMatrix2D(f.x,f.y,k*f.Ex,k*f.Ey,k*f.Ez,k*f.Hx,k*f.Hy,k*f.Hz);
Base.:/(f::VectorFieldMatrix2D,k::Number) = VectorFieldMatrix2D(f.x,f.y,f.Ex/k,f.Ey/k,f.Ez/k,f.Hx/k,f.Hy/k,f.Hz/k);
Base.:+(f::VectorFieldMatrix2D) = VectorFieldMatrix2D(f.x,f.y,f.Ex,f.Ey,f.Ez,f.Hx,f.Hy,f.Hz);
Base.:-(f::VectorFieldMatrix2D) = VectorFieldMatrix2D(f.x,f.y,-f.Ex,-f.Ey,-f.Ez,-f.Hx,-f.Hy,-f.Hz);
Base.:real(f::VectorFieldMatrix2D) = VectorFieldMatrix2D(f.x,f.y,real(f.Ex),real(f.Ey),real(f.Ez),real(f.Hx),real(f.Hy),real(f.Hz));
Base.:imag(f::VectorFieldMatrix2D) = VectorFieldMatrix2D(f.x,f.y,imag(f.Ex),imag(f.Ey),imag(f.Ez),imag(f.Hx),imag(f.Hy),imag(f.Hz));
Base.:conj(f::VectorFieldMatrix2D) = VectorFieldMatrix2D(f.x,f.y,conj(f.Ex),conj(f.Ey),conj(f.Ez),conj(f.Hx),conj(f.Hy),conj(f.Hz));

"""
    mutable struct VectorFieldMatrix2D <: VectorField2D
- dΩ :: `Gridap.CellData.Measure`
- Ex :: `Gridap.CellField`
- Ey :: `Gridap.CellField`
- Ez :: `Gridap.CellField`
- Hx :: `Gridap.CellField`
- Hy :: `Gridap.CellField`
- Hz :: `Gridap.CellField`
"""
mutable struct VectorFieldFEM2D <: VectorField2D
    dΩ::Gridap.CellData.Measure
    Ex::Gridap.CellField
    Ey::Gridap.CellField
    Ez::Gridap.CellField
    Hx::Gridap.CellField
    Hy::Gridap.CellField
    Hz::Gridap.CellField
end

Base.:+(f1::VectorFieldFEM2D, f2::VectorFieldFEM2D) = VectorFieldFEM2D(f1.dΩ,f1.Ex+f2.Ex,f1.Ey+f2.Ey,f1.Ez+f2.Ez,f1.Hx+f2.Hx,f1.Hy+f2.Hy,f1.Hz+f2.Hz)
Base.:-(f1::VectorFieldFEM2D, f2::VectorFieldFEM2D) = VectorFieldFEM2D(f1.dΩ,f1.Ex-f2.Ex,f1.Ey-f2.Ey,f1.Ez-f2.Ez,f1.Hx-f2.Hx,f1.Hy-f2.Hy,f1.Hz-f2.Hz) 
Base.:*(k::Number, f::VectorFieldFEM2D) = VectorFieldFEM2D(f.dΩ,k*f.Ex,k*f.Ey,k*f.Ez,k*f.Hx,k*f.Hy,k*f.Hz);
Base.:*(f::VectorFieldFEM2D,k::Number) = VectorFieldFEM2D(f.dΩ,k*f.Ex,k*f.Ey,k*f.Ez,k*f.Hx,k*f.Hy,k*f.Hz);
Base.:/(f::VectorFieldFEM2D,k::Number) = VectorFieldFEM2D(f.dΩ,f.Ex/k,f.Ey/k,f.Ez/k,f.Hx/k,f.Hy/k,f.Hz/k);
Base.:+(f::VectorFieldFEM2D) = VectorFieldFEM2D(f.dΩ,f.Ex,f.Ey,f.Ez,f.Hx,f.Hy,f.Hz);
Base.:-(f::VectorFieldFEM2D) = VectorFieldFEM2D(f.dΩ,-f.Ex,-f.Ey,-f.Ez,-f.Hx,-f.Hy,-f.Hz);
Base.:real(f::VectorFieldFEM2D) = VectorFieldFEM2D(f.dΩ,real(f.Ex),real(f.Ey),real(f.Ez),real(f.Hx),real(f.Hy),real(f.Hz));
Base.:imag(f::VectorFieldFEM2D) = VectorFieldFEM2D(f.dΩ,imag(f.Ex),imag(f.Ey),imag(f.Ez),imag(f.Hx),imag(f.Hy),imag(f.Hz));
Base.:conj(f::VectorFieldFEM2D) = VectorFieldFEM2D(f.dΩ,conj(f.Ex),conj(f.Ey),conj(f.Ez),conj(f.Hx),conj(f.Hy),conj(f.Hz));

"""
    isvalidField(f::Field)
"""
function isvalidField(f::Field)
    if isa(f,ScalarFieldFunction1D)
        return (1 in nb_args(f.E)) ? true : false
    elseif isa(f,ScalarFieldFunction2D)
        return (1 in nb_args(f.E)) ? true : false
    elseif isa(f,VectorFieldFunction2D)
        return ((1 in nb_args(f.Ex)) && (1 in nb_args(f.Ey)) && (1 in nb_args(f.Ez)) && (1 in nb_args(f.Hx)) && (1 in nb_args(f.Hy)) && (1 in nb_args(f.Hz))) ? true : false
     elseif isa(f,ScalarFieldMatrix1D)
        return (length(f.r)==length(f.E)) ? true : false
    elseif isa(f,ScalarFieldMatrix2D)
        return ((length(f.x)==size(f.E,1)) && (length(f.y)==size(f.E,2))) ? true : false
    elseif isa(f,VectorFieldMatrix2D)
        return ((size(f.Ex)==size(f.Ey)==size(f.Ez)==size(f.Hx)==size(f.Hy)==size(f.Hz)) && (length(f.x)==size(f.Ex,1)) && (length(f.y)==size(f.Ex,2))) ? true : false
    elseif isa(f,ScalarFieldFEM1D)
        return (get_triangulation(f.dΩ.quad)==get_triangulation(f.E)) ? true : false
    elseif isa(f,ScalarFieldFEM2D)
        return (get_triangulation(f.dΩ.quad)==get_triangulation(f.E)) ? true : false
    elseif isa(f,VectorFieldFEM2D)
        return (get_triangulation(f.dΩ.quad)==get_triangulation(f.Ex)==get_triangulation(f.Ey)==get_triangulation(f.Ez)==get_triangulation(f.Hx)==get_triangulation(f.Hy)==get_triangulation(f.Hz)) ? true : false
    end
end

"""
    struct Mode{T<:Union{Field,Nothing}}

- Name :: `String`
- neff :: `Number`
- lambda :: `Real`
- field :: `Field` or `Nothing`
"""
struct Mode{T<:Union{Field,Nothing}}
    Name::String
    neff::Number
    lambda::Real
    field::T
    Mode(Name,neff,lambda,field)=(isnothing(field) || isvalidField(field)) ? new{typeof(field)}(Name,neff,lambda,field) : new{Nothing}(Name,neff,lambda,nothing)
    Mode(Name,neff,lambda)=new{Nothing}(Name,neff,lambda,nothing)
end

Base.:transpose(m::Mode) = m
Broadcast.:broadcastable(m::Mode)=Ref(m)

function Base.show(io::IO, ::MIME"text/plain",m::Mode) 
    print(io,"Name = ",m.Name,"\nneff = ",m.neff,"\nlambda = ",m.lambda,"\nfield = ",typeof(m.field));
end
Base.show(io::IO, m::Mode) = print(io,"[",m.Name,",",m.neff,",",m.lambda,",",typeof(m.field),"]") 

"""
    isvalidField(m::Mode)
"""
function isvalidMode(m::Mode)
    return isnothing(m.field) ? true : isvalidField(m.field)
end

############################# getValue for FEM Fields #############################
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
    result=zeros(ComplexF64,size(p))
    dim=num_dims(get_triangulation(f));
    if (dim==1)
        cache1,cache2=Gridap.Fields.return_cache(f,Point(0.0,));
    else
        cache1,cache2=Gridap.Fields.return_cache(f,Point(0.0,0.0));
    end
    searchmethod, kdtree, vertex_to_cells, cell_to_ctype, ctype_to_polytope, cell_map, table_cache = cache1
    cell_f_cache, f_cache, cell_f, f₀ = cache2
    nt=Threads.nthreads()
    cache1t=Vector{typeof(cache1)}(undef,nt)
    cell_f_cachet=Vector{typeof(cell_f_cache)}(undef,nt)
    f_cachet=Vector{typeof(f_cache)}(undef,nt)
    inv_cmap = inverse_map.(cell_map)
    polytope = ctype_to_polytope[cell_to_ctype]
    for i=1:Threads.nthreads()
        cache1t[i]=(KDTreeSearch(num_nearest_vertices=dim), kdtree, vertex_to_cells, cell_to_ctype, ctype_to_polytope, cell_map, deepcopy(table_cache))
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
    if (num_dims(get_triangulation(f))==2)
        cache1,cache2=Gridap.Fields.return_cache(f,Point(0.0,0.0));
    else
        cache1,cache2=Gridap.Fields.return_cache(f,Point(0.0,));
    end
    cell_f_cache, f_cache, cell_f, f₀ = cache2
    cell=get_cell(cache1,p)
    if (cell!=0)
        cf = getindex!(cell_f_cache, cell_f, cell)
        return ComplexF64(evaluate!(f_cache,cf,p))
    else
        return zero(ComplexF64)
    end
end

############################# computeField ######################
"""
    computeField(f::ScalarField1D,r::Union{Real,AbstractArray{<:Real}})
"""
function computeField(f::ScalarField1D,r::Union{Real,AbstractArray{<:Real}})
    Gridap.Helpers.@notimplemented
end

"""
    computeField(f::ScalarField2D,x::Union{Real,AbstractArray{<:Real}},y::Union{Real,AbstractArray{<:Real}})
"""
function computeField(f::ScalarField2D,x::Union{Real,AbstractArray{<:Real}},y::Union{Real,AbstractArray{<:Real}})
    Gridap.Helpers.@notimplemented
end

"""
    computeField(f::VectorField2D,x::Union{Real,AbstractArray{<:Real}},y::Union{Real,AbstractArray{<:Real}},component::Symbol=:Ex)
"""
function computeField(f::VectorField2D,x::Union{Real,AbstractArray{<:Real}},y::Union{Real,AbstractArray{<:Real}},component::Symbol=:Ex)
    Gridap.Helpers.@notimplemented
end  


function computeField(f::ScalarFieldFEM1D,r::Union{Real,AbstractArray{<:Real}})
    return getValue(f.E,Point.(r,))
end

function computeField(f::ScalarFieldFEM2D,x::Union{Real,AbstractArray{<:Real}},y::Union{Real,AbstractArray{<:Real}})
    return getValue(f.E,Point.(x,y))
end

function computeField(f::VectorFieldFEM2D,x::Union{Real,AbstractArray{<:Real}},y::Union{Real,AbstractArray{<:Real}},component::Symbol=:Ex)
    if component==:Ex
        return getValue(f.Ex,Point.(x,y))
    elseif component==:Ey
        return getValue(f.Ey,Point.(x,y))
    elseif component==:Ez
        return getValue(f.Ez,Point.(x,y))
    elseif component==:Hx
        return getValue(f.Hx,Point.(x,y))
    elseif component==:Hy
        return getValue(f.Hy,Point.(x,y))
    elseif component==:Hz
        return getValue(f.Hz,Point.(x,y))
    else
        throw(ArgumentError("component must be :Ex, :Ey, :Ez, :Hx, :Hy or :Hz"))
    end
end

function computeField(f::ScalarFieldFunction1D,r::Union{Real,AbstractArray{<:Real}})
    return f.E.(Point.(r,))
end


function computeField(f::ScalarFieldFunction2D,x::Union{Real,AbstractArray{<:Real}},y::Union{Real,AbstractArray{<:Real}})
    return f.E.(Point.(x,y))
end

function computeField(f::VectorFieldFunction2D,x::Union{Real,AbstractArray{<:Real}},y::Union{Real,AbstractArray{<:Real}},component::Symbol=:Ex)
    if component==:Ex
        return f.Ex.(Point.(x,y))
    elseif component==:Ey
        return f.Ey.(Point.(x,y))
    elseif component==:Ez
        return f.Ez.(Point.(x,y))
    elseif component==:Hx
        return f.Hx.(Point.(x,y))
    elseif component==:Hy
        return f.Hy.(Point.(x,y))
    elseif component==:Hz
        return f.Hz.(Point.(x,y))
    else
        throw(ArgumentError("component must be :Ex, :Ey, :Ez, :Hx, :Hy or :Hz"))
    end
end


function computeField(f::ScalarFieldMatrix1D,r::Union{Real,AbstractArray{<:Real}})
    interp=LinearInterpolation(f.r,f.E,extrapolation_bc=0);
    return interp.(r);
end

function computeField(f::ScalarFieldMatrix2D,x::Union{Real,AbstractArray{<:Real}},y::Union{Real,AbstractArray{<:Real}})
    interp=LinearInterpolation((f.y,f.x),f.E,extrapolation_bc=0);
    return interp.(x,y);
end

function computeField(f::VectorFieldMatrix2D,x::Union{Real,AbstractArray{<:Real}},y::Union{Real,AbstractArray{<:Real}},component::Symbol=:Ex)
    if component==:Ex
        interp=LinearInterpolation((f.y,f.x),f.Ex,extrapolation_bc=0);
        return interp.(x,y);
    elseif component==:Ey
        interp=LinearInterpolation((f.y,f.x),f.Ey,extrapolation_bc=0);
        return interp.(x,y);
    elseif component==:Ez
        interp=LinearInterpolation((f.y,f.x),f.Ez,extrapolation_bc=0);
        return interp.(x,y);
    elseif component==:Hx
        interp=LinearInterpolation((f.y,f.x),f.Hx,extrapolation_bc=0);
        return interp.(x,y);
    elseif component==:Hy
        interp=LinearInterpolation((f.y,f.x),f.Hy,extrapolation_bc=0);
        return interp.(x,y);
    elseif component==:Hz
        interp=LinearInterpolation((f.y,f.x),f.Hz,extrapolation_bc=0);
        return interp.(x,y);
    else
        throw(ArgumentError("component must be :Ex, :Ey, :Ez, :Hx, :Hy or :Hz"))
    end
end

"""
    computeField(m::Mode{<:ScalarField1D},r::Union{Real,AbstractArray{<:Real}},z::Real=0)
"""
function computeField(m::Mode{<:ScalarField1D},r::Union{Real,AbstractArray{<:Real}},z::Real=0)
    if (z==0)
        return computeField(m.field,r)
    else
        computeField(m.field,r).*exp(im*z*2*pi/m.lambda*m.neff)
    end
end

"""
    computeField(m::Mode{<:ScalarField2D},x::Union{Real,AbstractArray{<:Real}},y::Union{Real,AbstractArray{<:Real}},z::Real=0)
"""
function computeField(m::Mode{<:ScalarField2D},x::Union{Real,AbstractArray{<:Real}},y::Union{Real,AbstractArray{<:Real}},z::Real=0)
    if (z==0)
        return computeField(m.field,x,y)
    else
        return computeField(m.field,x,y).*exp(im*z*2*pi/m.lambda*m.neff)
    end
end

"""
    computeField(m::Mode{<:VectorField2D},x::Union{Real,AbstractArray{<:Real}},y::Union{Real,AbstractArray{<:Real}},component::Symbol=:Ex,z::Real=0)
"""
function computeField(m::Mode{<:VectorField2D},x::Union{Real,AbstractArray{<:Real}},y::Union{Real,AbstractArray{<:Real}},component::Symbol=:Ex,z::Real=0)
    if (z==0)
        return computeField(m.field,x,y,component)
    else
        return computeField(m.field,x,y,component).*exp(im*z*2*pi/m.lambda*m.neff)
    end
end

############################# Field conversion #############################
ScalarFieldMatrix1D(f::ScalarField1D,r::AbstractVector{<:Real}) = ScalarFieldMatrix1D(f.nu,r,computeField(f,r))
Mode{ScalarFieldMatrix1D}(m::Mode{<:ScalarField1D},r::AbstractVector{<:Real}) = Mode(m.Name,m.neff,m.lambda,ScalarFieldMatrix1D(m.field,r))

ScalarFieldMatrix2D(f::ScalarField2D,x::AbstractVector{<:Real},y::AbstractVector{<:Real}) = ScalarFieldMatrix2D(x,y,computeField(f,x,y'))
Mode{ScalarFieldMatrix2D}(m::Mode{<:ScalarField2D},x::AbstractVector{<:Real},y::AbstractVector{<:Real}) = Mode(m.Name,m.neff,m.lambda,ScalarFieldMatrix2D(m.field,x,y))

VectorFieldMatrix2D(f::VectorField2D,x::AbstractVector{<:Real},y::AbstractVector{<:Real}) = VectorFieldMatrix2D(x,y,computeField(f,x,y',:Ex),computeField(f,x,y',:Ey),computeField(f,x,y',:Ez),computeField(f,x,y',:Hx),computeField(f,x,y',:Hy),computeField(f,x,y',:Hz))
Mode{VectorFieldMatrix2D}(m::Mode{<:VectorField2D},x::AbstractVector{<:Real},y::AbstractVector{<:Real}) = Mode(m.Name,m.neff,m.lambda,VectorFieldMatrix2D(m.field,x,y))

"""
    convertToMatrix(f::ScalarField1D,r::AbstractVector{<:Real})
"""
convertToMatrix(f::ScalarField1D,r::AbstractVector{<:Real}) = return ScalarFieldMatrix1D(f,r)
"""
    convertToMatrix(m::Mode{<:ScalarField1D},r::AbstractVector{<:Real})
"""
convertToMatrix(m::Mode{<:ScalarField1D},r::AbstractVector{<:Real}) = return Mode{ScalarFieldMatrix1D}(m,r)
"""
    convertToMatrix(f::ScalarField2D,x::AbstractVector{<:Real},y::AbstractVector{<:Real})
"""
convertToMatrix(f::ScalarField2D,x::AbstractVector{<:Real},y::AbstractVector{<:Real}) = return ScalarFieldMatrix2D(f,x,y)
"""
    convertToMatrix(m::Mode{<:ScalarField2D},x::AbstractVector{<:Real},y::AbstractVector{<:Real})
"""
convertToMatrix(m::Mode{<:ScalarField2D},x::AbstractVector{<:Real},y::AbstractVector{<:Real}) = return Mode{ScalarFieldMatrix2D}(m,x,y)
"""
    convertToMatrix(f::VectorField2D,x::AbstractVector{<:Real},y::AbstractVector{<:Real})
"""
convertToMatrix(f::VectorField2D,x::AbstractVector{<:Real},y::AbstractVector{<:Real}) = return VectorFieldMatrix2D(f,x,y)
"""
    convertToMatrix(m::Mode{<:VectorField2D},x::AbstractVector{<:Real},y::AbstractVector{<:Real})
"""
convertToMatrix(m::Mode{<:VectorField2D},x::AbstractVector{<:Real},y::AbstractVector{<:Real}) = return Mode{VectorFieldMatrix2D}(m,x,y)

ScalarFieldFunction2D(f::ScalarFieldFunction1D,angle::Real=0) = ScalarFieldFunction2D(x->f.E(hypot(x[1],x[2]))*cos(f.nu*atan(x[2],x[1])+angle/180*pi))
Mode{ScalarFieldFunction2D}(m::Mode{ScalarFieldFunction1D},angle::Real=0) = Mode(m.Name,m.neff,m.lambda,ScalarFieldFunction2D(m.field,angle))

ScalarFieldMatrix2D(f::ScalarFieldMatrix1D,angle::Real=0) = ScalarFieldMatrix2D(vcat(-reverse(f.r[2:end]),f.r),vcat(-reverse(f.r[2:end]),f.r),computeField(f,hypot.(vcat(-reverse(f.r[2:end]),f.r),vcat(-reverse(f.r[2:end]),f.r)')).*(cos.(angle/180*pi+f.nu*atan.(vcat(-reverse(f.r[2:end]),f.r)',vcat(-reverse(f.r[2:end]),f.r)))))
Mode{ScalarFieldMatrix2D}(m::Mode{ScalarFieldMatrix1D},angle::Real=0) = Mode(m.Name,m.neff,m.lambda,ScalarFieldMatrix2D(m.field,angle))

"""
    convertTo2D(f::ScalarFieldFunction1D,angle::Real=0)
"""
convertTo2D(f::ScalarFieldFunction1D,angle::Real=0) = return ScalarFieldFunction2D(f,angle)
"""
    convertTo2D(m::Mode{ScalarFieldFunction1D},angle::Real=0)
"""
convertTo2D(m::Mode{ScalarFieldFunction1D},angle::Real=0) = return Mode{ScalarFieldFunction2D}(m,angle)
"""
    convertTo2D(f::ScalarFieldMatrix1D,angle::Real=0)
"""
convertTo2D(f::ScalarFieldMatrix1D,angle::Real=0) = return ScalarFieldMatrix2D(f,angle)
"""
    convertTo2D(m::Mode{ScalarFieldMatrix1D},angle::Real=0)
"""
convertTo2D(m::Mode{ScalarFieldMatrix1D},angle::Real=0) = return Mode{ScalarFieldMatrix2D}(m,angle)


Mode{VectorFieldFunction2D}(m::Mode{ScalarFieldFunction2D},angle::Real=0) = return Mode(m.Name,m.neff,m.lambda,VectorFieldFunction2D(x->m.field.E(x)*cosd(angle),x->m.field.E(x)*sind(angle),x->0,x->-m.field.E(x)*sind(angle)*m.neff/c/mu0,x->m.field.E(x)*cosd(angle)*m.neff/c/mu0,x->0))

Mode{VectorFieldMatrix2D}(m::Mode{ScalarFieldMatrix2D},angle::Real=0) = return Mode(m.Name,m.neff,m.lambda,VectorFieldMatrix2D(m.field.x,m.field.y,m.field.E*cosd(angle),m.field.E*sind(angle),zero(m.field.E),-m.field.E*sind(angle)*m.neff/c/mu0,m.field.E*cosd(angle)*m.neff/c/mu0,zero(m.field.E)))

Mode{VectorFieldFEM2D}(m::Mode{ScalarFieldFEM2D},angle::Real=0) = return Mode(m.Name,m.neff,m.lambda,VectorFieldFEM2D(m.field.dΩ,m.field.E*cosd(angle),m.field.E*sind(angle),CellField(0,get_triangulation(m.field.E)),-m.field.E*sind(angle)*m.neff/c/mu0,m.field.E*cosd(angle)*m.neff/c/mu0,CellField(0,get_triangulation(m.field.E))))

"""
    convertToVector(m::Mode{<:ScalarField2D},angle::Real=0)
"""
function convertToVector(m::Mode{<:ScalarField2D},angle::Real=0)
    Gridap.Helpers.@notimplemented
end

convertToVector(m::Mode{ScalarFieldFunction2D},angle::Real=0) = return Mode{VectorFieldFunction2D}(m,angle)
convertToVector(m::Mode{ScalarFieldMatrix2D},angle::Real=0) = return Mode{VectorFieldMatrix2D}(m,angle)
convertToVector(m::Mode{ScalarFieldFEM2D},angle::Real=0) = return Mode{VectorFieldFEM2D}(m,angle)


############################# GetField->Propagation ######################
"""
    getField(m::Mode,z::Real=0)
"""
function getField(m::Mode,z::Real=0)
    if (z==0)
        return m.field
    else
        return m.field*exp(im*z*2*pi/m.lambda*m.neff)
    end
end

############################# Poynting Vector #############################
"""
    PoyntingVector(f::VectorField2D)
"""
function PoyntingVector(f::VectorField2D)
    Gridap.Helpers.@notimplemented
end

"""
    PoyntingVector(m::Mode{<:ScalarField2D})
"""
function PoyntingVector(m::Mode{<:ScalarField2D})
    Gridap.Helpers.@notimplemented
end


function PoyntingVector(f::VectorFieldFEM2D)
    return 0.5*real(f.Ey*conj(f.Hz)-f.Ez*conj(f.Hy)),0.5*real(f.Ez*conj(f.Hx)-f.Ex*conj(f.Hz)),0.5*real(f.Ex*conj(f.Hy)-f.Ey*conj(f.Hx))
end

function PoyntingVector(f::VectorFieldMatrix2D)
    return 0.5*real(f.Ey.*conj(f.Hz)-f.Ez.*conj(f.Hy)),0.5*real(f.Ez.*conj(f.Hx)-f.Ex.*conj(f.Hz)),0.5*real(f.Ex.*conj(f.Hy)-f.Ey.*conj(f.Hx))
end

function PoyntingVector(f::VectorFieldFunction2D)
    return x->0.5*real(f.Ey(x)*conj(f.Hz(x))-f.Ez(x)*conj(f.Hy(x))),x->0.5*real(f.Ez(x)*conj(f.Hx(x))-f.Ex(x)*conj(f.Hz(x))),x->0.5*real(f.Ex(x)*conj(f.Hy(x))-f.Ey(x)*conj(f.Hx(x)))
end

function PoyntingVector(m::Mode{ScalarFieldFEM2D})
    Ω=get_triangulation(m[1].field.E);
    return CellField(0,Ω),CellField(0,Ω),0.5*real(m.neff/c/mu0*abs2(m.field.E))
end

function PoyntingVector(m::Mode{ScalarFieldMatrix2D})
    return zeros(size(m.field.E)),zeros(size(m.field.E)),0.5*real(m.neff/c/mu0*abs2.(m.field.E))
end

function PoyntingVector(m::Mode{ScalarFieldFunction2D})
    return x->0.0,x->0.0,x->0.5*real(m.neff/c/mu0*abs2(m.field.E(x)))
end

"""
    PoyntingVector(m::Mode{<:VectorField2D})
"""
PoyntingVector(m::Mode{<:VectorField2D}) = PoyntingVector(m.field)

############################# Effective area #############################
"""
    Aeff(m::Mode{ScalarFieldFunction1D};kwargs...)
"""
function Aeff(m::Mode{ScalarFieldFunction1D};kwargs...)
    f2=x->2*pi*x[1]*abs2(m.field.E(x))
    f4=x->2*pi*x[1]*abs2(m.field.E(x))^2
    E2=integrate1D(f2;kwargs...)[1]
    E4=integrate1D(f4;kwargs...)[1]
    if (m.field.nu!=0)
        E2=E2*0.5;
        E4=E4*3.0/8.0;
    end
    return E2^2/E4;
end

"""
    Aeff(m::Mode{ScalarFieldFunction2D};kwargs...)
"""
function Aeff(m::Mode{ScalarFieldFunction2D};kwargs...)
    f2=x->abs2(m.field.E(x))
    f4=x->(abs2(m.field.E(x)))^2
    E2=integrate2D(f2;kwargs...)[1]
    E4=integrate2D(f4;kwargs...)[1]
    return E2^2/E4;
end

"""
    Aeff(m::Mode{VectorFieldFunction2D},n0::Union{Real,Function}=0;kwargs...)
"""
function Aeff(m::Mode{VectorFieldFunction2D},n0::Union{Real,Function}=0;kwargs...)
    if isa(n0,Function)
        if (first(methods(n0)).nargs!=2)
            throw(DomainError(n0, "The refractive index function n0 must have 1 argument (x,y)"));
        end
        f2=x->m.field.Ex(x)*conj(m.field.Hy(x))-m.field.Ey(x)*conj(m.field.Hx(x));
        f4_1=x->(abs2(abs2(m.field.Ex(x))+abs2.(m.field.Ey(x))+abs2.(m.field.Ez(x))))*(real(n0(x))^2);
        f4_2=x->(abs2(m.field.Ex(x)^2+m.field.Ey(x)^2+m.field.Ez(x)^2))*(real(n0(x))^2);
        E2=integrate2D(f2;kwargs...)[1]
        E4_1=integrate2D(f4_1;kwargs...)[1]
        E4_2=integrate2D(f4_2;kwargs...)[1]
        return real(E2)^2/E4_1*mu0/eps0,real(E2)^2/E4_2*mu0/eps0;
    else
        f2=x->m.field.Ex(x)*conj(m.field.Hy(x))-m.field.Ey(x)*conj(m.field.Hx(x));
        f4_1=x->(abs2(abs2(m.field.Ex(x))+abs2.(m.field.Ey(x))+abs2.(m.field.Ez(x))));
        f4_2=x->(abs2(m.field.Ex(x)^2+m.field.Ey(x)^2+m.field.Ez(x)^2));
        E2=integrate2D(f2;kwargs...)[1]
        E4_1=integrate2D(f4_1;kwargs...)[1]
        E4_2=integrate2D(f4_2;kwargs...)[1]
        if (n0==0)
            return real(E2)^2/E4_1*mu0/eps0/(real(m.neff))^2,real(E2)^2/E4_2*mu0/eps0/(real(m.neff))^2;
        else
            return real(E2)^2/E4_1*mu0/eps0/real(n0)^2,real(E2)^2/E4_2*mu0/eps0/real(n0)^2;
        end
    end
end

"""
    Aeff(m::Mode{ScalarFieldMatrix1D})
"""
function Aeff(m::Mode{ScalarFieldMatrix1D})
    E2=trapz(m.field.r,2*pi*m.field.r.*abs2.(m.field.E));
    E4=trapz(m.field.r,2*pi*m.field.r.*(abs2.(m.field.E)).^2);
    if (m.field.nu!=0)
        E2=E2*0.5;
        E4=E4*3.0/8.0;
    end
    return E2^2/E4;
end

"""
    Aeff(m::Mode{ScalarFieldMatrix2D})
"""
function Aeff(m::Mode{ScalarFieldMatrix2D})
    E2=trapz((m.field.x,m.field.y),abs2.(m.field.E));
    E4=trapz((m.field.x,m.field.y),(abs2.(m.field.E)).^2);
    return E2^2/E4;
end

"""
    Aeff(m::Mode{VectorFieldMatrix2D},n0::Union{Real,Function}=0)
"""
function Aeff(m::Mode{VectorFieldMatrix2D},n0::Union{Real,Function}=0)
    if isa(n0,Function)
        if (first(methods(n0)).nargs!=2)
            throw(DomainError(n0, "The refractive index function n0 must have 1 argument (x,y)"));
        end
        E2=trapz((m.field.x,m.field.y),m.field.Ex.*conj(m.field.Hy)-m.field.Ey.*conj(m.field.Hx));
        E4_1=trapz((m.field.x,m.field.y),(abs2.(abs2.(m.field.Ex)+abs2.(m.field.Ey)+abs2.(m.field.Ez))).*(real(n0.(Point.(m.field.x,m.field.y')).^2)));
        E4_2=trapz((m.field.x,m.field.y),(abs2.(m.field.Ex.^2+m.field.Ey.^2+m.field.Ez.^2)).*(real(n0.(Point.(m.field.x,m.field.y')).^2)));
        return real(E2)^2/E4_1*mu0/eps0,real(E2)^2/E4_2*mu0/eps0;
    else
        E2=trapz((m.field.x,m.field.y),m.field.Ex.*conj(m.field.Hy)-m.field.Ey.*conj(m.field.Hx));
        E4_1=trapz((m.field.x,m.field.y),(abs2.(abs2.(m.field.Ex)+abs2.(m.field.Ey)+abs2.(m.field.Ez))));
        E4_2=trapz((m.field.x,m.field.y),(abs2.(m.field.Ex.^2+m.field.Ey.^2+m.field.Ez.^2)));
        if (n0==0)
            return real(E2)^2/E4_1*mu0/eps0/(real(m.neff))^2,real(E2)^2/E4_2*mu0/eps0/(real(m.neff))^2;
        else
            return real(E2)^2/E4_1*mu0/eps0/real(n0)^2,real(E2)^2/E4_2*mu0/eps0/real(n0)^2;
        end
    end
end

"""
    Aeff(m::Mode{ScalarFieldFEM1D})
"""
function Aeff(m::Mode{ScalarFieldFEM1D})
    r=x->x[1];
    E2_task=Threads.@spawn sum(integrate(abs2(m.field.E)*2*pi*r,m.field.dΩ));
    E4_task=Threads.@spawn sum(integrate(abs2(abs2(m.field.E))*2*pi*r,m.field.dΩ));
    E2=fetch(E2_task);
    E4=fetch(E4_task);
    if (m.field.nu!=0)
        E2=E2*0.5;
        E4=E4*3.0/8.0;
    end
    return E2^2/E4;
end

"""
    Aeff(m::Mode{ScalarFieldFEM2D})
"""
function Aeff(m::Mode{ScalarFieldFEM2D})
    E2_task=Threads.@spawn sum(integrate(abs2(m.field.E),m.field.dΩ));
    E4_task=Threads.@spawn sum(integrate(abs2(abs2(m.field.E)),m.field.dΩ));
    E2=fetch(E2_task);
    E4=fetch(E4_task);
    return E2^2/E4;
end

"""
    Aeff(m::Mode{VectorFieldFEM2D},n0::Union{Real,Function}=0)
"""
function Aeff(m::Mode{VectorFieldFEM2D},n0::Union{Real,Function}=0)
    if isa(n0,Function)
        if (first(methods(n0)).nargs!=2)
            throw(DomainError(n0, "The refractive index function n0 must have 1 argument (x,y)"));
        end
        n02=x->real(n0(x))^2;
        E2_task=Threads.@spawn sum(integrate(m.field.Ex*conj(m.field.Hy)-m.Ey*conj(m.field.Hx),m.field.dΩ));
        E4_1_task=Threads.@spawn sum(integrate(n02*abs2(abs2(m.field.Ex)+abs2(m.field.Ey)+abs2(m.field.Ez)),m.field.dΩ));
        E4_2_task=Threads.@spawn sum(integrate(n02*abs2(m.field.Ex*m.field.Ex+m.field.Ey*m.field.Ey+m.field.Ez*m.field.Ez),m.field.dΩ));
        E2=fetch(E2_task);
        E4_1=fetch(E4_1_task);
        E4_2=fetch(E4_2_task);
        return real(E2)^2/E4_1*mu0/eps0,real(E2)^2/E4_2*mu0/eps0;
    else
        E2_task=Threads.@spawn sum(integrate(m.field.Ex*conj(m.field.Hy)-m.field.Ey*conj(m.field.Hx),m.field.dΩ));
        E4_1_task=Threads.@spawn sum(integrate(abs2(abs2(m.field.Ex)+abs2(m.field.Ey)+abs2(m.field.Ez)),m.field.dΩ));
        E4_2_task=Threads.@spawn sum(integrate(abs2(m.field.Ex*m.field.Ex+m.field.Ey*m.field.Ey+m.field.Ez*m.field.Ez),m.field.dΩ));
        E2=fetch(E2_task);
        E4_1=fetch(E4_1_task);
        E4_2=fetch(E4_2_task);
        if (n0==0)
            return real(E2)^2/E4_1*mu0/eps0/(real(m.neff))^2,real(E2)^2/E4_2*mu0/eps0/(real(m.neff))^2;
        else
            return real(E2)^2/E4_1*mu0/eps0/real(n0)^2,real(E2)^2/E4_2*mu0/eps0/real(n0)^2;
        end
    end
end

############################# Normalization #############################
"""
    normalize!(m::Mode{ScalarFieldFunction1D};unitIntegral::Bool=true,kwargs...)
"""
function normalize!(m::Mode{ScalarFieldFunction1D};unitIntegral::Bool=true,kwargs...)
    f=x->2*pi*x[1]*abs2(m.field.E(x))
    integral=integrate1D(f;kwargs...)[1]
    if (m.field.nu!=0)
        integral=integral*0.5;
    end
    E0=deepcopy(m.field.E)
    if unitIntegral
        m.field.E=x->E0(x)/sqrt(integral);
    else
        m.field.E=x->E0(x)*sqrt(2*mu0*c/real(m.neff)/integral);
    end
    m
end

"""
    normalize!(m::Mode{ScalarFieldFunction2D};unitIntegral::Bool=true,kwargs...)
"""
function normalize!(m::Mode{ScalarFieldFunction2D};unitIntegral::Bool=true,kwargs...)
    f=x->abs2(m.field.E(x));
    integral=integrate2D(f;kwargs...)[1]
    E0=deepcopy(m.field.E)
    if unitIntegral
        m.field.E=x->E0(x)/sqrt(integral);
    else
        m.field.E=x->E0(x)*sqrt(2*mu0*c/real(m.neff)/integral);
    end
    m
end

"""
    normalize!(m::Mode{VectorFieldFunction2D};kwargs...)
"""
function normalize!(m::Mode{VectorFieldFunction2D};kwargs...)
    f=x->(m.field.Ex(x)*conj(m.field.Hy(x))-m.field.Ey(x)*conj(m.field.Hx(x)))*0.5;
    integral=abs(integrate2D(f;kwargs...)[1])
    Ex0=deepcopy(m.field.Ex)
    Ey0=deepcopy(m.field.Ey)
    Ez0=deepcopy(m.field.Ez)
    Hx0=deepcopy(m.field.Hx)
    Hy0=deepcopy(m.field.Hy)
    Hz0=deepcopy(m.field.Hz)
    m.field.Ex=x->Ex0(x)/sqrt(integral);
    m.field.Ey=x->Ey0(x)/sqrt(integral);
    m.field.Ez=x->Ez0(x)/sqrt(integral);
    m.field.Hx=x->Hx0(x)/sqrt(integral);
    m.field.Hy=x->Hy0(x)/sqrt(integral);
    m.field.Hz=x->Hz0(x)/sqrt(integral);
    m
end

"""
    normalize!(m::Mode{ScalarFieldMatrix1D};unitIntegral::Bool=true)
"""
function normalize!(m::Mode{ScalarFieldMatrix1D};unitIntegral::Bool=true)
    integral=trapz(m.field.r,2*pi*m.field.r.*(abs2.(m.field.E)));
    if (m.field.nu!=0)
        integral=integral*0.5;
    end
    if unitIntegral
        (m.field.E).=(m.field.E)/sqrt(integral);
    else
        (m.field.E).=(m.field.E)*sqrt(2*mu0*c/real(m.neff)/integral);
    end
    m
end

"""
    normalize!(m::Mode{ScalarFieldMatrix2D};unitIntegral::Bool=true)
"""
function normalize!(m::Mode{ScalarFieldMatrix2D};unitIntegral::Bool=true)
    integral=trapz((m.field.x,m.field.y),(abs2.(m.field.E)));
    if unitIntegral
        (m.field.E).=(m.field.E)/sqrt(integral);
    else
        (m.field.E).=(m.field.E)*sqrt(2*mu0*c/real(m.neff)/integral);
    end
    m
end

"""
    normalize!(m::Mode{VectorFieldMatrix2D})
"""
function normalize!(m::Mode{VectorFieldMatrix2D})
    integral=abs(trapz((m.field.x,m.field.y),m.field.Ex.*conj(m.field.Hy)-m.field.Ey.*conj(m.field.Hx))*0.5);
    (m.field.Ex).=(m.field.Ex)/sqrt(integral);
    (m.field.Ey).=(m.field.Ey)/sqrt(integral);
    (m.field.Ez).=(m.field.Ez)/sqrt(integral);
    (m.field.Hx).=(m.field.Hx)/sqrt(integral);
    (m.field.Hy).=(m.field.Hy)/sqrt(integral);
    (m.field.Hz).=(m.field.Hz)/sqrt(integral);
    m
end

"""
    normalize!(m::Mode{ScalarFieldFEM1D};unitIntegral::Bool=true)
"""
function normalize!(m::Mode{ScalarFieldFEM1D};unitIntegral::Bool=true)
    r=x->x[1]
    integral=sum(integrate(abs2(m.field.E)*2*pi*r,m.field.dΩ));
    if (m.field.nu!=0)
        integral=integral*0.5;
    end
    if unitIntegral
        (m.field.E)=(m.field.E)/sqrt(integral);
    else
        (m.field.E)=(m.field.E)*sqrt(2*mu0*c/real(m.neff)/integral);
    end
    m
end


"""
    normalize!(m::Mode{ScalarFieldFEM2D};unitIntegral::Bool=true)
"""
function normalize!(m::Mode{ScalarFieldFEM2D};unitIntegral::Bool=true)
    integral=sum(integrate(abs2(m.field.E),m.field.dΩ));
    if unitIntegral
        m.field.E=m.field.E/sqrt(integral);
    else
        m.field.E=m.field.E*sqrt(2*mu0*c/real(m.neff)/integral);
    end
    m
end

"""
    normalize!(m::Mode{VectorFieldFEM2D})
"""
function normalize!(m::Mode{VectorFieldFEM2D})
    integral=0.5*abs(sum(integrate(m.field.Ex*conj(m.field.Hy)-m.field.Ey*conj(m.field.Hx),m.field.dΩ)));
    m.field.Ex=m.field.Ex/sqrt(integral);
    m.field.Ey=m.field.Ey/sqrt(integral);
    m.field.Ez=m.field.Ez/sqrt(integral);
    m.field.Hx=m.field.Hx/sqrt(integral);
    m.field.Hy=m.field.Hy/sqrt(integral);
    m.field.Hz=m.field.Hz/sqrt(integral);
    m
end

############################# Overlap #############################
"""
    overlap(m1::ScalarFieldFunction1D,m2::ScalarFieldFunction1D;kwargs...)
"""
function overlap(m1::ScalarFieldFunction1D,m2::ScalarFieldFunction1D;kwargs...)
    if m1.nu!=m2.nu
        return 0.0;
    end
        f=x->m1.E(x)*conj(m2.E(x))*x[1]
        return integrate1D(f;kwargs...)[1]*(1+(m1.nu==0))*pi;
end

function overlap(m1::ScalarFieldFunction2D,m2::ScalarFieldFunction2D;kwargs...)
    f=x->m1.E(x)*conj(m2.E(x))
    return integrate2D(f;kwargs...)[1];
end

function overlap(f1::VectorFieldFunction2D,f2::VectorFieldFunction2D;kwargs...)
    f=x->(f1.Ex(x)*conj(f2.Hy(x))-f1.Ey(x)*conj(f2.Hx(x)))*0.5;
    return integrate2D(f;kwargs...)[1];
end

"""
    overlap(f1::Union{ScalarField,ScalarMode},f2::Union{ScalarField,ScalarMode})
"""
function overlap(f1::ScalarFieldMatrix2D,f2::ScalarFieldMatrix2D)
    if ((f1.x==f2.x) && (f1.y==f2.y))
        return trapz((f1.x,f1.y),f1.E.*(conj.(f2.E)));
    else
        x=unique(sort(vcat(f1.x,f2.x)));
        y=unique(sort(vcat(f1.y,f2.y)));
        interp=LinearInterpolation((f1.x,f1.y),f1.E,extrapolation_bc=0);
        E1=interp.(x,y');
        interp=LinearInterpolation((f2.x,f2.y),f2.E,extrapolation_bc=0);
        E2=interp.(x,y');
        return trapz((y,x),E1.*(conj.(E2)));
    end
end

"""
    overlap(f1::Union{VectorField,VectorMode},f2::Union{VectorField,VectorMode})
"""
function overlap(f1::VectorFieldMatrix2D,f2::VectorFieldMatrix2D)
    if ((f1.x==f2.x) && (f1.y==f2.y))
        return trapz((f2.x,f2.y),f1.Ex.*conj.(f2.Hy)-f1.Ey.*conj.(f2.Hx))*0.5;
    else
        x=unique(sort(vcat(f1.x,f2.x)));
        y=unique(sort(vcat(f1.y,f2.y)));
        interp=LinearInterpolation((f1.x,f1.y),f1.Ex,extrapolation_bc=0);
        Ex1=interp.(x,y');
        interp=LinearInterpolation((f1.x,f1.y),f1.Ey,extrapolation_bc=0);
        Ey1=interp.(x,y');
        interp=LinearInterpolation((f2.x,f2.y),f2.Hx,extrapolation_bc=0);
        Hx2=interp.(x,y');
        interp=LinearInterpolation((f2.x,f2.y),f2.Hy,extrapolation_bc=0);
        Hy2=interp.(x,y');
        return trapz((y,x),Ex1.*conj.(Hy2)-Ey1.*conj.(Hx2))*0.5;
    end
end

"""
    overlap(m1::ScalarMode1D,m2::ScalarMode1D)
"""
function overlap(m1::ScalarFieldMatrix1D,m2::ScalarFieldMatrix1D)
    if m1.nu!=m2.nu
        return 0.0;
    end
    if (m1.r==m2.r)
        return trapz(m1.r,m1.r.*m1.E.*conj(m2.E))*(1+(m1.nu==0))*pi;
    else
        r=unique(sort(vcat(m1.r,m2.r)));
        interp1=LinearInterpolation(m1.r,m1.E,extrapolation_bc=0);
        E1=interp1.(r);
        interp2=LinearInterpolation(m2.r,m2.E,extrapolation_bc=0);
        E2=interp2.(r);
        return trapz(r,r.*E1.*conj(E2))*(1+(m1.nu==0))*pi;
    end
end

"""
    overlap(m1::ScalarModeFEM1D,m2::ScalarModeFEM1D)
"""
function overlap(m1::ScalarFieldFEM1D,m2::ScalarFieldFEM1D)
    if (get_triangulation(m1.E)!=get_triangulation(m2.E))
        throw(ArgumentError("The modes must have the same triangulation Ω"));
    end
    if m1.nu!=m2.nu
        return 0.0;
    else
        r=x->x[1]
        integral=sum(integrate(m1.E*conj(m2.E)*r,m2.dΩ));
        return integral*(1+(m1.nu==0))*pi;
    end
end

"""
    overlap(m1::Union{ScalarFieldFEM,ScalarModeFEM},m2::Union{ScalarFieldFEM,ScalarModeFEM})
"""
function overlap(m1::ScalarFieldFEM2D,m2::ScalarFieldFEM2D)
    if (get_triangulation(m1.E)!=get_triangulation(m2.E))
        throw(ArgumentError("The modes must have the same triangulation Ω"));
    else
        integral=sum(integrate(m1.E*conj(m2.E),m2.dΩ));
        return integral
    end
end

"""
    overlap(m1::Union{VectorFieldFEM,VectorModeFEM},m2::Union{VectorFieldFEM,VectorModeFEM})
"""
function overlap(m1::VectorFieldFEM2D,m2::VectorFieldFEM2D)
    if (get_triangulation(m1.Ex)!=get_triangulation(m2.Ex))
        throw(ArgumentError("The modes must have the same triangulation Ω"));
    else
        integral=0.5*sum(integrate(m1.Ex*conj(m2.Hy)-m1.Ey*conj(m2.Hx),m1.dΩ));
        return integral
    end
end

function overlap(m1::ScalarFieldFEM1D,m2::ScalarFieldFunction1D)
    if m1.nu!=m2.nu
        return 0.0;
    else
        r2=x->x[1]*conj(m2.E(x))
        integral=sum(integrate(m1.E*r2,m1.dΩ));
        return integral*(1+(m1.nu==0))*pi;
    end
end

function overlap(m1::ScalarFieldFunction1D,m2::ScalarFieldFEM1D)
    return conj(overlap(m2,m1))
end

function overlap(m1::ScalarFieldFEM2D,m2::ScalarFieldFunction2D)
    f=x->conj(m2.E(x))
    integral=sum(integrate(m1.E*f,m1.dΩ));
    return integral
end

function overlap(m1::ScalarFieldFunction2D,m2::ScalarFieldFEM2D)
    return conj(overlap(m2,m1))
end

function overlap(m1::VectorFieldFEM2D,m2::VectorFieldFunction2D)
    f1=x->conj(m2.Hy(x))
    f2=x->conj(m2.Hx(x))
    integral=0.5*sum(integrate(m1.Ex*f1-m1.Ey*f2,m1.dΩ));
    return integral
end

function overlap(m1::VectorFieldFunction2D,m2::VectorFieldFEM2D)
    return conj(overlap(m2,m1))
end

overlap(m::Mode,f::Field) = overlap(m.field,f)
overlap(f::Field,m::Mode) = overlap(f,m.field)
overlap(m1::Mode,m2::Mode) = overlap(m1.field,m2.field)
overlap(m1::Mode{VectorFieldFunction2D},m2::Mode{VectorFieldFunction2D};kwargs...) = overlap(m1.field,m2.field;kwargs...)
overlap(m1::Mode{ScalarFieldFunction2D},m2::Mode{ScalarFieldFunction2D};kwargs...) = overlap(m1.field,m2.field;kwargs...)
overlap(m1::Mode{ScalarFieldFunction1D},m2::Mode{ScalarFieldFunction1D};kwargs...) = overlap(m1.field,m2.field;kwargs...)

############################# Nonlinear Coefficient #############################
"""
    nonLinearCoefficient(m::Mode{ScalarFieldFunction1D},n2::Union{Real,Function};kwargs...)
"""
function nonLinearCoefficient(m::Mode{ScalarFieldFunction1D},n2::Union{Real,Function};kwargs...)
    if (isa(n2,Real))
        return n2*2*pi/m.lambda/Aeff(m;kwargs...);
    else
        if (first(methods(n2)).nargs!=2)
            throw(DomainError(n2, "The non-linear index function must have 1 argument"));
        end
        f2=x->2*pi*x[1]*abs2.(m.field.E(x));
        f4=x->2*pi*x[1]*(abs2(m.field.E(x)))^2*n2(x);
        E2=integrate1D(f2;kwargs...)[1]
        E4=integrate1D(f4;kwargs...)[1]
        if (m.field.nu!=0)
            E2=E2*0.5;
            E4=E4*3.0/8.0;
        end
        k0=2*pi/m.lambda;
        return E4/(E2^2)*k0;
    end
end

"""
    nonLinearCoefficient(m::Mode{ScalarFieldFunction2D},n2::Union{Real,Function};kwargs...)
"""
function nonLinearCoefficient(m::Mode{ScalarFieldFunction2D},n2::Union{Real,Function};kwargs...)
    if (isa(n2,Real))
        return n2*2*pi/m.lambda/Aeff(m;kwargs...);
    else
        if (first(methods(n2)).nargs!=2)
            throw(DomainError(n2, "The non-linear index function must have 1 argument"));
        end
        f2=x->abs2(m.field.E(x));
        f4=x->n2(x)*(abs2(m.field.E(x)))^2;
        E2=integrate2D(f2;kwargs...)[1]
        E4=integrate2D(f4;kwargs...)[1]
        k0=2*pi/m.lambda;
        return E4/(E2^2)*k0;
    end
end

"""
    nonLinearCoefficient(m::Mode{VectorFieldFunction2D},n2::Union{Real,Function},n0::Union{Real,Function}=0;kwargs...)
"""
function nonLinearCoefficient(m::Mode{VectorFieldFunction2D},n2::Union{Real,Function},n0::Union{Real,Function}=0;kwargs...)
    if (isa(n2,Real))
        return n2*2*pi/m.lambda./Aeff(m,n0;kwargs...);
    else
        if (first(methods(n2)).nargs!=2)
            throw(DomainError(n2, "The non-linear index function must have 2 argument"));
        end
        f2=x->m.Ex(x)*conj(m.Hy(x))-m.Ey(x)*conj(m.Hx(x))
        E2=real(integrate2D(f2;kwargs...)[1])
        k0=2*pi/m.lambda;
        if (isa(n0,Real))
            if (n0==0)
                n0n2=x->real(m.neff)^2*n2(x);
            else
                n0n2=x->real(n0)^2*n2(x);
            end
        else
            n0n2=x->real(n0(x))^2*n2(x);
        end
        f4_1=x->abs2(abs2(m.field.Ex(x))+abs2(m.field.Ey(x))+abs2(m.field.Ez(x)))*n0n2(x)
        f4_2=x->abs2(m.field.Ex(x)^2+m.field.Ey(x)^2+m.field.Ez(x)^2)*n2(x)
        E4_1=integrate2D(f4_1;kwargs...)[1]
        E4_2=integrate2D(f4_2;kwargs...)[1]
        return 1.0/(E2^2/E4_1*mu0/eps0)*k0,1.0/(E2^2/E4_2*mu0/eps0)*k0;
    end
end

"""
    nonLinearCoefficient(m::Mode{ScalarFieldMatrix1D},n2::Union{Real,Function})
"""
function nonLinearCoefficient(m::Mode{ScalarFieldMatrix1D},n2::Union{Real,Function})
    if (isa(n2,Real))
        return n2*2*pi/m.lambda./Aeff(m);
    else
        if (first(methods(n2)).nargs!=2)
            throw(DomainError(n2, "The non-linear index function must have 1 argument"));
        end
        E2=trapz(m.field.r,2*pi*m.field.r.*abs2.(m.field.E));
        E4=trapz(m.field.r,2*pi*(m.r.*(abs2.(m.field.E)).^2).*n2.(m.field.r));
        if (m.field.nu!=0)
            E2=E2*0.5;
            E4=E4*3.0/8.0;
        end
        k0=2*pi/m.lambda;
        return E4/(E2^2)*k0;
    end
end

"""
    nonLinearCoefficient(m::Mode{ScalarFieldMatrix2D},n2::Union{Real,Function})
"""
function nonLinearCoefficient(m::Mode{ScalarFieldMatrix2D},n2::Union{Real,Function})
    if (isa(n2,Real))
        return n2*2*pi/m.lambda./Aeff(m);
    else
        if (first(methods(n2)).nargs!=2)
            throw(DomainError(n2, "The non-linear index function must have 1 argument"));
        end
        E2=trapz((m.field.x,m.field.y),abs2.(m.field.E));
        E4=trapz((m.field.x,m.field.y),((abs2.(m.field.E)).^2).*n2.(Point.(m.field.x,m.field.y')));
        k0=2*pi/m.lambda;
        return E4/(E2^2)*k0;
    end
end

"""
    nonLinearCoefficient(m::Mode{VectorFieldMatrix2D},n2::Union{Real,Function},n0::Union{Real,Function}=0)
"""
function nonLinearCoefficient(m::Mode{VectorFieldMatrix2D},n2::Union{Real,Function},n0::Union{Real,Function}=0)
    if (isa(n2,Real))
        return n2*2*pi/m.lambda./Aeff(m,n0);
    else
        if (first(methods(n2)).nargs!=2)
            throw(DomainError(n2, "The non-linear index function must have 2 argument"));
        end
        E2=real(trapz((m.x,m.y),m.Ex.*conj(m.Hy)-m.Ey.*conj(m.Hx)));
        k0=2*pi/m.lambda;
        if (isa(n0,Real))
            if (n0==0)
                n0n2=x->real(m.neff)^2*n2(x);
            else
                n0n2=x->real(n0)^2*n2(x);
            end
        else
            n0n2=x->real(n0(x))^2*n2(x);
        end
        E4_1=trapz((m.field.x,m.field.y),(abs2.(abs2.(m.field.Ex)+abs2.(m.field.Ey)+abs2.(m.field.Ez))).*n0n2.(Point.(m.field.x,m.field.y')));
        E4_2=trapz((m.field.x,m.field.y),(abs2.(m.field.Ex.^2+m.field.Ey.^2+m.field.Ez.^2)).*n2.(Point.(m.field.x,m.field.y')));
        return 1.0/(E2^2/E4_1*mu0/eps0)*k0,1.0/(E2^2/E4_2*mu0/eps0)*k0;
    end
end

"""
    nonLinearCoefficient(m::Mode{ScalarFieldFEM1D},n2::Union{Real,Function})
"""
function nonLinearCoefficient(m::Mode{ScalarFieldFEM1D},n2::Union{Real,Function})
    if (isa(n2,Real))
        return n2*2*pi/m.lambda./Aeff(m);
    else
        if (first(methods(n2)).nargs!=2)
            throw(DomainError(n2, "The non-linear index function must have 1 argument"));
        end
        r=x->x[1];
        n2bis=x->n2(x);
        E2_task=Threads.@spawn sum(integrate(abs2(m.field.E)*2*pi*r,m.field.dΩ));
        E4_task=Threads.@spawn sum(integrate(n2bis*abs2(abs2(m.field.E))*2*pi*r,m.field.dΩ));
        E2=fetch(E2_task);
        E4=fetch(E4_task);
        if (m.field.nu!=0)
            E2=E2*0.5;
            E4=E4*3.0/8.0;
        end
        k0=2*pi/m.lambda;
        return E4/(E2^2)*k0;
    end
end

"""
    nonLinearCoefficient(m::Mode{ScalarFieldFEM2D},n2::Union{Real,Function})
"""
function nonLinearCoefficient(m::Mode{ScalarFieldFEM2D},n2::Union{Real,Function})
    if (isa(n2,Real))
        return n2*2*pi/m.lambda./Aeff(m);
    else
        if (first(methods(n2)).nargs!=2)
            throw(DomainError(n2, "The non-linear index function must have 1 argument"));
        end
        n2bis=x->n2(x);
        E2_task=Threads.@spawn sum(integrate(abs2(m.field.E),m.field.dΩ));
        E4_task=Threads.@spawn sum(integrate(n2bis*abs2(abs2(m.field.E)),m.field.dΩ));
        E2=fetch(E2_task);
        E4=fetch(E4_task);
        k0=2*pi/m.lambda;
        return E4/(E2^2)*k0;
    end
end

"""
    nonLinearCoefficient(m::Mode{VectorFieldFEM2D},n2::Union{Real,Function},n0::Union{Real,Function}=0)
"""
function nonLinearCoefficient(m::Mode{VectorFieldFEM2D},n2::Union{Real,Function},n0::Union{Real,Function}=0)
    if (isa(n2,Real))
        return n2*2*pi/m.lambda./Aeff(m,n0);
    else
        if (first(methods(n2)).nargs!=2)
            throw(DomainError(n2, "The non-linear index n2 function must have 1 argument"));
        end
        k0=2*pi/m.lambda;
        if (isa(n0,Real))
            if (n0==0)
                n02n2=x->n2(x)*real(m.neff)^2;
            else
                n02n2=x->n2(x)*n0^2;
            end
        else
            if (first(methods(n0)).nargs!=2)
                throw(DomainError(n0, "The refractive index n0 function must have 1 argument"));
            end
            n02n2=x->n2(x)*n0(x)^2;
        end
        E2_task=Threads.@spawn real(sum(integrate(m.field.Ex*conj(m.field.Hy)-m.field.Ey*conj(m.field.Hx),m.field.dΩ)));
        E4_1_task=Threads.@spawn sum(integrate(n02n2*abs2(abs2(m.field.Ex)+abs2(m.field.Ey)+abs2(m.field.Ez)),m.field.dΩ));
        E4_2_task=Threads.@spawn sum(integrate(n02n2*abs2(m.field.Ex*m.field.Ex+m.field.Ey*m.field.Ey+m.field.Ez*m.field.Ez),m.field.dΩ));
        E2=fetch(E2_task);
        E4_1=fetch(E4_1_task);
        E4_2=fetch(E4_2_task);
        return 1.0/(E2^2/E4_1*mu0/eps0)*k0,1.0/(E2^2/E4_2*mu0/eps0)*k0;
    end
end

############################# MFD #############################
"""
    MFD(m::ScalarFieldFunction1D)
"""
function MFD(m::ScalarFieldFunction1D)
    f=x->m.E(x)-m.E(Point(0.0,))/ℯ
    return 2*abs(find_zero(f,0.0))
end

"""
    MFD(m::ScalarFieldFunction2D,theta::Real=0)
"""
function MFD(m::ScalarFieldFunction2D,theta::Real=0)
    f=x->m.E(Point(x*cosd(theta),x*sind(theta)))-m.E(Point(0.0,0.0))/ℯ
    return 2*abs(find_zero(f,0.0))
end

"""
    MFD(m::VectorFieldFunction2D,theta::Real=0)
"""
function MFD(m::VectorFieldFunction2D,theta::Real=0)
    Pz=x->real(m.Ex(x)*conj(m.Hy(x))-m.Ey(x)*conj(m.Hx(x)))
    f=x->Pz(Point(x*cosd(theta),x*sind(theta)))-Pz(Point(0.0,0.0))/ℯ^2
    return 2*abs(find_zero(f,0.0))
end

"""
    MFD(m::ScalarFieldMatrix1D)
"""
function MFD(m::ScalarFieldMatrix1D)
    pos=argmin(abs.(abs.(m.E).-(maximum(abs.(m.E))/ℯ)));
    return 2*m.r[pos];
end

"""
    MFD(m::ScalarFieldMatrix2D,theta::Real=0)
"""
function MFD(m::ScalarFieldMatrix2D,theta::Real=0)
    rmax=sqrt(maximum(m.x)^2+maximum(m.y)^2);
    dr=min(minimum(diff(m.x)),minimum(diff(m.y)));
    nb=Integer(2*ceil(rmax/dr)+1);
    r=collect(LinRange(-rmax,rmax,nb));
    x=r*cosd(theta);
    y=r*sind(theta);
    interp=LinearInterpolation((m.x,m.y),m.E,extrapolation_bc=0);
    E=interp.(x,y);
    posmax=argmax(abs.(E));
    pos=argmin(abs.(abs.(E).-(maximum(abs.(E))/ℯ)));
    return abs(2*(r[pos]-r[posmax]));
end

"""
    MFD(m::VectorFieldMatrix2D,theta::Real=0)
"""
function MFD(m::VectorFieldMatrix2D,theta::Real=0)
    rmax=sqrt(maximum(m.x)^2+maximum(m.y)^2);
    dr=min(minimum(diff(m.x)),minimum(diff(m.y)));
    nb=Integer(2*ceil(rmax/dr)+1);
    r=collect(LinRange(-rmax,rmax,nb));
    x=r*cosd(theta);
    y=r*sind(theta);
    Pz=real.(m.Ex.*conj(m.Hy)-m.Ey.*conj(m.Hx));
    interp=LinearInterpolation((m.x,m.y),Pz,extrapolation_bc=0);
    E2=interp.(x,y);
    posmax=argmax(abs.(E2));
    pos=argmin(abs.(abs.(E2).-(maximum(abs.(E2))/(ℯ^2))));
    return abs(2*(r[pos]-r[posmax]));
end

"""
    MFD(m::ScalarFieldFEM1D)
"""
function MFD(m::ScalarFieldFEM1D)
    Ω=get_triangulation(m.E)
    r=getindex.(Ω.grid.node_coords,1)
    rmax=maximum(r);
    r2=LinRange(0.0,rmax,100*length(r))
    E=getValue(m.E,Point.(r2,))
    pos=argmin(abs.(abs.(E).-(maximum(abs.(E))/ℯ)));
    return 2*r2[pos];
end

"""
    MFD(m::ScalarFieldFEM2D,theta::Real=0)
"""
function MFD(m::ScalarFieldFEM2D,theta::Real=0)
    Ω=get_triangulation(m.E)
    x=getindex.(Ω.grid.node_coordinates,1)
    y=getindex.(Ω.grid.node_coordinates,2)
    rmin,rmax=extrema(x*cosd(theta)+y*sind(theta))
    r=LinRange(rmin,rmax,10000)
    E=getValue(m.E,Point.(r*cosd(theta),r*sind(theta)))
    posmax=argmax(abs.(E));
    pos=argmin(abs.(abs.(E).-(maximum(abs.(E))/ℯ)));
    return abs(2*(r[pos]-r[posmax]));
end

"""
    MFD(m::VectorFieldFEM2D,theta::Real=0)
"""
function MFD(m::VectorFieldFEM2D,theta::Real=0)
    Ω=get_triangulation(m.Ex)
    x=getindex.(Ω.grid.node_coordinates,1)
    y=getindex.(Ω.grid.node_coordinates,2)
    rmax=maximum(x*cosd(theta)+y*sind(theta))
    rmin=minimum(x*cosd(theta)+y*sind(theta))
    r=LinRange(rmin,rmax,10000)
    Pzf=real(m.Ex*conj(m.Hy)-m.Ey*conj(m.Hx));
    Pz=getValue(Pzf,Point.(r*cosd(theta),r*sind(theta)))
    posmax=argmax(abs.(Pz));
    pos=argmin(abs.(abs.(Pz).-(maximum(abs.(Pz))/(ℯ^2))));
    return abs(2*(r[pos]-r[posmax]));
end

"""
    MFD(m::Mode{<:ScalarField1D})
"""
MFD(m::Mode{<:ScalarField1D}) = MFD(m.field)
"""
    MFD(m::Union{Mode{<:ScalarField2D},Mode{<:VectorField2D}},theta::Real=0)
"""
MFD(m::Union{Mode{<:ScalarField2D},Mode{<:VectorField2D}},theta::Real=0) = MFD(m.field,theta)

####################writevtk############################
"""
    Gridap.:writevtk(name::String,f::Union{ScalarFieldFEM2D,ScalarFieldFEM1D})
"""
function Gridap.:writevtk(name::String,f::Union{ScalarFieldFEM2D,ScalarFieldFEM1D})
    Ω=get_triangulation(f.E)
    Gridap.writevtk(Ω,name,cellfields=["real(E)"=>real(f.E),"imag(E)"=>imag(f.E),"abs2(E)"=>abs2(f.E)]);
    return nothing;
end

"""
    Gridap.:writevtk(name::String,f::VectorFieldFEM2D)
"""
function Gridap.:writevtk(name::String,f::VectorFieldFEM2D)
    Ω=get_triangulation(f.Ex)
    Gridap.writevtk(Ω,name,cellfields=["real(Ex)"=>real(f.Ex),"imag(Ex)"=>imag(f.Ex),"real(Ey)"=>real(f.Ey),"imag(Ey)"=>imag(f.Ey),"real(Ez)"=>real(f.Ez),"imag(Ez)"=>imag(f.Ez),"real(Hx)"=>real(f.Hx),"imag(Hx)"=>imag(f.Hx),"real(Hy)"=>real(f.Hy),"imag(Hy)"=>imag(f.Hy),"real(Hz)"=>real(f.Hz),"imag(Hz)"=>imag(f.Hz)]);
    return nothing;
end

"""
    Gridap.:writevtk(name::String,m::Union{Mode{ScalarFieldFEM2D},Mode{ScalarFieldFEM1D},Mode{VectorFieldFEM2D}})
"""
Gridap.:writevtk(name::String,m::Union{Mode{ScalarFieldFEM2D},Mode{ScalarFieldFEM1D},Mode{VectorFieldFEM2D}}) = Gridap.writevtk(name,m.field)


############################# Losses #############################
"""
    losses(m::Mode)
"""
function losses(m::Mode)
    return 4*pi/m.lambda*imag(m.neff)*10/log(10)*1000
end

############################# Triangulation #############################
Gridap.:get_triangulation(f::Union{ScalarFieldFEM2D,ScalarFieldFEM1D}) = Gridap.get_triangulation(f.E)
Gridap.:get_triangulation(f::VectorFieldFEM2D) = Gridap.get_triangulation(f.Ex)
"""
    Gridap.:get_triangulation(m::Union{Mode{ScalarFieldFEM2D},Mode{ScalarFieldFEM1D},Mode{VectorFieldFEM2D}})
"""
Gridap.:get_triangulation(m::Union{Mode{ScalarFieldFEM2D},Mode{ScalarFieldFEM1D},Mode{VectorFieldFEM2D}}) = Gridap.get_triangulation(m.field)
