############################# Fields #############################

"""
Abstract structure to describe an electromagnetic field

`abstract type Field end`
"""
abstract type Field end

"""
Structure describing a 2D scalar field in cartesian coordinates

- x :: `Vector{Float64}` - Horizontal coordinate
- y :: `Vector{Float64}` - Vertical coordinate
- E :: `Matrix{ComplexF64}` - Electric field
"""
struct ScalarField <: Field
    x::Vector{Float64}
    y::Vector{Float64}
    E::Matrix{ComplexF64}
    ScalarField(x,y,E) = new(Float64.(x),Float64.(y),ComplexF64.(E))
end

function Base.:+(f1::ScalarField, f2::ScalarField)
    if ((f1.x==f2.x) && (f1.y==f2.y))
        ScalarField(f1.x,f1.y,f1.E+f2.E)
    else
        throw(ArgumentError("Both fields must have the same coordinates"))
    end
end

function Base.:-(f1::ScalarField, f2::ScalarField)
    if ((f1.x==f2.x) && (f1.y==f2.y))
        ScalarField(f1.x,f1.y,f1.E-f2.E)
    else
        throw(ArgumentError("Both fields must have the same coordinates"))
    end
end

Base.:*(k::Number, f::ScalarField) = ScalarField(f.x,f.y,k*f.E);
Base.:*(f::ScalarField,k::Number) = ScalarField(f.x,f.y,k*f.E);
Base.:/(f::ScalarField,k::Number) = ScalarField(f.x,f.y,f.E/k);
Base.:+(f::ScalarField) = ScalarField(f.x,f.y,f.E);
Base.:-(f::ScalarField) = ScalarField(f.x,f.y,-f.E);
Base.:real(f::ScalarField) = ScalarField(f.x,f.y,real(f.E));
Base.:imag(f::ScalarField) = ScalarField(f.x,f.y,imag(f.E));
Base.:conj(f::ScalarField) = ScalarField(f.x,f.y,conj(f.E));

Base.show(io::IO, f::ScalarField) = print(io,"x = ",f.x[1]," : ",f.x[end],"\ny = ",f.y[1]," : ",f.y[end],"\n|E| ∈ [",minimum(abs.(f.E)),",",maximum(abs.(f.E)),"]");

"""
Structure describing a 2D vector field in cartesian coordinates

- x :: `Vector{Float64}` - Horizontal coordinate
- y :: `Vector{Float64}` - Vertical coordinate
- Ex :: `Matrix{ComplexF64}` - x-component of the electric field
- Ey :: `Matrix{ComplexF64}` - y-component of the electric field
- Ez :: `Matrix{ComplexF64}` - z-component of the electric field
- Hx :: `Matrix{ComplexF64}` - x-component of the magnetic field
- Hy :: `Matrix{ComplexF64}` - y-component of the magnetic field
- Hz :: `Matrix{ComplexF64}` - z-component of the magnetic field
"""
struct VectorField <: Field
    x::Vector{Float64}
    y::Vector{Float64}
    Ex::Matrix{ComplexF64}
    Ey::Matrix{ComplexF64}
    Ez::Matrix{ComplexF64}
    Hx::Matrix{ComplexF64}
    Hy::Matrix{ComplexF64}
    Hz::Matrix{ComplexF64}
    VectorField(x,y,Ex,Ey,Ez,Hx,Hy,Hz) = new(Float64.(x),Float64.(y),ComplexF64.(Ex),ComplexF64.(Ey),ComplexF64.(Ez),ComplexF64.(Hx),ComplexF64.(Hy),ComplexF64.(Hz))
end

function Base.:+(f1::VectorField, f2::VectorField)
    if ((f1.x==f2.x) && (f1.y==f2.y))
        VectorField(f1.x,f1.y,f1.Ex+f2.Ex,f1.Ey+f2.Ey,f1.Ez+f2.Ez,f1.Hx+f2.Hx,f1.Hy+f2.Hy,f1.Hz+f2.Hz)
    else
        throw(ArgumentError("Both fields must have the same coordinates"))
    end
end

function Base.:-(f1::VectorField, f2::VectorField)
    if ((f1.x==f2.x) && (f1.y==f2.y))
        VectorField(f1.x,f1.y,f1.Ex-f2.Ex,f1.Ey-f2.Ey,f1.Ez-f2.Ez,f1.Hx-f2.Hx,f1.Hy-f2.Hy,f1.Hz-f2.Hz)
    else
        throw(ArgumentError("Both fields must have the same coordinates"))
    end
end

Base.:*(k::Number, f::VectorField) = VectorField(f.x,f.y,k*f.Ex,k*f.Ey,k*f.Ez,k*f.Hx,k*f.Hy,k*f.Hz);
Base.:*(f::VectorField,k::Number) = VectorField(f.x,f.y,k*f.Ex,k*f.Ey,k*f.Ez,k*f.Hx,k*f.Hy,k*f.Hz);
Base.:/(f::VectorField,k::Number) = VectorField(f.x,f.y,f.Ex/k,f.Ey/k,f.Ez/k,f.Hx/k,f.Hy/k,f.Hz/k);
Base.:+(f::VectorField) = VectorField(f.x,f.y,f.Ex,f.Ey,f.Ez,f.Hx,f.Hy,f.Hz);
Base.:-(f::VectorField) = VectorField(f.x,f.y,-(f.Ex),-(f.Ey),-(f.Ez),-(f.Hx),-(f.Hy),-(f.Hz));
Base.:real(f::VectorField) = VectorField(f.x,f.y,real(f.Ex),real(f.Ey),real(f.Ez),real(f.Hx),real(f.Hy),real(f.Hz));
Base.:imag(f::VectorField) = VectorField(f.x,f.y,imag(f.Ex),imag(f.Ey),imag(f.Ez),imag(f.Hx),imag(f.Hy),imag(f.Hz));
Base.:conj(f::VectorField) = VectorField(f.x,f.y,conj(f.Ex),conj(f.Ey),conj(f.Ez),conj(f.Hx),conj(f.Hy),conj(f.Hz));

############################# Modes #############################

"""
Abstract structure to describe an optical fiber mode

`abstract type Mode end`
"""
abstract type Mode end

"""
Structure describing a scalar mode in an optical fiber with a revolution symmetry

- Name :: `String` - Name of the mode
- neff :: `Float64` - Effective index
- lambda :: `Float64` - the wavelength at which the mode was calculated
- nu :: `Int64` - Azimuthal number
- r :: `Vector{Float64}` - Radial coordinate
- E :: `Vector{Float64}` - Electric field
"""
struct ScalarMode1D <: Mode
    Name::String
    neff::Float64
    lambda::Float64
    nu::Int64
    r::Vector{Float64}
    E::Vector{Float64}
    ScalarMode1D(Name,neff,lambda,nu,r,E) = new(Name,neff,lambda,nu,Float64.(r),Float64.(E))
    ScalarMode1D(Name,neff,lambda,nu)= new(Name,neff,lambda,nu,[],[])
end

function Base.show(io::IO, ::MIME"text/plain",f::ScalarMode1D) 
    if isempty(f.r) 
        print(io,"Name = ",f.Name,"\nneff = ",f.neff,"\nlambda = ",f.lambda,"\nnu = ",f.nu);
    else
        print(io,"Name = ",f.Name,"\nneff = ",f.neff,"\nlambda = ",f.lambda,"\nnu = ",f.nu,"\nr ∈ [",f.r[1],",",f.r[end],"]\nE ∈ [",minimum(f.E),",",maximum(f.E),"]");
    end
end
Base.show(io::IO, f::ScalarMode1D) = isempty(f.r) ? print(io,"[\"",f.Name,"\",",f.neff,",",f.lambda,",",f.nu,"]") : print(io,"[\"",f.Name,"\",",f.neff,",",f.lambda,",",f.nu,",[",minimum(f.r),",",maximum(f.r),"],[",minimum(f.E),",",maximum(f.E),"]]");

"""
Structure describing a scalar mode in an optical fiber in cartesian coordinates

- Name :: `String` - Name of the mode
- neff :: `Float64` - Effective index
- lambda :: `Float64` - the wavelength at which the mode was calculated
- x :: `Vector{Float64}` - Horizontal coordinate
- y :: `Vector{Float64}` - Vertical coordinate
- E :: `Matrix{Float64}` - Electric field
"""
struct ScalarMode2D <: Mode
    Name::String
    neff::Float64
    lambda::Float64
    x::Vector{Float64}
    y::Vector{Float64}
    E::Matrix{Float64}
    ScalarMode2D(Name,neff,lambda,x,y,E)=new(Name,neff,lambda,Float64.(x),Float64.(y),Float64.(E))
    ScalarMode2D(Name,neff,lambda)=new(Name,neff,lambda,[],[],Matrix{Float64}(undef,0,0))
end

function Base.show(io::IO, ::MIME"text/plain",f::ScalarMode2D) 
    if isempty(f.x) 
        print(io,"Name = ",f.Name,"\nneff = ",f.neff,"\nlambda = ",f.lambda);
    else
        print(io,"Name = ",f.Name,"\nneff = ",f.neff,"\nlambda = ",f.lambda,"\nx ∈ [",f.x[1],",",f.x[end],"]\ny ∈ [",f.y[1],",",f.y[end],"]\nE ∈ [",minimum(f.E),",",maximum(f.E),"]");
    end
end
Base.show(io::IO, f::ScalarMode2D) = isempty(f.x) ? print(io,"[\"",f.Name,"\",",f.neff,",",f.lambda,"]") : print(io,"[\"",f.Name,"\",",f.neff,",",f.lambda,",[",minimum(f.x),",",maximum(f.x),"],[",minimum(f.y),",",maximum(f.y),"],[",minimum(f.E),",",maximum(f.E),"]]");

"""
Structure describing a vector mode in an optical fiber in cartesian coordinates.  

This structure assumes that the fiber is made of a non-lossy materials. In this case, the x- and y-components of the electric/magnetic fields can be chosen real and their z-components are then pure imaginary numbers.

- Name :: `String` - Name of the mode
- neff :: `Float64` - Effective index
- lambda :: `Float64` - the wavelength at which the mode was calculated
- x :: `Vector{Float64}` - Horizontal coordinate
- y :: `Vector{Float64}` - Vertical coordinate
- Ex :: `Matrix{Float64}` - Real part of the x-component of the electric field
- Ey :: `Matrix{Float64}` - Real part of the y-component of the electric field
- Ez :: `Matrix{Float64}` - Imaginary part of the z-component of the electric field
- Hx :: `Matrix{Float64}` - Real part of the x-component of the magnetic field
- Hy :: `Matrix{Float64}` - Real part of the y-component of the magnetic field
- Hz :: `Matrix{Float64}` - Imaginary part of the z-component of the magnetic field
"""
struct VectorMode <: Mode
    Name::String
    neff::Float64
    lambda::Float64
    x::Vector{Float64}
    y::Vector{Float64}
    Ex::Matrix{Float64}
    Ey::Matrix{Float64}
    Ez::Matrix{Float64}
    Hx::Matrix{Float64}
    Hy::Matrix{Float64}
    Hz::Matrix{Float64}
    VectorMode(Name,neff,lambda,x,y,Ex,Ey,Ez,Hx,Hy,Hz)=new(Name,neff,lambda,Float64.(x),Float64.(y),Float64.(Ex),Float64.(Ey),Float64.(Ez),Float64.(Hx),Float64.(Hy),Float64.(Hz))
    VectorMode(Name,neff,lambda)=new(Name,neff,lambda,[],[],Matrix{Float64}(undef,0,0),Matrix{Float64}(undef,0,0),Matrix{Float64}(undef,0,0),Matrix{Float64}(undef,0,0),Matrix{Float64}(undef,0,0),Matrix{Float64}(undef,0,0))
end

function Base.show(io::IO, ::MIME"text/plain",f::VectorMode) 
    if isempty(f.x) 
        print(io,"Name = ",f.Name,"\nneff = ",f.neff,"\nlambda = ",f.lambda);
    else
        print(io,"Name = ",f.Name,"\nneff = ",f.neff,"\nlambda = ",f.lambda,"\nx ∈ [",f.x[1],",",f.x[end],"]\ny ∈ [",f.y[1],",",f.y[end],"]\nEx ∈ [",minimum(f.Ex),",",maximum(f.Ex),"]\nEy ∈ [",minimum(f.Ey),",",maximum(f.Ey),"]","]\nEz ∈ [",minimum(f.Ez),",",maximum(f.Ez),"]","]\nHx ∈ [",minimum(f.Hx),",",maximum(f.Hx),"]","]\nHy ∈ [",minimum(f.Hy),",",maximum(f.Hy),"]","]\nHz ∈ [",minimum(f.Hz),",",maximum(f.Hz),"]","]");
    end
end
Base.show(io::IO, f::VectorMode) = isempty(f.x) ? print(io,"[\"",f.Name,"\",",f.neff,",",f.lambda,"]") : print(io,"[\"",f.Name,"\",",f.neff,",",f.lambda,",[",minimum(f.x),",",maximum(f.x),"],[",minimum(f.y),",",maximum(f.y),"],[",minimum(f.Ex),",",maximum(f.Ex),"],[",minimum(f.Ey),",",maximum(f.Ey),"],[",minimum(f.Ez),",",maximum(f.Ez),"],[",minimum(f.Hx),",",maximum(f.Hx),"],[",minimum(f.Hy),",",maximum(f.Hy),"],[",minimum(f.Hz),",",maximum(f.Hz),"]]");

mutable struct ScalarModeFEM <: Mode
    Name::String
    neff::ComplexF64
    lambda::Float64
    reffe#::ReferenceFE
    U#::FESpace
    Ω#::Triangulation
    dΩ#::Measure
    E#::CellField
    ScalarModeFEM(Name,neff,lambda,reffe,U,Ω,dΩ,E)=new(Name,neff,lambda,reffe,U,Ω,dΩ,E)
    ScalarModeFEM(Name,neff,lambda)=new(Name,neff,lambda,[],[],[],[],[])
end

function Base.show(io::IO, ::MIME"text/plain",f::ScalarModeFEM) 
    if isempty(f.reffe) 
        print(io,"Name = ",f.Name,"\nneff = ",f.neff,"\nlambda = ",f.lambda,"\nfield: no");
    else
        print(io,"Name = ",f.Name,"\nneff = ",f.neff,"\nlambda = ",f.lambda,"\nfield: yes");
    end
end
Base.show(io::IO, f::ScalarModeFEM) = isempty(f.reffe) ? print(io,"[",f.Name,",",f.neff,",",f.lambda,",no]") : print(io,"[",f.Name,",",f.neff,",",f.lambda,",yes]");

function writevtk(name::String,m::ScalarModeFEM)
    if isempty(m.reffe)
        throw(DomainError(m, "The mode does not contain a field"));
    else
        Gridap.writevtk(m.Ω,name,cellfields=["real(E)"=>real(m.E),"imag(E)"=>imag(m.E)]);
    end
    return nothing;
end

mutable struct VectorModeFEM <: Mode
    Name::String
    neff::ComplexF64
    lambda::Float64
    reffe1#::ReferenceFE
    reffe2#::ReferenceFE
    U1#::FESpace
    U2#::FESpace
    Ω#::Triangulation
    dΩ#::Measure
    Ex#::CellField
    Ey#::CellField
    Ez#::CellField
    Hx#::CellField
    Hy#::CellField
    Hz#::CellField
    VectorModeFEM(Name,neff,lambda,reffe1,reffe2,U1,U2,Ω,dΩ,Ex,Ey,Ez,Hx,Hy,Hz)=new(Name,neff,lambda,reffe1,reffe2,U1,U2,Ω,dΩ,Ex,Ey,Ez,Hx,Hy,Hz)
    VectorModeFEM(Name,neff,lambda)=new(Name,neff,lambda,[],[],[],[],[],[],[],[],[],[],[],[])
end

function Base.show(io::IO, ::MIME"text/plain",f::VectorModeFEM) 
    if isempty(f.reffe1) 
        print(io,"Name = ",f.Name,"\nneff = ",f.neff,"\nlambda = ",f.lambda,"\nfield: no");
    else
        print(io,"Name = ",f.Name,"\nneff = ",f.neff,"\nlambda = ",f.lambda,"\nfield: yes");
    end
end
Base.show(io::IO, f::VectorModeFEM) = isempty(f.reffe1) ? print(io,"[",f.Name,",",f.neff,",",f.lambda,",no]") : print(io,"[",f.Name,",",f.neff,",",f.lambda,",yes]");

function writevtk(name::String,m::VectorModeFEM)
    if isempty(m.reffe1)
        throw(DomainError(m, "The mode does not contain a field"));
    else
        Gridap.writevtk(m.Ω,name,cellfields=["real(Ex)"=>real(m.Ex),"imag(Ex)"=>imag(m.Ex),"real(Ey)"=>real(m.Ey),"imag(Ey)"=>imag(m.Ey),"real(Ez)"=>real(m.Ez),"imag(Ez)"=>imag(m.Ez),"real(Hx)"=>real(m.Hx),"imag(Hx)"=>imag(m.Hx),"real(Hy)"=>real(m.Hy),"imag(Hy)"=>imag(m.Hy),"real(Hz)"=>real(m.Hz),"imag(Hz)"=>imag(m.Hz)]);
    end
    return nothing;
end
############################# Mode conversion #############################

"""
    ScalarMode2D(m::ScalarMode1D;sincos::Char='c')

Convert a ScalarMode1D into a ScalarMode2D

- sincos must be 'c' for a field in cos(nu.θ) or 's' for a field in sin(nu.θ) if m.nu≠0
"""
function ScalarMode2D(m::ScalarMode1D;sincos::Char='c')
    if (!in(sincos,['s','c']))
        throw(DomainError(sincos, "sincos must be 's' or 'c'"))
    end
    if ((m.nu==0) && (sincos=='s'))
        throw(DomainError(sincos, "sincos must be 'c' if nu=0"))
    end
    if isempty(m.r)
        return ScalarMode2D(m.Name,m.neff,m.lambda);
    else
        xmax=maximum(m.r);
        nbPoints=2*length(m.r)-1;
        x=collect(LinRange(-xmax,xmax,nbPoints));
        R=hypot.(x,x')
        theta=atan.(x',x);
        interp=LinearInterpolation(m.r,m.E,extrapolation_bc=0);
        E=interp.(R);
        if sincos=='c'
            @. E=E*cos(theta*m.nu);
        else
            @. E=E*sin(theta*m.nu);
        end
        return ScalarMode2D(m.Name,m.neff,m.lambda,x,x,E);
    end
end

"""
    VectorMode(m::ScalarMode2D;polar::Char='x')

Convert a ScalarMode2D into a x- or y-polarized VectorMode

- polar must be 'x' or 'y'
"""
function VectorMode(m::ScalarMode2D;polar::Char='x')
    if (!in(polar,['x','y']))
        throw(DomainError(polar, "polar must be 'x' or 'y'"))
    end
    if isempty(m.x)
        return VectorMode(m.Name,m.neff,m.lambda);
    else
        z=zeros(size(m.E));
        if polar=='x'
            return VectorMode(m.Name,m.neff,m.lambda,m.x,m.y,m.E,z,z,z,m.E*m.neff/c/mu0,z);
        else
            return VectorMode(m.Name,m.neff,m.lambda,m.x,m.y,z,m.E,z,-m.E*m.neff/c/mu0,z,z);
        end
    end
end

"""
    VectorMode(m::ScalarMode1D;polar::Char='x',sincos::Char='c')

Convert a ScalarMode1D into a x- or y-polarized VectorMode

- polar must be 'x' or 'y'
- sincos must be 'c' for a field in cos(nu.θ) or 's' for a field in sin(nu.θ) if m.nu≠0
"""
function VectorMode(m::ScalarMode1D;polar::Char='x',sincos::Char='c')
    VectorMode(ScalarMode2D(m,sincos=sincos),polar=polar)
end

############################# Poynting Vector #############################

"""
    PoyntingVector(f::VectorField)

Return a tuple of 3 matrix that describes the Poynting Vector (Px,Py,Pz)
"""
function PoyntingVector(f::VectorField)
    0.5*real(f.Ey.*conj(f.Hz)-f.Ez.*conj(f.Hy)),0.5*real(f.Ez.*conj(f.Hx)-f.Ex.*conj(f.Hz)),0.5*real(f.Ex.*conj(f.Hy)-f.Ey.*conj(f.Hx))
end

"""
    PoyntingVector(m::VectorMode)

Return a tuple of 3 matrix that describes the Poynting Vector (Px,Py,Pz)
"""
function PoyntingVector(m::VectorMode)
    zeros(size(m.Ex)),zeros(size(m.Ex)),0.5*(m.Ex.*m.Hy-m.Ey.*m.Hx)
end

"""
    PoyntingVector(m::ScalarMode1D;sincos::Char='c')

Return a tuple of 3 matrix that describes the Poynting Vector (Px,Py,Pz)

- sincos must be 'c' for a field in cos(nu.θ) or 's' for a field in sin(nu.θ) if m.nu≠0
"""
function PoyntingVector(m::ScalarMode1D;sincos::Char='c')
    PoyntingVector(VectorMode(m,sincos=sincos))
end

"""
    PoyntingVector(m::ScalarMode2D)

Return a tuple of 3 matrix that describes the Poynting Vector (Px,Py,Pz)
"""
function PoyntingVector(m::ScalarMode2D)
    PoyntingVector(VectorMode(m))
end

############################# Conversion from mode to field i.e. propagation #############################

"""
    ScalarField(m::ScalarMode2D;z::Real=0)

Returns the scalar field due to the mode m after a propagation distance z 

- z: Propagation distance - Must be in the same unit as m.lambda
"""
function ScalarField(m::ScalarMode2D;z::Real=0)
    ScalarField(m.x,m.y,m.E*exp(im*z*2*pi/m.lambda*m.neff));
end

"""
    ScalarField(m::ScalarMode1D;sincos::Char='c',z::Real=0)

Returns the scalar field due to the mode m after a propagation distance z 

- sincos must be 'c' for a field in cos(nu.θ) or 's' for a field in sin(nu.θ) if m.nu≠0
- z: Propagation distance - Must be in the same unit as m.lambda
"""
function ScalarField(m::ScalarMode1D;sincos::Char='c',z::Real=0)
    ScalarField(ScalarMode2D(m,sincos=sincos),z=z)
end

"""
    VectorField(m::VectorMode;z::Real=0)

Returns the vector field due to the mode m after a propagation distance z 

- z: Propagation distance - Must be in the same unit as m.lambda
"""
function VectorField(m::VectorMode;z::Real=0)
    p=exp(im*z*2*pi/m.lambda*m.neff);
    VectorField(m.x,m.y,m.Ex*p,m.Ey*p,m.Ez*im*p,m.Hx*p,m.Hy*p,m.Hz*im*p);
end

"""
    VectorField(m::ScalarMode1D;polar::Char='x',sincos::Char='c',z::Real=0)

Returns the vector field due to the mode m after a propagation distance z 

- polar must be 'x' or 'y'
- sincos must be 'c' for a field in cos(nu.θ) or 's' for a field in sin(nu.θ) if m.nu≠0
- z: Propagation distance - Must be in the same unit as m.lambda
"""
function VectorField(m::ScalarMode1D;polar::Char='x',sincos::Char='c',z::Real=0)
    VectorField(VectorMode(m,polar=polar,sincos=sincos),z=z)
end

"""
    VectorField(m::ScalarMode2D;polar::Char='x',z::Real=0)

Returns the vector field due to the mode m after a propagation distance z 

- polar must be 'x' or 'y'
- z: Propagation distance - Must be in the same unit as m.lambda
"""
function VectorField(m::ScalarMode2D;polar::Char='x',z::Real=0)
    VectorField(VectorMode(m,polar=polar),z=z)
end


############################# Field interpolation #############################
"""
    Interpolation(f::ScalarField,x::Vector{<:Real},y::Vector{<:Real})

Returns a scalar field obtained by the linear interpolation of the field f.  

This function assumes that the electric field is null outside the box defined by f.x and f.y
"""
function Interpolation(f::ScalarField,x::Vector{<:Real},y::Vector{<:Real})
    if isempty(f.x)
        return ScalarField([],[],Matrix{Float64}(undef,0,0));
    else
        interp=LinearInterpolation((f.y,f.x),f.E,extrapolation_bc=0);
        E=interp.(x,y');
        return ScalarField(x,y,E);
    end
end

"""
    Interpolation(f::VectorField,x::Vector{<:Real},y::Vector{<:Real})

Returns a vector field obtained by the linear interpolation of the field f.  

This function assumes that the electric/magnetic field is null outside the box defined by f.x and f.y
"""
function Interpolation(f::VectorField,x::Vector{<:Real},y::Vector{<:Real})
    if isempty(f.x)
        return VectorField([],[],Matrix{Float64}(undef,0,0),Matrix{Float64}(undef,0,0),Matrix{Float64}(undef,0,0),Matrix{Float64}(undef,0,0),Matrix{Float64}(undef,0,0),Matrix{Float64}(undef,0,0));
    else
        interp=LinearInterpolation((f.x,f.y),f.Ex,extrapolation_bc=0);
        Ex=interp.(x,y');
        interp=LinearInterpolation((f.x,f.y),f.Ey,extrapolation_bc=0);
        Ey=interp.(x,y');
        interp=LinearInterpolation((f.x,f.y),f.Ez,extrapolation_bc=0);
        Ez=interp.(x,y');
        interp=LinearInterpolation((f.x,f.y),f.Hx,extrapolation_bc=0);
        Hx=interp.(x,y');
        interp=LinearInterpolation((f.x,f.y),f.Hy,extrapolation_bc=0);
        Hy=interp.(x,y');
        interp=LinearInterpolation((f.x,f.y),f.Hz,extrapolation_bc=0);
        Hz=interp.(x,y');
        return VectorField(x,y,Ex,Ey,Ez,Hx,Hy,Hz);
    end
end

function Interpolation(f::Field,x::AbstractRange,y::Vector{<:Real})
    Interpolation(f,collect(x),y)
end

function Interpolation(f::Field,x::Vector{<:Real},y::AbstractRange)
    Interpolation(f,x,collect(y))
end

function Interpolation(f::Field,x::AbstractRange,y::AbstractRange)
    Interpolation(f,collect(x),collect(y))
end

############################# Normalization #############################
"""
    normalize!(m::ScalarMode1D;unitIntegral::Bool=true)

Normalize the mode m with the method given by the boolean parameter `unitIntegral`.
"""
function normalize!(m::ScalarMode1D;unitIntegral::Bool=true)
    integral=trapz(m.r,2*pi*m.r.*((m.E).^2));
    if (m.nu!=0)
        integral=integral*0.5;
    end
    if unitIntegral
        (m.E).=(m.E)/sqrt(integral);
    else
        (m.E).=(m.E)*sqrt(2*mu0*c/m.neff/integral);
    end
    return nothing
end

"""
    normalize!(m::ScalarMode2D;unitIntegral::Bool=true)

Normalize the mode m with the method given by the boolean parameter `unitIntegral`.
"""
function normalize!(m::ScalarMode2D;unitIntegral::Bool=true)
    integral=trapz((m.x,m.y),((m.E).^2));
    if unitIntegral
        (m.E).=(m.E)/sqrt(integral);
    else
        (m.E).=(m.E)*sqrt(2*mu0*c/m.neff/integral);
    end
    return nothing
end

"""
    normalize!(m::VectorMode)

Normalize the mode m.
"""
function normalize!(m::VectorMode)
    integral=abs(trapz((m.x,m.y),m.Ex.*m.Hy-m.Ey.*m.Hx)*0.5);
    (m.Ex).=(m.Ex)/sqrt(integral);
    (m.Ey).=(m.Ey)/sqrt(integral);
    (m.Ez).=(m.Ez)/sqrt(integral);
    (m.Hx).=(m.Hx)/sqrt(integral);
    (m.Hy).=(m.Hy)/sqrt(integral);
    (m.Hz).=(m.Hz)/sqrt(integral);
    return nothing
end

function normalize!(m::ScalarModeFEM;unitIntegral::Bool=true)
    if !isempty(m.reffe)
        integral=sum(integrate(abs2(m.E),m.dΩ));
        if unitIntegral
            m.E=m.E/sqrt(integral);
        else
            m.E=m.E*sqrt(2*mu0*c/m.neff/integral);
        end
    end
end

function normalize!(m::VectorModeFEM)
    if !isempty(m.reffe1)
        integral=0.5*abs(sum(integrate(m.Ex*conj(m.Hy)-m.Ey*conj(m.Hx),m.dΩ)));
        m.Ex=m.Ex/sqrt(integral);
        m.Ey=m.Ey/sqrt(integral);
        m.Ez=m.Ez/sqrt(integral);
        m.Hx=m.Hx/sqrt(integral);
        m.Hy=m.Hy/sqrt(integral);
        m.Hz=m.Hz/sqrt(integral);
    end
end

############################# Overlap #############################
"""
    overlap(f1::ScalarField,f2::ScalarField)
"""
function overlap(f1::ScalarField,f2::ScalarField)
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
    overlap(f1::VectorField,f2::VectorField)
"""
function overlap(f1::VectorField,f2::VectorField)
    if ((f1.x==f2.x) && (f1.y==f2.y))
        return trapz((f2.x,f2.y),f1.Ex.*conj.(f2.Hy)-f1.Ey.*conj.(f2.Hx))*0.5;
    else
        x=unique(sort(vcat(f1.x,f2.x)));
        y=unique(sort(vcat(f1.y,f2.y)));
        interp=LinearInterpolation((f1.x,f1.y),f1.Ex,extrapolation_bc=0);
        Ex1=interp.(x,y');
        interp=LinearInterpolation((f1.x,f1.y),f1.Ey,extrapolation_bc=0);
        Ey1=interp.(x,y');
        interp=LinearInterpolation((f1.x,f1.y),f1.Hx,extrapolation_bc=0);
        Hx1=interp.(x,y');
        interp=LinearInterpolation((f1.x,f1.y),f1.Hy,extrapolation_bc=0);
        Hy1=interp.(x,y');
        interp=LinearInterpolation((f2.x,f2.y),f2.Ex,extrapolation_bc=0);
        Ex2=interp.(x,y');
        interp=LinearInterpolation((f2.x,f2.y),f2.Ey,extrapolation_bc=0);
        Ey2=interp.(x,y');
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
function overlap(m1::ScalarMode1D,m2::ScalarMode1D)
    if ((isempty(m1.r)) || (isempty(m2.r)))
        throw(ArgumentError("At least one of the two modes does not have a field"));
    end
    if m1.nu!=m2.nu
        return 0.0;
    end
    if (m1.r==m2.r)
        integral1=trapz(m1.r,m1.r.*(m1.E.^2));
        integral2=trapz(m2.r,m2.r.*(m2.E.^2));
        integral=trapz(m1.r,m1.r.*m1.E.*m2.E);
        return integral/sqrt(integral1*integral2);
    else
        r=unique(sort(vcat(m1.r,m2.r)));
        interp1=LinearInterpolation(m1.r,m1.E,extrapolation_bc=0);
        E1=interp1.(r);
        interp2=LinearInterpolation(m2.r,m2.E,extrapolation_bc=0);
        E2=interp2.(r);
        integral1=trapz(r,r.*(E1.^2));
        integral2=trapz(r,r.*(E2.^2));
        integral=trapz(r,r.*E1.*E2);
        return integral/sqrt(integral1*integral2);
    end
end

"""
    overlap(m1::ScalarMode2D,m2::ScalarMode2D)
"""
function overlap(m1::ScalarMode2D,m2::ScalarMode2D)
    if ((isempty(m1.x)) || (isempty(m2.x)))
        throw(ArgumentError("At least one of the two modes does not have a field"));
    end
    if ((m1.x==m2.x) && (m1.y==m2.y))
        integral1=trapz((m1.x,m1.y),((m1.E).^2));
        integral2=trapz((m2.x,m2.y),((m2.E).^2));
        integral=trapz((m1.x,m1.y),(m2.E.*m1.E));
        return integral/sqrt(integral1*integral2);
    else
        x=unique(sort(vcat(m1.x,m2.x)));
        y=unique(sort(vcat(m1.y,m2.y)));
        interp=LinearInterpolation((m1.x,m1.y),m1.E,extrapolation_bc=0);
        E1=interp.(x,y');
        interp=LinearInterpolation((m2.x,m2.y),m2.E,extrapolation_bc=0);
        E2=interp.(x,y');
        integral1=trapz((x,y),E1.^2);
        integral2=trapz((x,y),E2.^2);
        integral=trapz((x,y),E1.*E2);
        return integral/sqrt(integral1*integral2);
    end
end

"""
    overlap(f1::ScalarField,m2::ScalarMode2D)
"""
function overlap(f1::ScalarField,m2::ScalarMode2D)
    if ((isempty(f1.x)) || (isempty(m2.x)))
        throw(ArgumentError("At least one of the two arguments does not have a field"));
    end
    if ((f1.x==m2.x) && (f1.y==m2.y))
        integral2=trapz((m2.x,m2.y),((m2.E).^2));
        integral=trapz((f1.x,f1.y),(m2.E.*f1.E));
        return integral/sqrt(integral2);
    else
        x=unique(sort(vcat(f1.x,m2.x)));
        y=unique(sort(vcat(f1.y,m2.y)));
        interp=LinearInterpolation((f1.x,f1.y),f1.E,extrapolation_bc=0);
        E1=interp.(x,y');
        interp=LinearInterpolation((m2.x,m2.y),m2.E,extrapolation_bc=0);
        E2=interp.(x,y');
        integral2=trapz((x,y),E2.^2);
        integral=trapz((x,y),E1.*E2);
        return integral/sqrt(integral2);
    end
end

"""
    overlap(m1::ScalarMode2D,f2::ScalarField)
"""
function overlap(m1::ScalarMode2D,f2::ScalarField)
    return conj(overlap(f2,m1));
end

"""
    overlap(m1::VectorMode,m2::VectorMode)
"""
function overlap(m1::VectorMode,m2::VectorMode)
    if ((isempty(m1.x)) || (isempty(m2.x)))
        throw(ArgumentError("At least one of the two modes does not have a field"));
    end
    if ((m1.x==m2.x) && (m1.y==m2.y))
        integral1=trapz((m1.x,m1.y),m1.Ex.*m1.Hy-m1.Ey.*m1.Hx)*0.5;
        integral2=trapz((m2.x,m2.y),m2.Ex.*m2.Hy-m2.Ey.*m2.Hx)*0.5;
        integral=trapz((m2.x,m2.y),m1.Ex.*m2.Hy-m1.Ey.*m2.Hx)*0.5;
        return integral/sqrt(integral1*integral2);
    else
        x=unique(sort(vcat(m1.x,m2.x)));
        y=unique(sort(vcat(m1.y,m2.y)));
        interp=LinearInterpolation((m1.x,m1.y),m1.Ex,extrapolation_bc=0);
        Ex1=interp.(x,y');
        interp=LinearInterpolation((m1.x,m1.y),m1.Ey,extrapolation_bc=0);
        Ey1=interp.(x,y');
        interp=LinearInterpolation((m1.x,m1.y),m1.Hx,extrapolation_bc=0);
        Hx1=interp.(x,y');
        interp=LinearInterpolation((m1.x,m1.y),m1.Hy,extrapolation_bc=0);
        Hy1=interp.(x,y');
        interp=LinearInterpolation((m2.x,m2.y),m2.Ex,extrapolation_bc=0);
        Ex2=interp.(x,y');
        interp=LinearInterpolation((m2.x,m2.y),m2.Ey,extrapolation_bc=0);
        Ey2=interp.(x,y');
        interp=LinearInterpolation((m2.x,m2.y),m2.Hx,extrapolation_bc=0);
        Hx2=interp.(x,y');
        interp=LinearInterpolation((m2.x,m2.y),m2.Hy,extrapolation_bc=0);
        Hy2=interp.(x,y');
        integral1=trapz((x,y),Ex1.*Hy1-Ey1.*Hx1)*0.5;
        integral2=trapz((x,y),Ex2.*Hy2-Ey2.*Hx2)*0.5;
        integral=trapz((x,y),Ex1.*Hy2-Ey1.*Hx2)*0.5;
        return integral/sqrt(integral1*integral2);
    end
end

"""
    overlap(f1::VectorField,m2::VectorMode)
"""
function overlap(f1::VectorField,m2::VectorMode)
    if ((isempty(f1.x)) || (isempty(m2.x)))
        throw(ArgumentError("At least one of the two modes does not have a field"));
    end
    if ((f1.x==m2.x) && (f1.y==m2.y))
        integral2=trapz((m2.x,m2.y),m2.Ex.*m2.Hy-m2.Ey.*m2.Hx)*0.5;
        integral=trapz((m2.x,m2.y),f1.Ex.*m2.Hy-f1.Ey.*m2.Hx)*0.5;
        return integral/sqrt(integral2);
    else
        x=unique(sort(vcat(f1.x,m2.x)));
        y=unique(sort(vcat(f1.y,m2.y)));
        interp=LinearInterpolation((f1.x,f1.y),f1.Ex,extrapolation_bc=0);
        Ex1=interp.(x,y');
        interp=LinearInterpolation((f1.x,f1.y),f1.Ey,extrapolation_bc=0);
        Ey1=interp.(x,y');
        interp=LinearInterpolation((f1.x,f1.y),f1.Hx,extrapolation_bc=0);
        Hx1=interp.(x,y');
        interp=LinearInterpolation((f1.x,f1.y),f1.Hy,extrapolation_bc=0);
        Hy1=interp.(x,y');
        interp=LinearInterpolation((m2.x,m2.y),m2.Ex,extrapolation_bc=0);
        Ex2=interp.(x,y');
        interp=LinearInterpolation((m2.x,m2.y),m2.Ey,extrapolation_bc=0);
        Ey2=interp.(x,y');
        interp=LinearInterpolation((m2.x,m2.y),m2.Hx,extrapolation_bc=0);
        Hx2=interp.(x,y');
        interp=LinearInterpolation((m2.x,m2.y),m2.Hy,extrapolation_bc=0);
        Hy2=interp.(x,y');
        integral2=trapz((x,y),Ex2.*Hy2-Ey2.*Hx2)*0.5;
        integral=trapz((x,y),Ex1.*Hy2-Ey1.*Hx2)*0.5;
        return integral/sqrt(integral2);
    end
end

"""
    overlap(m1::VectorMode,f2::VectorField)
"""
function overlap(m1::VectorMode,f2::VectorField)
    return conj(overlap(f2,m1));
end

############################# Effective area #############################
"""
    Aeff(m::ScalarMode1D)
"""
function Aeff(m::ScalarMode1D)
    E2=trapz(m.r,2*pi*m.r.*(m.E).^2);
    E4=trapz(m.r,2*pi*m.r.*(m.E).^4);
    if (m.nu!=0)
        E2=E2*0.5;
        E4=E4*3.0/8.0;
    end
    return E2^2/E4;
end

"""
    Aeff(m::ScalarMode2D)
"""
function Aeff(m::ScalarMode2D)
    E2=trapz((m.x,m.y),(m.E).^2);
    E4=trapz((m.x,m.y),(m.E).^4);
    return E2^2/E4;
end

"""
    Aeff(m::VectorMode)

Compute effective areas by assuming that n0≃neff (true in a weakly-guiding fiber)
"""
function Aeff(m::VectorMode)
    #Approximation (neff=n0)
    E2=trapz((m.x,m.y),m.Ex.*m.Hy-m.Ey.*m.Hx);
    E4_1=trapz((m.x,m.y),abs2.(abs2.(m.Ex)+abs2.(m.Ey)+abs2.(m.Ez)));
    E4_2=trapz((m.x,m.y),abs2.(abs2.(m.Ex)+abs2.(m.Ey)-abs2.(m.Ez)));
    return E2^2/E4_1*mu0/eps0/(m.neff)^2,E2^2/E4_2*mu0/eps0/(m.neff)^2;
end

"""
    Aeff(m::VectorMode)

Compute effective areas by assuming that n0 is uniform
"""
function Aeff(m::VectorMode,n0::Real)
    E2=trapz((m.x,m.y),m.Ex.*m.Hy-m.Ey.*m.Hx);
    E4_1=trapz((m.x,m.y),abs2.(abs2.(m.Ex)+abs2.(m.Ey)+abs2.(m.Ez)));
    E4_2=trapz((m.x,m.y),abs2.(abs2.(m.Ex)+abs2.(m.Ey)-abs2.(m.Ez)));
    return E2^2/E4_1*mu0/eps0/n0^2,E2^2/E4_2*mu0/eps0/n0^2;
end

"""
    Aeff(m::VectorMode)

n0 must be a function of the cartesian coordinates x and y
"""
function Aeff(m::VectorMode,n0::Function)
    if (first(methods(n0)).nargs!=3)
        throw(DomainError(n0, "The refractive index function must have 2 argument"));
    end
    E2=trapz((m.x,m.y),m.Ex.*m.Hy-m.Ey.*m.Hx);
    E4_1=trapz((m.x,m.y),(abs2.(abs2.(m.Ex)+abs2.(m.Ey)+abs2.(m.Ez))).*(n0.(m.x',m.y).^2));
    E4_2=trapz((m.x,m.y),(abs2.(abs2.(m.Ex)+abs2.(m.Ey)-abs2.(m.Ez))).*(n0.(m.x',m.y).^2));
    return E2^2/E4_1*mu0/eps0,E2^2/E4_2*mu0/eps0;
end

function Aeff(m::ScalarModeFEM)
    E2=sum(integrate(abs2(m.E),m.dΩ));
    E4=sum(integrate(abs2(abs2(m.E)),m.dΩ));
    return E2^2/E4;
end

function Aeff(m::VectorModeFEM)
    #Approximation (real(neff)=n0)
    E2=sum(integrate(m.Ex*conj(m.Hy)-m.Ey*conj(m.Hx),m.dΩ));
    E4_1=sum(integrate(abs2(abs2(m.Ex)+abs2(m.Ey)+abs2(m.Ez)),m.dΩ));
    E4_2=sum(integrate(abs2(m.Ex*m.Ex+m.Ey*m.Ey+m.Ez*m.Ez),m.dΩ));
    return real(E2)^2/E4_1*mu0/eps0/(real(m.neff))^2,real(E2)^2/E4_2*mu0/eps0/(real(m.neff))^2;
end

############################# Nonlinear Coefficient #############################

"""
    nonLinearCoefficient(m::Mode,n2::Real)
"""
function nonLinearCoefficient(m::Mode,n2::Real)
    return n2*2*pi/m.lambda./Aeff(m);
end

"""
    nonLinearCoefficient(m::VectorMode,n0::Real,n2::Real)
"""
function nonLinearCoefficient(m::VectorMode,n0::Real,n2::Real)
    return n2*2*pi/m.lambda./Aeff(m,n0);
end

"""
    nonLinearCoefficient(m::ScalarMode1D,n2::Function)

n2 must be a function of the radial coordinate r
"""
function nonLinearCoefficient(m::ScalarMode1D,n2::Function)
    if (first(methods(n2)).nargs!=2)
        throw(DomainError(n2, "The non-linear index function must have 1 argument"));
    end
    E2=trapz(m.r,2*pi*m.r.*(m.E).^2);
    E4=trapz(m.r,2*pi*(m.r.*(m.E).^4).*n2.(m.r));
    if (m.nu!=0)
        E2=E2*0.5;
        E4=E4*3.0/8.0;
    end
    omega=2*pi*c/m.lambda;
    return E4/(E2^2)*omega/c;
end

"""
    nonLinearCoefficient(m::ScalarMode2D,n2::Function)

n2 must be a function of the cartesian coordinates x and y
"""
function nonLinearCoefficient(m::ScalarMode2D,n2::Function)
    if (first(methods(n2)).nargs!=3)
        throw(DomainError(n2, "The non-linear index function must have 2 argument"));
    end
    #X,Y=meshgrid(m.x,m.y);
    E2=trapz((m.x,m.y),(m.E).^2);
    E4=trapz((m.x,m.y),((m.E).^4).*n2.(m.x,m.y'));
    omega=2*pi*c/m.lambda;
    return E4/(E2^2)*omega/c;
end

"""
    nonLinearCoefficient(m::VectorMode,n0::Function,n2::Function)

n0 and n2 must be functions of the cartesian coordinates x and y
"""
function nonLinearCoefficient(m::VectorMode,n0::Function,n2::Function)
    if (first(methods(n0)).nargs!=3)
        throw(DomainError(n0, "The refractive index function must have 2 argument"));
    end
    if (first(methods(n2)).nargs!=3)
        throw(DomainError(n2, "The non-linear index function must have 2 argument"));
    end
    #X,Y=meshgrid(m.x,m.y);
    E2=trapz((m.x,m.y),m.Ex.*m.Hy-m.Ey.*m.Hx);
    E4_1=trapz((m.x,m.y),(abs2.(abs2.(m.Ex)+abs2.(m.Ey)+abs2.(m.Ez))).*(n0.(m.x,m.y').^2).*n2.(m.x,m.y'));
    E4_2=trapz((m.x,m.y),(abs2.(abs2.(m.Ex)+abs2.(m.Ey)-abs2.(m.Ez))).*(n0.(m.x,m.y').^2).*n2.(m.x,m.y'));
    omega=2*pi*c/m.lambda;
    return 1.0/(E2^2/E4_1*mu0/eps0)*omega/c,1.0/(E2^2/E4_2*mu0/eps0)*omega/c;
end

############################# MFD #############################
"""
    MFD(m::ScalarMode1D)
"""
function MFD(m::ScalarMode1D)
    pos=argmin(abs.(abs.(m.E).-(maximum(abs.(m.E))/exp(1))));
    return 2*m.r[pos];
end

"""
    MFD(m::ScalarMode2D;angle::Real=0)

Compute the MFD in the direction given by the angle (angle with the x-axis)
"""
function MFD(m::ScalarMode2D;angle::Real=0)
    rmax=sqrt(maximum(m.x)^2+maximum(m.y)^2);
    dr=minimum([minimum(diff(m.x)),minimum(diff(m.y))]);
    nb=Integer(2*ceil(rmax/dr)+1);
    r=collect(LinRange(-rmax,rmax,nb));
    x=r*cosd(angle);
    y=r*sind(angle);
    interp=LinearInterpolation((m.x,m.y),m.E,extrapolation_bc=0);
    E=interp.(x,y);
    posmax=argmax(abs.(E));
    pos=argmin(abs.(abs.(E).-(maximum(abs.(E))/exp(1))));
    return abs(2*(r[pos]-r[posmax]));
end

"""
    MFD(m::VectorMode;angle::Real=0)

Compute the MFD in the direction given by the angle (angle with the x-axis)
"""
function MFD(m::VectorMode;angle::Real=0)
    rmax=sqrt(maximum(m.x)^2+maximum(m.y)^2);
    dr=minimum([minimum(diff(m.x)),minimum(diff(m.y))]);
    nb=Integer(2*ceil(rmax/dr)+1);
    r=collect(LinRange(-rmax,rmax,nb));
    x=r*cosd(angle);
    y=r*sind(angle);
    Pz=m.Ex.*m.Hy-m.Ey.*m.Hx;
    interp=LinearInterpolation((m.x,m.y),Pz,extrapolation_bc=0);
    E2=interp.(x,y);
    posmax=argmax(abs.(E2));
    pos=argmin(abs.(abs.(E2).-(maximum(abs.(E2))/exp(2))));
    return abs(2*(r[pos]-r[posmax]));
end






