############################# Fields #############################

"""
Abstract structure to describe an electromagnetic field

`abstract type Field end`
"""
abstract type Field end
Broadcast.:broadcastable(f::Field)=Ref(f)

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

"""
Structure describing a scalar field defined with a Triangulation  

- Ω :: `Gridap.Triangulation`
- dΩ :: `Gridap.CellData.Measure`
- E :: `Gridap.CellField` - Electric field
"""
mutable struct ScalarFieldFEM <: Field
    Ω::Gridap.Triangulation
    dΩ::Gridap.CellData.Measure
    E::Gridap.CellField
end

function Base.:+(f1::ScalarFieldFEM, f2::ScalarFieldFEM)
    if (f1.Ω==f2.Ω)
        ScalarFieldFEM(f1.Ω,f1.dΩ,f1.E+f2.E)
    else
        throw(ArgumentError("Both fields must have the same triangulation"))
    end
end

function Base.:-(f1::ScalarFieldFEM, f2::ScalarFieldFEM)
    if (f1.Ω==f2.Ω)
        ScalarFieldFEM(f1.Ω,f1.dΩ,f1.E-f2.E)
    else
        throw(ArgumentError("Both fields must have the same triangulation"))
    end
end

Base.:*(k::Number, f::ScalarFieldFEM) = ScalarFieldFEM(f.Ω,f.dΩ,k*f.E);
Base.:*(f::ScalarFieldFEM,k::Number) = ScalarFieldFEM(f.Ω,f.dΩ,k*f.E);
Base.:/(f::ScalarFieldFEM,k::Number) = ScalarFieldFEM(f.Ω,f.dΩ,f.E/k);
Base.:+(f::ScalarFieldFEM) = ScalarFieldFEM(f.Ω,f.dΩ,f.E);
Base.:-(f::ScalarFieldFEM) = ScalarFieldFEM(f.Ω,f.dΩ,-f.E);
Base.:real(f::ScalarFieldFEM) = ScalarFieldFEM(f.Ω,f.dΩ,real(f.E));
Base.:imag(f::ScalarFieldFEM) = ScalarFieldFEM(f.Ω,f.dΩ,imag(f.E));
Base.:conj(f::ScalarFieldFEM) = ScalarFieldFEM(f.Ω,f.dΩ,conj(f.E));

"""
Structure describing a vector field defined with a Triangulation   

- Ω :: `Gridap.Triangulation`
- dΩ :: `Gridap.CellData.Measure`
- Ex :: `Gridap.CellField` - x-component of the electric field
- Ey :: `Gridap.CellField` - y-component of the electric field
- Ez :: `Gridap.CellField` - z-component of the electric field
- Hx :: `Gridap.CellField` - x-component of the magnetic field
- Hy :: `Gridap.CellField` - y-component of the magnetic field
- Hz :: `Gridap.CellField` - z-component of the magnetic field
"""
mutable struct VectorFieldFEM <: Field
    Ω::Gridap.Triangulation
    dΩ::Gridap.CellData.Measure
    Ex::Gridap.CellField
    Ey::Gridap.CellField
    Ez::Gridap.CellField
    Hx::Gridap.CellField
    Hy::Gridap.CellField
    Hz::Gridap.CellField
end

function Base.:+(f1::VectorFieldFEM, f2::VectorFieldFEM)
    if (f1.Ω==f2.Ω)
        VectorFieldFEM(f1.Ω,f1.dΩ,f1.Ex+f2.Ex,f1.Ey+f2.Ey,f1.Ez+f2.Ez,f1.Hx+f2.Hx,f1.Hy+f2.Hy,f1.Hz+f2.Hz)
    else
        throw(ArgumentError("Both fields must have the same triangulation"))
    end
end

function Base.:-(f1::VectorFieldFEM, f2::VectorFieldFEM)
    if (f1.Ω==f2.Ω)
        VectorFieldFEM(f1.Ω,f1.dΩ,f1.Ex-f2.Ex,f1.Ey-f2.Ey,f1.Ez-f2.Ez,f1.Hx-f2.Hx,f1.Hy-f2.Hy,f1.Hz-f2.Hz)
    else
        throw(ArgumentError("Both fields must have the same triangulation"))
    end
end

Base.:*(k::Number, f::VectorFieldFEM) = VectorFieldFEM(f.Ω,f.dΩ,k*f.Ex,k*f.Ey,k*f.Ez,k*f.Hx,k*f.Hy,k*f.Hz);
Base.:*(f::VectorFieldFEM,k::Number) = VectorFieldFEM(f.Ω,f.dΩ,k*f.Ex,k*f.Ey,k*f.Ez,k*f.Hx,k*f.Hy,k*f.Hz);
Base.:/(f::VectorFieldFEM,k::Number) = VectorFieldFEM(f.Ω,f.dΩ,f.Ex/k,f.Ey/k,f.Ez/k,f.Hx/k,f.Hy/k,f.Hz/k);
Base.:+(f::VectorFieldFEM) = VectorFieldFEM(f.Ω,f.dΩ,f.Ex,f.Ey,f.Ez,f.Hx,f.Hy,f.Hz);
Base.:-(f::VectorFieldFEM) = VectorFieldFEM(f.Ω,f.dΩ,-f.Ex,-f.Ey,-f.Ez,-f.Hx,-f.Hy,-f.Hz);
Base.:real(f::VectorFieldFEM) = VectorFieldFEM(f.Ω,f.dΩ,real(f.Ex),real(f.Ey),real(f.Ez),real(f.Hx),real(f.Hy),real(f.Hz));
Base.:imag(f::VectorFieldFEM) = VectorFieldFEM(f.Ω,f.dΩ,imag(f.Ex),imag(f.Ey),imag(f.Ez),imag(f.Hx),imag(f.Hy),imag(f.Hz));
Base.:conj(f::VectorFieldFEM) = VectorFieldFEM(f.Ω,f.dΩ,conj(f.Ex),conj(f.Ey),conj(f.Ez),conj(f.Hx),conj(f.Hy),conj(f.Hz));

############################# Modes #############################

"""
Abstract structure to describe an optical fiber mode

`abstract type Mode end`
"""
abstract type Mode end
Broadcast.:broadcastable(f::Mode)=Ref(f)

"""
Structure describing a scalar mode of an optical fiber with a revolution symmetry

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
Structure describing a scalar mode of an optical fiber in cartesian coordinates

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
Structure describing a vector mode of an optical fiber in cartesian coordinates.  

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

"""
Structure describing a scalar mode of an optical fiber computed with FEM.  

- Name :: `String` - Name of the mode
- neff :: `Float64` - Effective index
- lambda :: `Float64` - the wavelength at which the mode was calculated
- Ω :: `Gridap.Triangulation`
- dΩ :: `Gridap.CellData.Measure`
- E :: `Gridap.CellField` - Electric field
"""
mutable struct ScalarModeFEM <: Mode
    Name::String
    neff::ComplexF64
    lambda::Float64
    Ω#::Gridap.Triangulation
    dΩ#::Gridap.CellData.Measure
    E#::Gridap.CellField
    ScalarModeFEM(Name,neff,lambda,Ω,dΩ,E)=new(Name,neff,lambda,Ω,dΩ,E)
    ScalarModeFEM(Name,neff,lambda)=new(Name,neff,lambda,[],[],[])
end

function Base.show(io::IO, ::MIME"text/plain",f::ScalarModeFEM) 
    if (f.E==[]) 
        print(io,"Name = ",f.Name,"\nneff = ",f.neff,"\nlambda = ",f.lambda,"\nfield: no");
    else
        print(io,"Name = ",f.Name,"\nneff = ",f.neff,"\nlambda = ",f.lambda,"\nfield: yes");
    end
end
Base.show(io::IO, f::ScalarModeFEM) = (f.E==[]) ? print(io,"[",f.Name,",",f.neff,",",f.lambda,",no]") : print(io,"[",f.Name,",",f.neff,",",f.lambda,",yes]");

function writevtk(name::String,m::ScalarModeFEM)
    if (m.E==[])
        throw(DomainError(m, "The mode does not contain any field"));
    else
        Gridap.writevtk(m.Ω,name,cellfields=["real(E)"=>real(m.E),"imag(E)"=>imag(m.E)]);
    end
    return nothing;
end

"""
Structure describing a vector mode of an optical fiber computed with FEM.  

- Name :: `String` - Name of the mode
- neff :: `Float64` - Effective index
- lambda :: `Float64` - the wavelength at which the mode was calculated
- Ω :: `Gridap.Triangulation`
- dΩ :: `Gridap.CellData.Measure`
- Ex :: `Gridap.CellField` - x-component of the electric field
- Ey :: `Gridap.CellField` - y-component of the electric field
- Ez :: `Gridap.CellField` - z-component of the electric field
- Hx :: `Gridap.CellField` - x-component of the magnetic field
- Hy :: `Gridap.CellField` - y-component of the magnetic field
- Hz :: `Gridap.CellField` - z-component of the magnetic field
"""
mutable struct VectorModeFEM <: Mode
    Name::String
    neff::ComplexF64
    lambda::Float64
    Ω#::Gridap.Triangulation
    dΩ#::Gridap.CellData.Measure
    Ex#::Gridap.CellField
    Ey#::Gridap.CellField
    Ez#::Gridap.CellField
    Hx#::Gridap.CellField
    Hy#::Gridap.CellField
    Hz#::Gridap.CellField
    VectorModeFEM(Name,neff,lambda,Ω,dΩ,Ex,Ey,Ez,Hx,Hy,Hz)=new(Name,neff,lambda,Ω,dΩ,Ex,Ey,Ez,Hx,Hy,Hz)
    VectorModeFEM(Name,neff,lambda)=new(Name,neff,lambda,[],[],[],[],[],[],[],[])
end

function Base.show(io::IO, ::MIME"text/plain",f::VectorModeFEM) 
    if (f.Ex==[]) 
        print(io,"Name = ",f.Name,"\nneff = ",f.neff,"\nlambda = ",f.lambda,"\nfield: no");
    else
        print(io,"Name = ",f.Name,"\nneff = ",f.neff,"\nlambda = ",f.lambda,"\nfield: yes");
    end
end
Base.show(io::IO, f::VectorModeFEM) = (f.Ex==[]) ? print(io,"[",f.Name,",",f.neff,",",f.lambda,",no]") : print(io,"[",f.Name,",",f.neff,",",f.lambda,",yes]");

function writevtk(name::String,m::VectorModeFEM)
    if (m.Ex==[])
        throw(DomainError(m, "The mode does not contain any field"));
    else
        Gridap.writevtk(m.Ω,name,cellfields=["real(Ex)"=>real(m.Ex),"imag(Ex)"=>imag(m.Ex),"real(Ey)"=>real(m.Ey),"imag(Ey)"=>imag(m.Ey),"real(Ez)"=>real(m.Ez),"imag(Ez)"=>imag(m.Ez),"real(Hx)"=>real(m.Hx),"imag(Hx)"=>imag(m.Hx),"real(Hy)"=>real(m.Hy),"imag(Hy)"=>imag(m.Hy),"real(Hz)"=>real(m.Hz),"imag(Hz)"=>imag(m.Hz)]);
    end
    return nothing;
end

############################# Mode conversion #############################
"""
    losses(m::Mode)

Return the losses of a mode in dB/km if the wavelength is in meters.  
If the wavelength is in microns, you have to multiply the result by 1E6 to obtain the losses in dB/km.  
"""
function losses(m::Mode)
    return -4*pi/m.lambda*imag(m.neff)*10/log(10)*1000
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

"""
    VectorModeFEM(m::ScalarModeFEM;polar::Char='x')

Convert a ScalarModeFEM into a x- or y-polarized VectorModeFEM

- polar must be 'x' or 'y'
"""
function VectorModeFEM(m::ScalarModeFEM;polar::Char='x')
    if (!in(polar,['x','y']))
        throw(DomainError(polar, "polar must be 'x' or 'y'"))
    end
    if (m.E==[])
        return VectorModeFEM(m.Name,m.neff,m.lambda);
    else
        z=CellField(0,m.Ω)
        if polar=='x'
            return VectorModeFEM(m.Name,m.neff,m.lambda,m.Ω,m.dΩ,m.E,z,z,z,m.E*m.neff/c/mu0,z);
        else
            return VectorModeFEM(m.Name,m.neff,m.lambda,m.Ω,m.dΩ,z,m.E,z,-m.E*m.neff/c/mu0,z,z);
        end
    end
end

############################# Poynting Vector #############################

"""
    PoyntingVector(f::VectorField)

Return a tuple of 3 matrices.
"""
function PoyntingVector(f::VectorField)
    0.5*real(f.Ey.*conj(f.Hz)-f.Ez.*conj(f.Hy)),0.5*real(f.Ez.*conj(f.Hx)-f.Ex.*conj(f.Hz)),0.5*real(f.Ex.*conj(f.Hy)-f.Ey.*conj(f.Hx))
end

"""
    PoyntingVector(m::VectorMode)

Return a tuple of 3 matrices.
"""
function PoyntingVector(m::VectorMode)
    zeros(size(m.Ex)),zeros(size(m.Ex)),0.5*(m.Ex.*m.Hy-m.Ey.*m.Hx)
end

"""
    PoyntingVector(m::ScalarMode1D;sincos::Char='c')

Return a tuple of 3 matrices.

- sincos must be 'c' for a field in cos(nu.θ) or 's' for a field in sin(nu.θ) if m.nu≠0
"""
function PoyntingVector(m::ScalarMode1D;sincos::Char='c')
    PoyntingVector(ScalarMode2D(m,sincos=sincos))
end

"""
    PoyntingVector(m::ScalarMode2D)

Return a tuple of 3 matrices.
"""
function PoyntingVector(m::ScalarMode2D)
    zeros(size(m.E)),zeros(size(m.E)),0.5*m.neff/c/mu0*abs2.(m.E)
end

"""
    PoyntingVector(f::VectorModeFEM)

Return a tuple of 3 CellFields.
"""
function PoyntingVector(m::VectorModeFEM)
    0.5*real(m.Ey.*conj(m.Hz)-m.Ez.*conj(m.Hy)),0.5*real(m.Ez.*conj(m.Hx)-m.Ex.*conj(m.Hz)),0.5*real(m.Ex.*conj(m.Hy)-m.Ey.*conj(m.Hx))
end

"""
    PoyntingVector(f::ScalarModeFEM)

Return a tuple of 3 CellFields.
"""
function PoyntingVector(m::ScalarModeFEM)
    CellField(0,m.Ω),CellField(0,m.Ω),0.5*real(m.neff/c/mu0*abs2(m.E))
end

############################# Conversion from mode to field i.e. propagation #############################

"""
    ScalarField(m::ScalarMode2D,z::Real=0)

Returns the scalar field due to the mode m after a propagation distance z 

- z: Propagation distance - Must be in the same unit as m.lambda
"""
function ScalarField(m::ScalarMode2D,z::Real=0)
    ScalarField(m.x,m.y,m.E*exp(im*z*2*pi/m.lambda*m.neff));
end

"""
    ScalarField(m::ScalarMode1D,z::Real=0;sincos::Char='c')

Returns the scalar field due to the mode m after a propagation distance z 

- sincos must be 'c' for a field in cos(nu.θ) or 's' for a field in sin(nu.θ) if m.nu≠0
- z: Propagation distance - Must be in the same unit as m.lambda
"""
function ScalarField(m::ScalarMode1D,z::Real=0;sincos::Char='c')
    ScalarField(ScalarMode2D(m,sincos=sincos),z)
end

"""
    VectorField(m::VectorMode,z::Real=0)

Returns the vector field due to the mode m after a propagation distance z 

- z: Propagation distance - Must be in the same unit as m.lambda
"""
function VectorField(m::VectorMode,z::Real=0)
    p=exp(im*z*2*pi/m.lambda*m.neff);
    VectorField(m.x,m.y,m.Ex*p,m.Ey*p,m.Ez*im*p,m.Hx*p,m.Hy*p,m.Hz*im*p);
end

"""
    VectorField(m::ScalarMode1D,z::Real=0;polar::Char='x',sincos::Char='c')

Returns the vector field due to the mode m after a propagation distance z 

- polar must be 'x' or 'y'
- sincos must be 'c' for a field in cos(nu.θ) or 's' for a field in sin(nu.θ) if m.nu≠0
- z: Propagation distance - Must be in the same unit as m.lambda
"""
function VectorField(m::ScalarMode1D,z::Real=0;polar::Char='x',sincos::Char='c')
    VectorField(VectorMode(m,polar=polar,sincos=sincos),z)
end

"""
    VectorField(m::ScalarMode2D,z::Real=0;polar::Char='x')

Returns the vector field due to the mode m after a propagation distance z 

- polar must be 'x' or 'y'
- z: Propagation distance - Must be in the same unit as m.lambda
"""
function VectorField(m::ScalarMode2D,z::Real=0;polar::Char='x')
    VectorField(VectorMode(m,polar=polar),z)
end


"""
    ScalarFieldFEM(m::ScalarModeFEM,z::Real=0)

Returns the scalar field due to the mode m after a propagation distance z 

- z: Propagation distance - Must be in the same unit as m.lambda
"""
function ScalarFieldFEM(m::ScalarModeFEM,z::Real=0)
    p=exp(im*z*2*pi/m.lambda*m.neff);
    ScalarFieldFEM(m.Ω,m.dΩ,p*m.E)
end

"""
    VectorFieldFEM(m::VectorModeFEM,z::Real=0)

Returns the vector field due to the mode m after a propagation distance z 

- z: Propagation distance - Must be in the same unit as m.lambda
"""
function VectorFieldFEM(m::VectorModeFEM,z::Real=0)
    p=exp(im*z*2*pi/m.lambda*m.neff);
    VectorFieldFEM(m.Ω,m.dΩ,m.Ex*p,m.Ey*p,m.Ez*p,m.Hx*p,m.Hy*p,m.Hz*p);
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

"""
    normalize!(m::ScalarModeFEM;unitIntegral::Bool=true)
"""
function normalize!(m::ScalarModeFEM;unitIntegral::Bool=true)
    if (m.E!=[])
        integral=sum(integrate(abs2(m.E),m.dΩ));
        if unitIntegral
            m.E=m.E/sqrt(integral);
        else
            m.E=m.E*sqrt(2*mu0*c/m.neff/integral);
        end
    end
    return
end

"""
    normalize!(m::VectorModeFEM)
"""
function normalize!(m::VectorModeFEM)
    if (m.Ex!=[])
        integral=0.5*abs(sum(integrate(m.Ex*conj(m.Hy)-m.Ey*conj(m.Hx),m.dΩ)));
        m.Ex=m.Ex/sqrt(integral);
        m.Ey=m.Ey/sqrt(integral);
        m.Ez=m.Ez/sqrt(integral);
        m.Hx=m.Hx/sqrt(integral);
        m.Hy=m.Hy/sqrt(integral);
        m.Hz=m.Hz/sqrt(integral);
    end
    return
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

"""
    overlap(m1::ScalarModeFEM,m2::ScalarModeFEM)
"""
function overlap(m1::ScalarModeFEM,m2::ScalarModeFEM)
    if (m1.Ω!=m2.Ω)
        throw(ArgumentError("The modes must have the same triangulation Ω"));
    else
        integral1=sum(integrate(abs2(m1.E),m1.dΩ));
        integral2=sum(integrate(abs2(m2.E),m2.dΩ));
        integral=sum(integrate(m1.E*conj(m2.E),m2.dΩ));
        return integral/sqrt(integral1*integral2);
    end
end

"""
    overlap(m1::ScalarModeFEM,m2::ScalarfieldFEM)
"""
function overlap(m1::ScalarModeFEM,m2::ScalarFieldFEM)
    if (m1.Ω!=m2.Ω)
        throw(ArgumentError("The modes must have the same triangulation Ω"));
    else
        integral1=sum(integrate(abs2(m1.E),m1.dΩ));
        integral=sum(integrate(m1.E*conj(m2.E),m2.dΩ));
        return integral/sqrt(integral1);
    end
end

"""
    overlap(m1::ScalarModeFEM,m2::ScalarfieldFEM)
"""
function overlap(m1::ScalarFieldFEM,m2::ScalarModeFEM)
    return conj(overlap(m2,m1))
end

"""
    overlap(m1::ScalarFieldFEM,m2::ScalarfieldFEM)
"""
function overlap(m1::ScalarFieldFEM,m2::ScalarFieldFEM)
    if (m1.Ω!=m2.Ω)
        throw(ArgumentError("The modes must have the same triangulation Ω"));
    else
        integral=sum(integrate(m1.E*conj(m2.E),m2.dΩ));
        return integral
    end
end

function overlap(m1::VectorModeFEM,m2::VectorModeFEM)
    if (m1.Ω!=m2.Ω)
        throw(ArgumentError("The modes must have the same triangulation Ω"));
    else
        integral1_task=Threads.@spawn 0.5*abs(sum(integrate(m1.Ex*conj(m1.Hy)-m1.Ey*conj(m1.Hx),m1.dΩ)));
        integral2_task=Threads.@spawn 0.5*abs(sum(integrate(m2.Ex*conj(m2.Hy)-m2.Ey*conj(m2.Hx),m2.dΩ)));
        integral_task=Threads.@spawn 0.5*sum(integrate(m1.Ex*conj(m2.Hy)-m1.Ey*conj(m2.Hx),m1.dΩ));
        integral1=fetch(integral1_task);
        integral2=fetch(integral2_task);
        integral=fetch(integral_task);
        return integral/sqrt(integral1*integral2);
    end
end

function overlap(m1::VectorModeFEM,m2::VectorFieldFEM)
    if (m1.Ω!=m2.Ω)
        throw(ArgumentError("The modes must have the same triangulation Ω"));
    else
        integral1_task=Threads.@spawn 0.5*abs(sum(integrate(m1.Ex*conj(m1.Hy)-m1.Ey*conj(m1.Hx),m1.dΩ)));
        integral_task=Threads.@spawn 0.5*sum(integrate(m1.Ex*conj(m2.Hy)-m1.Ey*conj(m2.Hx),m1.dΩ));
        integral1=fetch(integral1_task);
        integral=fetch(integral_task);
        return integral/sqrt(integral1);
    end
end

function overlap(m1::VectorFieldFEM,m2::VectorModeFEM)
    if (m1.Ω!=m2.Ω)
        throw(ArgumentError("The modes must have the same triangulation Ω"));
    else
        integral2_task=Threads.@spawn 0.5*abs(sum(integrate(m2.Ex*conj(m2.Hy)-m2.Ey*conj(m2.Hx),m2.dΩ)));
        integral_task =Threads.@spawn 0.5*sum(integrate(m1.Ex*conj(m2.Hy)-m1.Ey*conj(m2.Hx),m1.dΩ));
        integral2=fetch(integral2_task);
        integral=fetch(integral_task);
        return integral/sqrt(integral2);
    end
end

function overlap(m1::VectorFieldFEM,m2::VectorFieldFEM)
    if (m1.Ω!=m2.Ω)
        throw(ArgumentError("The modes must have the same triangulation Ω"));
    else
        integral=0.5*sum(integrate(m1.Ex*conj(m2.Hy)-m1.Ey*conj(m2.Hx),m1.dΩ));
        return integral
    end
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
    Aeff(m::VectorMode,n0::Union{Real,Function}=0)

n0 must be a constant or a function of x and y.

If n0=0, this function assumes that n0≃neff (correct for weakly-guiding fibers).
"""
function Aeff(m::VectorMode,n0::Union{Real,Function}=0)
    #Approximation (neff=n0)
    E2=trapz((m.x,m.y),m.Ex.*m.Hy-m.Ey.*m.Hx);
    if isa(n0,Function)
        if (first(methods(n0)).nargs!=3)
            throw(DomainError(n0, "The refractive index function n0 must have 2 argument (x,y)"));
        end
        E4_1=trapz((m.x,m.y),(abs2.(abs2.(m.Ex)+abs2.(m.Ey)+abs2.(m.Ez))).*(n0.(m.x',m.y).^2));
        E4_2=trapz((m.x,m.y),(abs2.(abs2.(m.Ex)+abs2.(m.Ey)-abs2.(m.Ez))).*(n0.(m.x',m.y).^2));
        return E2^2/E4_1*mu0/eps0,E2^2/E4_2*mu0/eps0;
    else
        E4_1=trapz((m.x,m.y),abs2.(abs2.(m.Ex)+abs2.(m.Ey)+abs2.(m.Ez)));
        E4_2=trapz((m.x,m.y),abs2.(abs2.(m.Ex)+abs2.(m.Ey)-abs2.(m.Ez)));
        if (n0==0)
            return E2^2/E4_1*mu0/eps0/(m.neff)^2,E2^2/E4_2*mu0/eps0/(m.neff)^2;
        else
            return E2^2/E4_1*mu0/eps0/n0^2,E2^2/E4_2*mu0/eps0/n0^2;
        end
    end
end

"""
    Aeff(m::ScalarModeFEM)

"""
function Aeff(m::ScalarModeFEM)
    E2_task=Threads.@spawn sum(integrate(abs2(m.E),m.dΩ));
    E4_task=Threads.@spawn sum(integrate(abs2(abs2(m.E)),m.dΩ));
    E2=fetch(E2_task);
    E4=fetch(E4_task);
    return E2^2/E4;
end

"""
    Aeff(m::VectorModeFEM,n0::Union{Real,Function}=0)

n0 must be a constant or a function of x and y.

If n0=0, this function assumes that n0≃neff (correct for weakly-guiding fibers).
"""
function Aeff(m::VectorModeFEM,n0::Union{Real,Function}=0)
    #Approximation (real(neff)=n0)
    if isa(n0,Function)
        if (first(methods(n0)).nargs!=3)
            throw(DomainError(n0, "The refractive index function n0 must have 2 argument (x,y)"));
        end
        n02=x->n0(x[1],x[2])^2;
        E2_task=Threads.@spawn sum(integrate(m.Ex*conj(m.Hy)-m.Ey*conj(m.Hx),m.dΩ));
        E4_1_task=Threads.@spawn sum(integrate(n02*abs2(abs2(m.Ex)+abs2(m.Ey)+abs2(m.Ez)),m.dΩ));
        E4_2_task=Threads.@spawn sum(integrate(n02*abs2(m.Ex*m.Ex+m.Ey*m.Ey+m.Ez*m.Ez),m.dΩ));
        E2=fetch(E2_task);
        E4_1=fetch(E4_1_task);
        E4_2=fetch(E4_2_task);
        return real(E2)^2/E4_1*mu0/eps0,real(E2)^2/E4_2*mu0/eps0;
    else
        E2_task=Threads.@spawn sum(integrate(m.Ex*conj(m.Hy)-m.Ey*conj(m.Hx),m.dΩ));
        E4_1_task=Threads.@spawn sum(integrate(abs2(abs2(m.Ex)+abs2(m.Ey)+abs2(m.Ez)),m.dΩ));
        E4_2_task=Threads.@spawn sum(integrate(abs2(m.Ex*m.Ex+m.Ey*m.Ey+m.Ez*m.Ez),m.dΩ));
        E2=fetch(E2_task);
        E4_1=fetch(E4_1_task);
        E4_2=fetch(E4_2_task);
        if (n0==0)
            return real(E2)^2/E4_1*mu0/eps0/(real(m.neff))^2,real(E2)^2/E4_2*mu0/eps0/(real(m.neff))^2;
        else
            return real(E2)^2/E4_1*mu0/eps0/n0^2,real(E2)^2/E4_2*mu0/eps0/n0^2;
        end
    end
end


############################# Nonlinear Coefficient #############################

"""
    nonLinearCoefficient(m::ScalarMode1D,n2::Union{Real,Function})

n2 must be a constant or a function of the radial coordinate r
"""
function nonLinearCoefficient(m::ScalarMode1D,n2::Union{Real,Function})
    if (isa(n2,Real))
        return n2*2*pi/m.lambda./Aeff(m);
    else
        if (first(methods(n2)).nargs!=2)
            throw(DomainError(n2, "The non-linear index function must have 1 argument"));
        end
        E2=trapz(m.r,2*pi*m.r.*(m.E).^2);
        E4=trapz(m.r,2*pi*(m.r.*(m.E).^4).*n2.(m.r));
        if (m.nu!=0)
            E2=E2*0.5;
            E4=E4*3.0/8.0;
        end
        k0=2*pi/m.lambda;
        return E4/(E2^2)*k0;
    end
end

"""
    nonLinearCoefficient(m::ScalarMode2D,n2::Union{Real,Function})

n2 must be a constant or a function of the cartesian coordinates x and y
"""
function nonLinearCoefficient(m::ScalarMode2D,n2::Union{Real,Function})
    if (isa(n2,Real))
        return n2*2*pi/m.lambda./Aeff(m);
    else
        if (first(methods(n2)).nargs!=3)
            throw(DomainError(n2, "The non-linear index function must have 2 argument"));
        end
        E2=trapz((m.x,m.y),(m.E).^2);
        E4=trapz((m.x,m.y),((m.E).^4).*n2.(m.x,m.y'));
        k0=2*pi/m.lambda;
        return E4/(E2^2)*k0;
    end
end

"""
    nonLinearCoefficient(m::VectorMode,n2::Union{Real,Function},n0::Union{Real,Function}=0)

n2 and n0 must be a constant or a function of x and y.

If n0=0, this function assumes that n0≃neff (correct for weakly-guiding fibers).
"""
function nonLinearCoefficient(m::VectorMode,n2::Union{Real,Function},n0::Union{Real,Function}=0)
    if (isa(n2,Real))
        return n2*2*pi/m.lambda./Aeff(m,n0);
    else
        if (first(methods(n2)).nargs!=3)
            throw(DomainError(n2, "The non-linear index function must have 2 argument"));
        end
        E2=trapz((m.x,m.y),m.Ex.*m.Hy-m.Ey.*m.Hx);
        k0=2*pi/m.lambda;
        if (isa(n0,Real))
            if (n0==0)
                nn=m.neff^2;
            else
                nn=n0^2;
            end
            E4_1=trapz((m.x,m.y),(abs2.(abs2.(m.Ex)+abs2.(m.Ey)+abs2.(m.Ez))).*n2.(m.x,m.y'));
            E4_2=trapz((m.x,m.y),(abs2.(abs2.(m.Ex)+abs2.(m.Ey)-abs2.(m.Ez))).*n2.(m.x,m.y'));
            return 1.0/(E2^2/E4_1*mu0/eps0)*k0*nn,1.0/(E2^2/E4_2*mu0/eps0)*k0*nn;
        else
            if (first(methods(n0)).nargs!=3)
                throw(DomainError(n0, "The refractive index n0 function must have 2 argument"));
            end
            E4_1=trapz((m.x,m.y),(abs2.(abs2.(m.Ex)+abs2.(m.Ey)+abs2.(m.Ez))).*(n0.(m.x,m.y').^2).*n2.(m.x,m.y'));
            E4_2=trapz((m.x,m.y),(abs2.(abs2.(m.Ex)+abs2.(m.Ey)-abs2.(m.Ez))).*(n0.(m.x,m.y').^2).*n2.(m.x,m.y'));
            return 1.0/(E2^2/E4_1*mu0/eps0)*k0,1.0/(E2^2/E4_2*mu0/eps0)*k0;
        end
    end
end

"""
    nonLinearCoefficient(m::ScalarModeFEM,n2::Union{Real,Function})

n2 must be a constant or a function of the cartesian coordinates x and y
"""
function nonLinearCoefficient(m::ScalarModeFEM,n2::Union{Real,Function})
    if (isa(n2,Real))
        return n2*2*pi/m.lambda./Aeff(m);
    else
        if (first(methods(n2)).nargs!=3)
            throw(DomainError(n2, "The non-linear index function must have 2 argument"));
        end
        n2bis=x->n2(x[1],x[2]);
        E2_task=Threads.@spawn sum(integrate(abs2(m.E),m.dΩ));
        E4_task=Threads.@spawn sum(integrate(n2bis*abs2(abs2(m.E)),m.dΩ));
        E2=fetch(E2_task);
        E4=fetch(E4_task);
        k0=2*pi/m.lambda;
        return E4/(E2^2)*k0;
    end
end

"""
    nonLinearCoefficient(m::VectorModeFEM,n2::Union{Real,Function},n0::Union{Real,Function}=0)

n2 and n0 must be constants or functions of x and y.

If n0=0, this function assumes that n0≃neff (correct for weakly-guiding fibers).
"""
function nonLinearCoefficient(m::VectorModeFEM,n2::Union{Real,Function},n0::Union{Real,Function}=0)
    if (isa(n2,Real))
        return n2*2*pi/m.lambda./Aeff(m,n0);
    else
        if (first(methods(n2)).nargs!=3)
            throw(DomainError(n2, "The non-linear index n2 function must have 2 argument"));
        end
        k0=2*pi/m.lambda;
        if (isa(n0,Real))
            if (n0==0)
                nn=m.neff^2;
            else
                nn=n0^2;
            end
            n2bis=x->n2(x[1],x[2]);
            E2_task=Threads.@spawn real(sum(integrate(m.Ex*conj(m.Hy)-m.Ey*conj(m.Hx),m.dΩ)));
            E4_1_task=Threads.@spawn sum(integrate(n2bis*abs2(abs2(m.Ex)+abs2(m.Ey)+abs2(m.Ez)),m.dΩ));
            E4_2_task=Threads.@spawn sum(integrate(n2bis*abs2(m.Ex*m.Ex+m.Ey*m.Ey+m.Ez*m.Ez),m.dΩ));
            E2=fetch(E2_task);
            E4_1=fetch(E4_1_task);
            E4_2=fetch(E4_2_task);
            return 1.0/(E2^2/E4_1*mu0/eps0)*k0*nn,1.0/(E2^2/E4_2*mu0/eps0)*k0*nn;
        else
            if (first(methods(n0)).nargs!=3)
                throw(DomainError(n0, "The refractive index n0 function must have 2 argument"));
            end
            n02n2=x->n2(x[1],x[2])*n0(x[1],x[2])^2;
            E2_task=Threads.@spawn real(sum(integrate(m.Ex*conj(m.Hy)-m.Ey*conj(m.Hx),m.dΩ)));
            E4_1_task=Threads.@spawn sum(integrate(n02n2*abs2(abs2(m.Ex)+abs2(m.Ey)+abs2(m.Ez)),m.dΩ));
            E4_2_task=Threads.@spawn sum(integrate(n02n2*abs2(m.Ex*m.Ex+m.Ey*m.Ey+m.Ez*m.Ez),m.dΩ));
            E2=fetch(E2_task);
            E4_1=fetch(E4_1_task);
            E4_2=fetch(E4_2_task);
            return 1.0/(E2^2/E4_1*mu0/eps0)*k0,1.0/(E2^2/E4_2*mu0/eps0)*k0;
        end
    end
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
    MFD(m::ScalarMode2D,theta::Real=0)
"""
function MFD(m::ScalarMode2D,theta::Real=0)
    rmax=sqrt(maximum(m.x)^2+maximum(m.y)^2);
    dr=min(minimum(diff(m.x)),minimum(diff(m.y)));
    nb=Integer(2*ceil(rmax/dr)+1);
    r=collect(LinRange(-rmax,rmax,nb));
    x=r*cosd(theta);
    y=r*sind(theta);
    interp=LinearInterpolation((m.x,m.y),m.E,extrapolation_bc=0);
    E=interp.(x,y);
    posmax=argmax(abs.(E));
    pos=argmin(abs.(abs.(E).-(maximum(abs.(E))/exp(1))));
    return abs(2*(r[pos]-r[posmax]));
end

"""
    MFD(m::VectorMode,theta::Real=0)
"""
function MFD(m::VectorMode,theta::Real=0)
    rmax=sqrt(maximum(m.x)^2+maximum(m.y)^2);
    dr=min(minimum(diff(m.x)),minimum(diff(m.y)));
    nb=Integer(2*ceil(rmax/dr)+1);
    r=collect(LinRange(-rmax,rmax,nb));
    x=r*cosd(theta);
    y=r*sind(theta);
    Pz=m.Ex.*m.Hy-m.Ey.*m.Hx;
    interp=LinearInterpolation((m.x,m.y),Pz,extrapolation_bc=0);
    E2=interp.(x,y);
    posmax=argmax(abs.(E2));
    pos=argmin(abs.(abs.(E2).-(maximum(abs.(E2))/exp(2))));
    return abs(2*(r[pos]-r[posmax]));
end

"""
    MFD(m::ScalarModeFEM,theta::Real=0)
"""
function MFD(m::ScalarModeFEM,theta::Real=0)
    x=getindex.(m.Ω.grid.node_coordinates,1)
    y=getindex.(m.Ω.grid.node_coordinates,2)
    rmax=maximum(x*cosd(theta)+y*sind(theta))
    rmin=minimum(x*cosd(theta)+y*sind(theta))
    r=LinRange(rmin,rmax,10000)
    E=zeros(10000);
    Threads.@threads for i=1:10000
        try
            E[i]=m.E(Point(r[i]*cosd(theta),r[i]*sind(theta)));
        catch e
            E[i]=0.0
        end
    end
    posmax=argmax(abs.(E));
    pos=argmin(abs.(abs.(E).-(maximum(abs.(E))/exp(1))));
    return abs(2*(r[pos]-r[posmax]));
end

"""
    MFD(m::VectorModeFEM,theta::Real=0)
"""
function MFD(m::VectorModeFEM,theta::Real=0)
    x=getindex.(m.Ω.grid.node_coordinates,1)
    y=getindex.(m.Ω.grid.node_coordinates,2)
    rmax=maximum(x*cosd(theta)+y*sind(theta))
    rmin=minimum(x*cosd(theta)+y*sind(theta))
    r=LinRange(rmin,rmax,10000)
    Pzf=real(m.Ex*conj(m.Hy)-m.Ey*conj(m.Hx));
    Pz=zeros(10000);
    Threads.@threads for i=1:10000
        try
            Pz[i]=Pzf(Point(r[i]*cosd(theta),r[i]*sind(theta)));
        catch e
            Pz[i]=0.0
        end
    end
    posmax=argmax(abs.(Pz));
    pos=argmin(abs.(abs.(Pz).-(maximum(abs.(Pz))/exp(2))));
    return abs(2*(r[pos]-r[posmax]));
end



