const realLength=Unitful.Quantity{<:Real,Unitful.𝐋}
const inverseRealLength=Unitful.Quantity{<:Real,Unitful.dimension(u"m^-1")}
const ElectricField=Unitful.Quantity{<:Number,Unitful.dimension(u"V/m")}
const MagneticField=Unitful.Quantity{<:Number,Unitful.dimension(u"A/m")}
const ElectricFieldDim=Unitful.dimension(u"V/m")
const MagneticFieldDim=Unitful.dimension(u"A/m")

"""
    abstract type FiberEMField end

Described an electromagnetic field in a fiber.
"""
abstract type FiberEMField end
"""
    abstract type ScalarFiberEMField <:FiberEMField end

Described an electromagnetic field in a fiber in the scalar approximation.
"""
abstract type ScalarFiberEMField <:FiberEMField end

Base.:transpose(f::FiberEMField) = f
Broadcast.:broadcastable(f::FiberEMField)=Ref(f)

"""
    struct ScalarFiberEMField1D <:ScalarFiberEMField

Described an electromagnetic field in a fiber with a cylindrical symmetry in the scalar approximation.
- nu :: `Int` - Azimuthal number
- E :: `UnitfulField{1,(),ElectricFieldDim}`

The additional property `normE` gives the norm of the electric field.
"""
struct ScalarFiberEMField1D <:ScalarFiberEMField
    nu::Int
    E::UnitfulField{1,(),ElectricFieldDim}
    function ScalarFiberEMField1D(nu,E::UnitfulField{1,(),ElectricFieldDim})
        if isa(E,FEMField)
            new(nu,E)
        else
            new(nu,deepcopy(E))
        end
    end
end

"""
    struct ScalarFiberEMField2D <:ScalarFiberEMField
    
Described an electromagnetic field in a fiber with a cylindrical symmetry in the scalar approximation.
- E :: `UnitfulField{2,(),ElectricFieldDim}`

The additional property `normE` gives the norm of the electric field.
"""
struct ScalarFiberEMField2D <:ScalarFiberEMField
    E::UnitfulField{2,(),ElectricFieldDim}
    function ScalarFiberEMField2D(E::UnitfulField{2,(),ElectricFieldDim})
        if isa(E,FEMField)
            new(E)
        else
            new(deepcopy(E))
        end
    end
    function ScalarFiberEMField2D(f::ScalarFiberEMField1D,orientation_angle::Real)
            new(FunctionField(2,x->f.E(Point(hypot(x[1],x[2]),))*cosd(f.nu*rad2deg(atan(x[2],x[1]))+orientation_angle)))
    end
end

#permettre des opérations avec des champs sans dimension?
Base.:*(f::ScalarFiberEMField1D,k::Number) = (unit(k)==NoUnits) ? ScalarFiberEMField1D(f.nu,k*f.E) : throw(ArgumentError("k cannot have a unit"));
Base.:*(k::Number,f::ScalarFiberEMField1D) = (unit(k)==NoUnits) ? ScalarFiberEMField1D(f.nu,k*f.E) : throw(ArgumentError("k cannot have a unit"));
Base.:/(f::ScalarFiberEMField1D,k::Number) = (unit(k)==NoUnits) ? ScalarFiberEMField1D(f.nu,f.E/k) : throw(ArgumentError("k cannot have a unit"));
Base.:*(f::ScalarFiberEMField2D,k::Number) = (unit(k)==NoUnits) ? ScalarFiberEMField2D(k*f.E) : throw(ArgumentError("k cannot have a unit"));
Base.:*(k::Number,f::ScalarFiberEMField2D) = (unit(k)==NoUnits) ? ScalarFiberEMField2D(k*f.E) : throw(ArgumentError("k cannot have a unit"));
Base.:/(f::ScalarFiberEMField2D,k::Number) = (unit(k)==NoUnits) ? ScalarFiberEMField2D(f.E/k) : throw(ArgumentError("k cannot have a unit"));
for op in (:+,:-)
    @eval begin
        function ($op)(f1::ScalarFiberEMField1D,f2::ScalarFiberEMField1D)
            if (f1.nu==f2.nu) 
                ScalarFiberEMField1D(f1.nu,($op)(f1.E,f2.E));
            else
                throw(ArgumentError("Both EMField must have the same nu"))
            end
        end
        function ($op)(f1::ScalarFiberEMField2D,f2::ScalarFiberEMField2D)
            ScalarFiberEMField2D(($op)(f1.E,f2.E));
        end
    end 
end

for op in (:+,:-,:real,:imag,:conj)
    @eval begin
        function ($op)(f::ScalarFiberEMField1D)
            ScalarFiberEMField1D(f.nu,($op)(f.E));
        end
        function ($op)(f::ScalarFiberEMField2D)
            ScalarFiberEMField2D(($op)(f.E));
        end
    end 
end

Base.:abs(f::ScalarFiberEMField2D) = ScalarFiberEMField(abs(f.E)) #norm?


function Base.show(io::IO, mime::MIME"text/plain", f::ScalarFiberEMField1D)
    print(io,typeof(f),"\n");
    print(io," - nu = ",f.nu);
end

function Base.show(io::IO, mime::MIME"text/plain", f::ScalarFiberEMField2D)
    print(io,typeof(f),"\n");
end

function Base.show(io::IO, f::ScalarFiberEMField1D)
    print(io,typeof(f));
    print(io,"(nu = ",f.nu,")");
end

function Base.show(io::IO, f::ScalarFiberEMField2D)
    print(io,typeof(f));
end

"""
    struct VectorFiberEMField <: FiberEMField
    
Described an electromagnetic field in a fiber.
- E :: `UnitfulField{2,(3,),ElectricFieldDim}`
- H :: `UnitfulField{2,(3,),ElectricFieldDim}`

Additional properties can be computed:
- normE : Norm of the electric field
- Ex, Ey, Ez : Electric field components
- normH : Norm of the magnetic field
- Hx, Hy, Hz : Magnetic field components
- P : Poynting vector
- Px, Py, Pz : Poynting vector components
"""
struct VectorFiberEMField <: FiberEMField
    E::UnitfulField{2,(3,),ElectricFieldDim}
    H::UnitfulField{2,(3,),MagneticFieldDim}
    function VectorFiberEMField(E::UnitfulField{2,(3,),ElectricFieldDim},H::UnitfulField{2,(3,),MagneticFieldDim})
        new(E,H)
    end
end

Base.:*(f::VectorFiberEMField,k::Number) = VectorFiberEMField(k*f.E,k*f.H);
Base.:*(k::Number,f::VectorFiberEMField) = VectorFiberEMField(k*f.E,k*f.H);
Base.:/(f::VectorFiberEMField,k::Number) = VectorFiberEMField(f.E/k,f.H/k);

for op in (:+,:-)
    @eval begin
        function ($op)(f1::VectorFiberEMField,f2::VectorFiberEMField) 
            VectorFiberEMField(($op)(f1.E,f2.E),($op)(f1.H,f2.H));
        end
    end 
end

for op in (:+,:-,:real,:imag,:conj)
    @eval begin
        function ($op)(f::VectorFiberEMField)
            VectorFiberEMField(($op)(f.E),($op)(f.H));
        end
    end 
end

function Base.show(io::IO, mime::MIME"text/plain", f::VectorFiberEMField)
    print(io,typeof(f));
end

function Base.show(io::IO, f::VectorFiberEMField)
    print(io,typeof(f));
end


Base.propertynames(::ScalarFiberEMField1D) =(:nu,:E,:normE)

function Base.getproperty(f::ScalarFiberEMField1D,s::Symbol)
    if s==:normE
        return abs(f.E)
    else
        return getfield(f,s)
    end
end


Base.propertynames(::ScalarFiberEMField2D) =(:E,:normE)

function Base.getproperty(f::ScalarFiberEMField2D,s::Symbol)
    if s==:normE
        return abs(f.E)
    else
        return getfield(f,s)
    end
end


Base.propertynames(::VectorFiberEMField) =(:E,:H,:Ex,:Ey,:Ez,:normE,:Hx,:Hy,:Hz,:normH,:P,:Px,:Py,:Pz,:normP)

function Base.getproperty(f::VectorFiberEMField,s::Symbol)#A optimiser
    if s==:Ex
        return dot(VectorValue(1,0,0),f.E)
    elseif s==:Ey
        return dot(VectorValue(0,1,0),f.E)
    elseif s==:Ez
        return dot(VectorValue(0,0,1),f.E)
    elseif s==:Hx
        return dot(VectorValue(1,0,0),f.H)
    elseif s==:Hy
        return dot(VectorValue(0,1,0),f.H)
    elseif s==:Hz
        return dot(VectorValue(0,0,1),f.H)
    elseif s==:normE    
        return sqrt(real(dot(f.E,conj(f.E))))
    elseif s==:normH
        return sqrt(real(dot(f.H,conj(f.H))))
    elseif s==:P
        return 0.5*real(cross(f.E,conj(f.H)))
    elseif s==:Px
        return dot(VectorValue(1,0,0),0.5*real(cross(f.E,conj(f.H))))
    elseif s==:Py
        return dot(VectorValue(0,1,0),0.5*real(cross(f.E,conj(f.H))))
    elseif s==:Pz
        return dot(VectorValue(0,0,1),0.5*real(cross(f.E,conj(f.H))))
    elseif s==:normP
        return norm(0.5*real(cross(f.E,conj(f.H))))
    else
        return getfield(f,s)
    end
end


"""
    struct Mode{T<:Union{Field,Nothing}}

- Name :: `String`
- neff :: `Number`
- lambda :: `realLength`
- field :: `FiberEMField` or `Nothing`

Additional properties can be computed:
- beta : Propagation constant 
- alpha : Attenuation constant
- losses : Losses in dB/km
"""
struct Mode{T<:Union{FiberEMField,Nothing}}
    Name::String
    neff::Number
    lambda::realLength
    EMField::T
    Mode(Name,neff,lambda,field)=isnothing(field) ? new{Nothing}(Name,neff,lambda,nothing) : new{typeof(field)}(Name,neff,lambda,field)
    Mode(Name,neff,lambda)=new{Nothing}(Name,neff,lambda,nothing)
end

"""
    convertTo2D(f::Union{ScalarFiberEMField1D,Mode{ScalarFiberEMField1D}},orientation_angle::Real=0)
"""
function convertTo2D(f::ScalarFiberEMField1D,orientation_angle::Real=0)
    return ScalarFiberEMField2D(FunctionField(2,x->f.E(Point(hypot(x[1],x[2]),))*cosd(f.nu*rad2deg(atan(x[2],x[1]))-orientation_angle)))
end

function convertTo2D(m::Mode{ScalarFiberEMField1D},orientation_angle::Real=0)
    return Mode(m.Name,m.neff,m.lambda,ScalarFiberEMField2D(FunctionField(2,x->m.EMField.E(Point(hypot(x[1],x[2]),))*cosd(m.EMField.nu*rad2deg(atan(x[2],x[1]))-orientation_angle))))
end

"""
    convertToVector(m::Mode{ScalarFiberEMField2D},polarization_angle::Real=0)
"""
function convertToVector(m::Mode{ScalarFiberEMField2D},polarization_angle::Real=0)
    E=m.EMField.E*VectorValue(cosd(polarization_angle),sind(polarization_angle),0)
    H=m.EMField.E*VectorValue(-sind(polarization_angle)*m.neff/c_Unitful/mu0_Unitful,cosd(polarization_angle)*m.neff/c_Unitful/mu0_Unitful,0u"s/H")
    return Mode(m.Name,m.neff,m.lambda,VectorFiberEMField(E,H))
end

"""
    convertToVector(m::Mode{ScalarFiberEMField1D},orientation_angle::Real=0,polarization_angle::Real=0)
"""
function convertToVector(m::Mode{ScalarFiberEMField1D},orientation_angle::Real=0,polarization_angle::Real=0)
    convertToVector(convertTo2D(m,orientation_angle),polarization_angle)
end

Base.propertynames(::Mode) =(:Name,:neff,:lambda,:EMField,:beta,:alpha,:losses)

function Base.getproperty(m::Mode,s::Symbol)
    if s==:losses
        return ustrip(u"m^-1",4*pi/m.lambda*imag(m.neff))*10/log(10)*1000*u"dB/km"
    elseif s==:beta
        return uconvert(u"m^-1",2*pi/m.lambda*m.neff)
    elseif s==:alpha
        return uconvert(u"m^-1",4*pi/m.lambda*imag(m.neff))
    else
        return getfield(m,s)
    end
end

Base.:transpose(m::Mode) = m
Broadcast.:broadcastable(m::Mode)=Ref(m)

function Base.show(io::IO, ::MIME"text/plain",m::Mode) 
    print(io,"Name = ",m.Name,"\nneff = ",m.neff,"\nlambda = ",m.lambda,"\nEMfield = ",typeof(m.EMField));
end
Base.show(io::IO, m::Mode) = print(io,"[",m.Name," - neff = ",m.neff," - λ = ",m.lambda," , EMField = ",typeof(m.EMField),"]") 


############################# GetField->Propagation ######################
"""
    getEMField(m::Mode,z::Real=0)
"""
function getEMField(m::Mode,z::realLength=0u"m")
    if (z==0u"m")
        return m.EMField
    else
        return m.EMField*exp(im*z*2*pi/m.lambda*m.neff)
    end
end

"""
    Aeff(m::Mode)
"""
function Aeff(::Mode{Nothing})
    throw(ArgumentError("The mode does not contain any EMField"))
end

function Aeff(m::Mode{ScalarFiberEMField1D};rtol::Number=1E-6,initdiv::Int=1)
    E2=integrate(abs2(m.EMField.E)*(x->x[1]),[0]*unit(m.lambda),[Inf]*unit(m.lambda),characteristic_length=[m.lambda],rtol=rtol,initdiv=initdiv)*2*pi
    E4=integrate(abs2(abs2(m.EMField.E))*(x->x[1]),[0]*unit(m.lambda),[Inf]*unit(m.lambda),characteristic_length=[m.lambda],rtol=rtol,initdiv=initdiv)*2*pi
    if (m.EMField.nu!=0)
        E2=E2*0.5;
        E4=E4*3.0/8.0;
    end
    return E2^2/E4;
end

function Aeff(m::Mode{ScalarFiberEMField2D};rtol::Number=1E-4,initdiv::Int=1)
    E2=integrate(abs2(m.EMField.E),characteristic_length=[m.lambda,m.lambda],rtol=rtol,initdiv=initdiv)
    E4=integrate(abs2(abs2(m.EMField.E)),characteristic_length=[m.lambda,m.lambda],rtol=rtol,initdiv=initdiv)
    return E2^2/E4;
end

function Aeff(m::Mode{VectorFiberEMField},n0::Union{Real,Function}=0;rtol::Number=1E-4,initdiv::Int=1)
    if (n0==0)
        n0=real(m.neff)
    end
    E2_task=Threads.@spawn integrate(m.EMField.Pz,characteristic_length=[m.lambda,m.lambda],rtol=rtol,initdiv=initdiv)*2
    E4_1_task=Threads.@spawn integrate(abs2(dot(m.EMField.E,conj(m.EMField.E)))*n0*n0,characteristic_length=[m.lambda,m.lambda],rtol=rtol,initdiv=initdiv)
    E4_2_task=Threads.@spawn integrate(abs2(dot(m.EMField.E,m.EMField.E))*n0*n0,characteristic_length=[m.lambda,m.lambda],rtol=rtol,initdiv=initdiv)
    E22=fetch(E2_task)^2*mu0_Unitful/eps0_Unitful
    return uconvert(unit(m.lambda)^2,E22/fetch(E4_1_task)),uconvert(unit(m.lambda^2),E22/fetch(E4_2_task))
end

function normalize(m::Mode{Nothing})
    return m
end

"""
    normalize(m::Mode{ScalarFiberEMField};unitIntegral::Bool=true)
"""
function normalize(m::Mode{ScalarFiberEMField1D};unitIntegral::Bool=true,rtol::Number=1E-6,initdiv::Int=1)
    E2=integrate(abs2(m.EMField.E)*(x->x[1]),[0]*unit(m.lambda),[Inf]*unit(m.lambda),characteristic_length=[m.lambda],rtol=rtol,initdiv=initdiv)*2*pi
    if (m.EMField.nu!=0)
        E2=E2*0.5;
    end
    if unitIntegral
        E=m.EMField.E/ustrip(u"V",sqrt(E2))
    else
        E=m.EMField.E*ustrip(u"W^(-1/2)",sqrt(2*mu0_Unitful*c_Unitful/real(m.neff)/E2));
    end
    return Mode(m.Name,m.neff,m.lambda,ScalarFiberEMField1D(m.EMField.nu,E))
end

function normalize(m::Mode{ScalarFiberEMField2D};unitIntegral::Bool=true,rtol::Number=1E-4,initdiv::Int=1)
    E2=integrate(abs2(m.EMField.E),characteristic_length=[m.lambda,m.lambda],rtol=rtol,initdiv=initdiv)
    if unitIntegral
        E=m.EMField.E/ustrip(u"V",sqrt(E2))
    else
        E=m.EMField.E*ustrip(u"W^(-1/2)",sqrt(2*mu0_Unitful*c_Unitful/real(m.neff)/E2));
    end
    return Mode(m.Name,m.neff,m.lambda,ScalarFiberEMField2D(E))
end

"""
    normalize(m::Mode{VectorFiberEMField})
"""
function normalize(m::Mode{VectorFiberEMField};rtol::Number=1E-4,initdiv::Int=1)
    E2=real(integrate(m.EMField.Pz,characteristic_length=[m.lambda,m.lambda],rtol=rtol,initdiv=initdiv));
    E=m.EMField.E/ustrip(u"sqrt(W)",sqrt(E2))
    H=m.EMField.H/ustrip(u"sqrt(W)",sqrt(E2))
    return Mode(m.Name,m.neff,m.lambda,VectorFiberEMField(E,H))
end


"""
    overlap(f::ScalarFiberEMField1D,m::Mode{ScalarFiberEMField1D})
"""
function overlap(f::ScalarFiberEMField1D,m::Mode{ScalarFiberEMField1D};rtol::Number=1E-6,atol::Number=0u"V",initdiv::Int=1)
    if (f.nu!= m.EMField.nu)
        return zero(ComplexF64)u"V"
    end
    m2=integrate(abs2(m.EMField.E)*(x->x[1]),[0]*unit(m.lambda),[Inf]*unit(m.lambda),characteristic_length=[m.lambda],rtol=rtol,initdiv=initdiv)*2*pi
    if atol==0u"V"
        atol=m2*rtol/2/pi
    else
        atol=atol*sqrt(m2)/2/pi
    end
    m1=integrate(f.E*conj(m.EMField.E)*(x->x[1]),[0]*unit(m.lambda),[Inf]*unit(m.lambda),characteristic_length=[m.lambda],rtol=rtol,atol=atol,intidiv=initdiv)*2*pi
    result=m1/sqrt(m2)
    return ComplexF64(ustrip(u"V",result))u"V"
end

"""
    overlap(m1::Mode{ScalarFiberEMField2D},m2::Mode{ScalarFiberEMField2D})
"""
function overlap(m1::Mode{ScalarFiberEMField2D},m2::Mode{ScalarFiberEMField2D};rtol::Number=1E-4,atol::Number=0,initdiv::Int=1)
    int2=integrate(abs2(m2.EMField.E),characteristic_length=[m2.lambda,m2.lambda],rtol=rtol,initdiv=initdiv)
    int1=integrate(abs2(m1.EMField.E),characteristic_length=[m1.lambda,m1.lambda],rtol=rtol,initdiv=initdiv)
    if atol==0
        atol=min(int2,int1)*rtol
    else
        atol=atol*sqrt(int2*int1)
    end
    int=integrate(m1.EMField.E*conj(m2.EMField.E),characteristic_length=[min(m1.lambda,m2.lambda),min(m1.lambda,m2.lambda)],rtol=rtol,atol=atol,intidiv=initdiv)
    result=int/sqrt(int2*int1)
    return ComplexF64(ustrip(NoUnits,result))
end

"""
    overlap(m::Mode{ScalarFiberEMField1D},f::ScalarFiberEMField1D)
"""
function overlap(m::Mode{ScalarFiberEMField1D},f::ScalarFiberEMField1D;rtol::Number=1E-6,atol::Number=0u"V",initdiv::Int=1)
    return conj(overlap(f,m,rtol=rtol,atol=atol,initdiv=initdiv))
end


"""
    overlap(function overlap(m1::Mode{ScalarFiberEMField1D},m2::Mode{ScalarFiberEMField1D})
"""
function overlap(m1::Mode{ScalarFiberEMField1D},m2::Mode{ScalarFiberEMField1D};rtol::Number=1E-6,atol::Number=0,initdiv::Int=1)
    if (m1.EMField.nu!= m2.EMField.nu)
        return zero(ComplexF64)
    end
    int2=integrate(abs2(m2.EMField.E)*(x->x[1]),[0]*unit(m2.lambda),[Inf]*unit(m2.lambda),characteristic_length=[m2.lambda],rtol=rtol,initdiv=initdiv)*2*pi
    int1=integrate(abs2(m1.EMField.E)*(x->x[1]),[0]*unit(m1.lambda),[Inf]*unit(m1.lambda),characteristic_length=[m1.lambda],rtol=rtol,initdiv=initdiv)*2*pi
    if atol==0
        atol=min(int2,int1)*rtol/2/pi
    else
        atol=atol*sqrt(int2*int1)/2/pi
    end
    int=integrate(m1.EMField.E*conj(m2.EMField.E)*(x->x[1]),[0]*unit(m1.lambda),[Inf]*unit(m1.lambda),characteristic_length=[min(m1.lambda,m2.lambda)],rtol=rtol,atol=atol,intidiv=initdiv)*2*pi
    result=int/sqrt(int2*int1)
    return ComplexF64(ustrip(NoUnits,result))
end

"""
    overlap(f::ScalarFiberEMField2D,m::Mode{ScalarFiberEMField2D})
"""
function overlap(f::ScalarFiberEMField2D,m::Mode{ScalarFiberEMField2D};rtol::Number=1E-4,atol::Number=0u"V",initdiv::Int=1)
    int2=integrate(abs2(m.EMField.E),characteristic_length=[m.lambda],rtol=rtol,initdiv=initdiv)
    if atol==0u"V"
        atol=int2*rtol
    else
        atol=atol*sqrt(int2)
    end
    int=integrate(f.E.*conj(m.EMField.E),characteristic_length=[m.lambda,m.lambda],rtol=rtol,atol=atol,intidiv=initdiv)
    result=int/sqrt(int2)
    return ComplexF64(ustrip(u"V",result))u"V"
end

"""
    overlap(m::Mode{ScalarFiberEMField2D},f::ScalarFiberEMField2D)
"""
function overlap(m::Mode{ScalarFiberEMField2D},f::ScalarFiberEMField2D;rtol::Number=1E-4,atol::Number=0u"V",initdiv::Int=1)
    return conj(overlap(f,m,rtol=rtol,atol=atol,initdiv=initdiv))
end

"""
    overlap(m1::Mode{VectorFiberEMField},m2::Mode{VectorFiberEMField})
"""
function overlap(m1::Mode{VectorFiberEMField},m2::Mode{VectorFiberEMField};rtol::Number=1E-4,atol::Number=0,initdiv::Int=1)
    int2_task=Threads.@spawn real(integrate(m2.EMField.Pz,characteristic_length=[m2.lambda,m2.lambda],rtol=rtol,initdiv=initdiv))
    int1_task=Threads.@spawn real(integrate(m1.EMField.Pz,characteristic_length=[m1.lambda,m1.lambda],rtol=rtol,initdiv=initdiv))
    int2=fetch(int2_task)
    int1=fetch(int1_task)
    if atol==0
        atol=min(int2,int1)*2*rtol
    else
        atol=atol*2*sqrt(int2*int1)
    end
    int=integrate(m1.EMField.Ex*conj(m2.EMField.Hy)-m1.EMField.Ey*conj(m2.EMField.Hx),characteristic_length=[min(m1.lambda,m2.lambda),min(m1.lambda,m2.lambda)],rtol=rtol,atol=atol,intidiv=initdiv)*0.5
    result=int/sqrt(int2*int1)
    return ComplexF64(ustrip(NoUnits,result))
end

"""
    overlap(m::Mode{VectorFiberEMField},f::VectorFiberEMField)
"""
function overlap(m1::Mode{VectorFiberEMField},f::VectorFiberEMField;rtol::Number=1E-4,atol::Number=0u"sqrt(W)",initdiv::Int=1)
    int1=real(integrate(m1.EMField.Pz,characteristic_length=[m1.lambda,m1.lambda],rtol=rtol,initdiv=initdiv))
    if atol==0u"sqrt(W)"
        atol=int1*2*rtol
    else
        atol=atol*2*sqrt(int1)
    end
    int=integrate(m1.EMField.Ex*conj(f.Hy)-m1.EMField.Ey*conj(f.Hx),characteristic_length=[m1.lambda,m1.lambda],rtol=rtol,atol=atol,intidiv=initdiv)*0.5
    result=int/sqrt(int1)
    return ComplexF64(ustrip(u"sqrt(W)",result))*u"sqrt(W)"
end

"""
    overlap(f::VectorFiberEMField,m::Mode{VectorFiberEMField})
"""
function overlap(f::VectorFiberEMField,m1::Mode{VectorFiberEMField};rtol::Number=1E-4,atol::Number=0u"sqrt(W)",initdiv::Int=1)
    int1=real(integrate(m1.EMField.Pz,characteristic_length=[m1.lambda,m1.lambda],rtol=rtol,initdiv=initdiv))
    if atol==0u"sqrt(W)"
        atol=int1*2*rtol
    else
        atol=atol*2*sqrt(int1)
    end
    int=integrate(f.Ex*conj(m1.EMField.Hy)-f.Ey*conj(m1.EMField.Hx),characteristic_length=[m1.lambda,m1.lambda],rtol=rtol,atol=atol,intidiv=initdiv)*0.5
    result=int/sqrt(int1)
    return ComplexF64(ustrip(u"sqrt(W)",result))*u"sqrt(W)"
end

"""
    nonLinearCoefficient(m::Mode{ScalarFiberEMField},n2::Union{Unitful.Quantity{<:Number,Unitful.dimension(u"m^2/W")},Function})
"""
function nonLinearCoefficient(m::Mode{ScalarFiberEMField1D},n2::Union{Unitful.Quantity{<:Number,Unitful.dimension(u"m^2/W")},Function};rtol::Number=1E-6,initdiv::Int=1)
    if (isa(n2,Quantity))
        return uconvert(u"W^-1*m^-1",n2*2*pi/m.lambda/Aeff(m,rtol=rtol,initdiv=initdiv));
    else
        x=n2(Point(0.0,)u"m")
        if dimension(x)!=Unitful.dimension(u"m^2/W")
            throw(DomainError(n2, "The non-linear index function must return a quantity with dimension m^2/W"));
        end
        E2=integrate(abs2(m.EMField.E)*(x->x[1]),[0]*unit(m.lambda),[Inf]*unit(m.lambda),characteristic_length=[m.lambda],rtol=rtol,intidiv=initdiv)*2*pi
        E4=integrate(abs2(abs2(m.EMField.E))*(x->x[1])*n2,[0]*unit(m.lambda),[Inf]*unit(m.lambda),characteristic_length=[m.lambda],rtol=rtol,initdiv=initdiv)*2*pi
        if (m.EMField.nu!=0)
            E2=E2*0.5;
            E4=E4*3.0/8.0;
        end
        k0=2*pi/m.lambda;
        return uconvert(u"W^-1*m^-1",E4/(E2^2)*k0);
    end
end

function nonLinearCoefficient(m::Mode{ScalarFiberEMField2D},n2::Union{Unitful.Quantity{<:Number,Unitful.dimension(u"m^2/W")},Function};rtol::Number=1E-4,initdiv::Int=1)
    if (isa(n2,Quantity))
        return uconvert(u"W^-1*m^-1",n2*2*pi/m.lambda/Aeff(m,rtol=rtol,initdiv=initdiv));
    else
        x=n2(Point(0.0,0.0)u"m")
        if dimension(x)!=Unitful.dimension(u"m^2/W")
            throw(DomainError(n2, "The non-linear index function must return a quantity with dimension m^2/W"));
        end
        E2=integrate(abs2(m.EMField.E),characteristic_length=[m.lambda,m.lambda],rtol=rtol,initdiv=initdiv)
        E4=integrate(abs2(abs2(m.EMField.E))*n2,characteristic_length=[m.lambda,m.lambda],rtol=rtol,initdiv=initdiv)
        k0=2*pi/m.lambda;
        return uconvert(u"W^-1*m^-1",E4/(E2^2)*k0);
    end
end

"""
    nonLinearCoefficient(m::Mode{VectorFiberEMField},n2::Union{Unitful.Quantity{<:Number,Unitful.dimension(u"m^2/W")},Function},n0::Union{Real,Function}=0)
"""
function nonLinearCoefficient(m::Mode{VectorFiberEMField},n2::Union{Unitful.Quantity{<:Number,Unitful.dimension(u"m^2/W")},Function},n0::Union{Real,Function}=0;rtol::Number=1E-4,initdiv::Int=1)
    if (isa(n2,Quantity))
        return uconvert.(u"m^-1*W^-1",(n2*2*pi/m.lambda./Aeff(m,n0,rtol=rtol,initdiv=initdiv)));
    else
        x=n2(Point(0.0,0.0)u"m")
        if dimension(x)!=Unitful.dimension(u"m^2/W")
            throw(DomainError(n2, "The non-linear index function must return a quantity with dimension m^2/W"));
        end
        k0=2*pi/m.lambda;
        if (isa(n0,Real))
            if (n0==0)
                n02n2=x->n2(x)*real(m.neff)^2;
            else
                n02n2=x->n2(x)*n0^2;
            end
        else
            x=n0(Point(0.0,0.0)u"m")
            if dimension(x)!=Unitful.NoDims
                throw(DomainError(n0, "The refractive index function must return a quantity with no dimension"));
            end
            n02n2=x->n2(x)*n0(x)^2;
        end
        E2_task=Threads.@spawn integrate(m.EMField.Pz,characteristic_length=[m.lambda,m.lambda],rtol=rtol,initdiv=initdiv)
        E4_1_task=Threads.@spawn integrate(abs2(dot(m.EMField.E,conj(m.EMField.E)))*n02n2,characteristic_length=[m.lambda,m.lambda],rtol=rtol,initdiv=initdiv)
        E4_2_task=Threads.@spawn integrate(abs2(dot(m.EMField.E,m.EMField.E))*n02n2,characteristic_length=[m.lambda,m.lambda],rtol=rtol,initdiv=initdiv)
        E2=fetch(E2_task);
        E4_1=fetch(E4_1_task);
        E4_2=fetch(E4_2_task);
        k0=2*pi/m.lambda;
        return uconvert.(u"m^-1*W^-1",(1.0/(E2^2/E4_1*mu0_Unitful/eps0_Unitful)*k0,1.0/(E2^2/E4_2*mu0_Unitful/eps0_Unitful)*k0));
    end
end

"""
    MFD(m::Mode{ScalarFiberEMField1D})
"""
function MFD(m::Mode{ScalarFiberEMField1D})
    f=x->m.EMField.E(x*m.lambda)-m.EMField.E(Point(0.0*m.lambda,))/ℯ
    return 2*abs(find_zero(f,0))*m.lambda
end

"""
    MFD(m::Mode{ScalarFiberEMField2D},theta::Real=0)
"""
function MFD(m::Mode{ScalarFiberEMField2D},theta::Real=0)
    f=x->m.EMField.E(Point(x*m.lambda*cosd(theta),x*m.lambda*sind(theta)))-m.EMField.E(Point(0.0*m.lambda,0.0*m.lambda))/ℯ
    return 2*abs(find_zero(f,0))*m.lambda
end

"""
    MFD(m::Mode{VectorFiberEMField},theta::Real=0)
"""
function MFD(m::Mode{VectorFiberEMField},theta::Real=0)
    f=x->m.EMField.Pz(Point(x*m.lambda*cosd(theta),x*m.lambda*sind(theta)))-m.EMField.Pz(Point(0.0*m.lambda,0.0*m.lambda))/ℯ^2
    return 2*abs(find_zero(f,0.0))*m.lambda
end

"""
    Gridap.:writevtk(name::String,f::Union{FiberEMField,Mode))
"""
function Gridap.:writevtk(name::String,f::ScalarFiberEMField)
    if isa(f.E,FEMField)
        Ω=get_triangulation(f.E.value)
        Gridap.writevtk(Ω,name,cellfields=["real(E)"=>real(f.E).value,"imag(E)"=>imag(f.E).value,"normE"=>f.normE.value]);
        return nothing;
    else
        throw(ArgumentError("Not implemented"))
    end
end

function Gridap.:writevtk(name::String,f::VectorFiberEMField)
    if isa(f.E,FEMField)
        Ω=get_triangulation(f.Ex.value)
        Gridap.writevtk(Ω,name,cellfields=["real(Ex)"=>real(f.Ex).value,"imag(Ex)"=>imag(f.Ex).value,"real(Ey)"=>real(f.Ey).value,"imag(Ey)"=>imag(f.Ey).value,"real(Ez)"=>real(f.Ez).value,"imag(Ez)"=>imag(f.Ez).value,"real(Hx)"=>real(f.Hx).value,"imag(Hx)"=>imag(f.Hx).value,"real(Hy)"=>real(f.Hy).value,"imag(Hy)"=>imag(f.Hy).value,"real(Hz)"=>real(f.Hz).value,"imag(Hz)"=>imag(f.Hz).value]);
        return nothing;
    else
        throw(ArgumentError("Not implemented"))
    end
end

Gridap.:writevtk(name::String,m::Union{Mode{ScalarFiberEMField1D},Mode{ScalarFiberEMField2D},Mode{VectorFiberEMField}}) = Gridap.writevtk(name,m.EMField)