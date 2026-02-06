module PhysicalData

using StaticArrays
using Unitful

export n_Ge_Doped_Silica
export n_F_Doped_Silica
export n_Silicon
export n_Germanium
export Sa_Ytterbium
export Se_Ytterbium
export c_Unitful, mu0_Unitful, eps0_Unitful, h_Unitful, Z0_Unitful		

"""
    c=299792458.0

Velocity of light in vacuum (in m/s)
"""
const c=299792458.0;

"""
    c_Unitful=299792458.0*u"m/s"

Velocity of light in vacuum
"""
const c_Unitful=c*u"m/s"
"""
    mu0=1.25663706212e-6

Vacuum permeability (in H/m)
"""
const mu0=4E-7*pi;
"""
    mu0_Unitful=1.25663706212e-6*u"H/m"

Vacuum permeability
"""
const mu0_Unitful=mu0*u"H/m"
"""
    eps0=1/(mu0*c^2)

Vacuum permittivity (in F/m)
"""
const eps0=1.0/(mu0*c^2);
"""
    eps0_Unitful=1/(mu0_Unitful*c_Unitful^2)

Vacuum permittivity
"""
const eps0_Unitful=eps0*u"F/m"
"""
    h=6.62607015e-34

Planck constant (in J.s)
"""
const h=6.626_070_15e-34;
"""
    h_Unitful=6.62607015e-34*u"J*s"

Planck constant
"""
const h_Unitful=h*u"J*s"
"""
    Z0=sqrt(mu0/eps0)

Vaccum impedance (in Ω)
"""
const Z0=mu0*c;
"""
    Z0_Unitful=sqrt(mu0_Unitful/eps0_Unitful)

Vaccum impedance
"""
const Z0_Unitful=Z0*u"Ω"

#Refractive index functions

"""
    n_Ge_Doped_Silica(lambda::Union{Real,Quantity{<:Real,Unitful.𝐋}},xGe::Real;author::Symbol=:Fleming)
 
Returns the refractive index of germanium-doped silica with a germanium fraction `xGe` (0 ≤ xGe ≤ 1) at the wavelength `lambda` (in meters if `lambda` is a `Real`).
Two Sellmeier formulas can be used : `author=:Fleming` [Fleming1984](@cite) and `author=:Sunak` [Sunak1989](@cite).

Domain of validity:
- 360 nm ≤`lambda`≤ 4300 nm for Fleming's formula
- 600 nm ≤`lambda`≤ 1800 nm and `xGe` ≤ 0.2 for Sunak's formula
"""
function n_Ge_Doped_Silica(lambda::Real,xGe::Real;author::Symbol=:Fleming)
    if (author==:Fleming)
        SiO = [0.69616630,0.068404300,0.40794260,0.11624140,0.89747940,9.8961610];
        GeO = [0.80686642,0.068972606,0.71815848,0.15396605,0.85416831,11.841931];
        result=0.0;
        for i=1:3
            @inbounds @fastmath result = result + (SiO[2*i-1] + xGe*(GeO[2*i-1]-SiO[2*i-1]))*(lambda*1E6)^2/ ((lambda*1e6)^2 - (SiO[2*i] + xGe*(GeO[2*i]-SiO[2*i]))^2);
        end
        return sqrt(1+result);
    elseif (author==:Sunak)
        ASi=[0.2045154578,0.06451676258,0.1311583151];
        zSi=[0.06130807320,0.1108859848,8.964441861];
        Bi=[-0.1011783769,0.1778934999,-0.1064179581];
        result=0.0;
        for i=1:3
            @inbounds @fastmath result=result+(ASi[i]+Bi[i]*xGe)*(lambda*1E6)^2/((lambda*1E6)^2-zSi[i]^2);
        end
        return sqrt((2*result+1)/(1-result));
    else
        return 0.0;
    end
end

function n_Ge_Doped_Silica(lambda::Unitful.Quantity{<:Real,Unitful.𝐋},xGe::Real;author::Symbol=:Fleming)
    return n_Ge_Doped_Silica(ustrip(u"m",lambda),xGe;author=author)
end

"""
    n_F_Doped_Silica(lambda::Union{Real,Quantity{<:Real,Unitful.𝐋}},xF::Real)

Returns the refractive index of fluorine-doped silica with a fluorine fraction `xF` (0 ≤ xGe ≤ 1) at the wavelength `lambda` (in meters if `lambda` is a `Real`) [Sunak1989](@cite).

Domain of validity: 600 nm ≤`lambda`≤ 1800 nm and `xF` ≤ 0.2
"""
function n_F_Doped_Silica(lambda::Real,xF::Real)
    ASi=[0.2045154578,0.06451676258,0.1311583151];
    zSi=[0.06130807320,0.1108859848,8.964441861];
    Bi=[-0.05413938039,-0.1788588824,-0.07445931332];
    result=0.0;
    for i=1:3
        @inbounds @fastmath result=result+(ASi[i]+Bi[i]*xF)*(lambda*1e6)^2/((lambda*1e6)^2-zSi[i]^2);
    end
    return sqrt((2*result+1)/(1-result));
end

function n_F_Doped_Silica(lambda::Unitful.Quantity{<:Real,Unitful.𝐋},xF::Real)
    return n_F_Doped_Silica(ustrip(u"m",lambda),xF)
end

"""
    n_Silicon(lambda::Union{Real,Quantity{<:Real,Unitful.𝐋}},T::Union{Real,Quantity{<:Real,Unitful.𝚯}})

Returns the refractive index of Silicon at wavelength `lambda` (in meters if `lambda` is a `Real`) and at temperature `T` (in Kelvin if `T` is a `Real`) [Frey2006](@cite)  

Domain of validity: 1100 nm ≤ `lambda` ≤ 5600 nm and 20 K ≤ `T` ≤ 300 K
"""
function n_Silicon(lambda::Real,T::Real)
    S1_coef=[10.4907,-2.08020E-4,4.21694E-6,-5.82298E-9,3.44688E-12];
    S2_coef=[-1346.61,29.1664,-0.278724,1.05939E-3,-1.35089E-6];
    S3_coef=[4.42827E7,-1.76213E6,-7.61575E4,678.414,103.243];
    lambda1_coef=[0.299713,-1.14234E-5,1.67134E-7,-2.51049E-10,2.32484E-14];
    lambda2_coef=[-3.51710E3,42.3892,-0.357957,1.17504E-3,-1.13212E-6];
    lambda3_coef=[1.71400E6,-1.44984E5,-6.90744E3,-39.3699,23.5770];
    TT=[1,T,T^2,T^3,T^4];
    S1=sum(S1_coef.*TT);
    S2=sum(S2_coef.*TT);
    S3=sum(S3_coef.*TT);
    lambda1=sum(lambda1_coef.*TT);
    lambda2=sum(lambda2_coef.*TT);
    lambda3=sum(lambda3_coef.*TT);
    return sqrt(1+S1*(lambda*1e6)^2/((lambda*1e6)^2-lambda1^2)+S2*(lambda*1e6)^2/((lambda*1e6)^2-lambda2^2)+S3*(lambda*1e6)^2/((lambda*1e6)^2-lambda3^2));
end

function n_Silicon(lambda::Unitful.Quantity{<:Real,Unitful.𝐋},T::Unitful.Quantity{<:Real,Unitful.𝚯})
    return n_Silicon(ustrip(u"m",lambda),ustrip(u"K",T))
end

"""
    n_Germanium(lambda::Real,T::Real)

Returns the refractive index of Germanium at wavelength `lambda` (in meters if `lambda` is a `Real`) and at temperature `T` (in Kelvin if `T` is a `Real`) [Frey2006](@cite)  

Domain of validity: 1900 nm ≤ `lambda` ≤ 5500 nm and 20 K ≤ `T` ≤ 300 K

"""
function n_Germanium(lambda::Real,T::Real)
    S1_coef=[13.9723,2.52809E-3,-5.02195E-6,2.22604E-8,-4.86238E-12];
    S2_coef=[0.452096,-3.09197E-3,2.16895E-5,-6.02290E-8,4.12038E-11];
    S3_coef=[751.447,-14.2843,-0.238093,2.96047E-3,-7.73454E-6];
    lambda1_coef=[0.386367,2.01871E-4,-5.93448E-7,-2.27923E-10,5.37423E-12];
    lambda2_coef=[1.08843,1.16510E-3,-4.97284E-6,1.12357E-8,9.40201E-12];
    lambda3_coef=[-2893.19,-0.967948,-0.527016,6.49364E-3,-1.95162E-5];
    TT=[1,T,T^2,T^3,T^4];
    S1=sum(S1_coef.*TT);
    S2=sum(S2_coef.*TT);
    S3=sum(S3_coef.*TT);
    lambda1=sum(lambda1_coef.*TT);
    lambda2=sum(lambda2_coef.*TT);
    lambda3=sum(lambda3_coef.*TT);
    return sqrt(1+S1*(1lambda*1e6)^2/((lambda*1e6)^2-lambda1^2)+S2*(lambda*1e6)^2/((lambda*1e6)^2-lambda2^2)+S3*(lambda*1e6)^2/((lambda*1e6)^2-lambda3^2));
end

function n_Germanium(lambda::Unitful.Quantity{<:Real,Unitful.𝐋},T::Unitful.Quantity{<:Real,Unitful.𝚯})
    return n_Germanium(ustrip(u"m",lambda),ustrip(u"K",T))
end


"""
    Sa_Ytterbium(lambda::Union{Real,Quantity{<:Real,Unitful.𝐋}})

Returns the absorption cross section (in m²) of ytterbium at wavelength `lambda` (in meters if `lambda` is a `Real`) [Marciante2006](@cite)  
"""
function Sa_Ytterbium(lambda::Real)
    A=[180.0,360.0,510.0,160.0,2325.0]*1E-27;
    l=[950.0,895.0,918.0,971.0,975.0]*1E-9;
    w=[70.0,24.0,22.0,12.0,4.0]*1E-9;
    Sa=0.0;
    for i=1:5
        @inbounds @fastmath Sa=Sa+A[i]*exp(-(Float64(lambda)-l[i])^2/w[i]^2);
    end
    return Sa;
end 

function Sa_Ytterbium(lambda::Unitful.Quantity{<:Real,Unitful.𝐋})
    return Sa_Ytterbium(ustrip(u"m",lambda))*u"m^2"
end

"""
    Se_Ytterbium(lambda::Real)

Returns the emission cross section (in m²) of ytterbium at wavelength `lambda` (in meters if `lambda` is a `Real`) [Marciante2006](@cite)  
"""


function Se_Ytterbium(lambda::Real)
    A=[2325.0,160.0,340.0,175.0,150.0]*1E-27;
    l=[975.0,978.0,1025.0,1050.0,1030.0]*1E-9;
    w=[4.0,12.0,20.0,60.0,90.0]*1E-9;
    Se=0.0;
    for i=1:5
        @inbounds @fastmath Se=Se+A[i]*exp(-(Float64(lambda)-l[i])^2/w[i]^2);
    end 
    return Se;
end 

function Se_Ytterbium(lambda::Unitful.Quantity{<:Real,Unitful.𝐋})
    return Se_Ytterbium(ustrip(u"m",lambda))*u"m^2"
end

end


# module PhysicalData