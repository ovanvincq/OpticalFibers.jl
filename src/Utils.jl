using LinearAlgebra
using LinearMaps
using SparseArrays
using ArnoldiMethod
using Gridap
using Gridap.Fields
using HCubature
using Unitful

export piecewiseIndex
export meshgrid
export derivative
export ring
export approx_neff_PCF
export approx_nFSM_PCF
export add_cylindrical_PML
export add_rectangular_PML
export add_twist_PML
export nb_args
export compute_kt
export function_integrate

const realLength=Unitful.Quantity{<:Real,Unitful.𝐋}
const realQuantity=Unitful.Quantity{<:Real}
#const inverseLength=Unitful.Quantity{<:Number,Unitful.𝐋^-1}

function Base.conj(g::Gain)
    return conj(g.val)*logunit(g)
end

"""
    piecewiseIndex(pos::T,r::Union{AbstractVector{<:T2},T2},n::AbstractVector) where {T<:Union{Real,realQuantity},T2<:Union{Real,realQuantity}}

Function that returns the refractive index at the position pos in the case of a step-index fiber profile defined by the vectors r and n.
"""
function piecewiseIndex(pos::T,r::Union{AbstractVector{<:T2},T2},n::AbstractVector) where {T<:Union{Real,realQuantity},T2<:Union{Real,realQuantity}}
    if (length(n) != (length(r)+1))
        throw(DimensionMismatch("dim(r) must be equal to dim(n)-1"));
    end
    if !issorted(r)
        throw(DimensionMismatch("r must be sorted"));
    end
    if !isa(r,AbstractVector)
        r=[r]
    end
    return n[searchsortedfirst(r,pos)]
end

"""
    meshgrid(v::AbstractVector)

Function equivalent to the matlab function
"""
meshgrid(v::AbstractVector) = meshgrid(v, v)

"""
    meshgrid(vx::AbstractVector{T}, vy::AbstractVector{T})

Function equivalent to the matlab function
"""
function meshgrid(vx::AbstractVector{T}, vy::AbstractVector{T}) where T
    m, n = length(vy), length(vx)
    vx = reshape(vx, 1, n)
    vy = reshape(vy, m, 1)
    (repeat(vx, m, 1), repeat(vy, 1, n))
end

function isequidistant(v0::AbstractVector{<:Real})
    v=diff(v0)
    return all(isapprox.(v,sum(v)/length(v)))
end

function derivative_iteration(t::Tuple{Vector{<:Union{Real,Quantity{<:Real}}},Vector{<:Union{Number,Quantity}}})
    if (length(t[1])!=length(t[2]))
        throw(ArgumentError("x and y must have the same length"));
    end
    xp=(t[1][1:end-1]+t[1][2:end])/2;
    yp=diff(t[2])./diff(t[1]);
    return xp,yp;
end

"""
    derivative(t::Tuple{Vector{<:Union{Real,Quantity{<:Real}}},Vector{<:Union{Number,Quantity}}},n::Int64=1)

Function that returns the nth order derivative ``\\frac{d^n x}{dy}`` of the tuple t=(x,y) 
"""
function derivative(t::Tuple{Vector{<:Union{Real,Quantity{<:Real}}},Vector{<:Union{Number,Quantity}}},n::Int64=1)
    if (n<0)
        throw(ArgumentError("The order must be positive"));
    end
    temp=t;
    for _=1:n
        temp=derivative_iteration(temp);
    end
    return temp;
end

"""
    ring(N::Integer)

Function that returns a tuple of vectors containing the coordinates (x,y) of the holes in the Nth ring of a PCF (you have to multiply by the pitch to get the real coordinates)
"""
function ring(N::Integer)
    z=zeros(2,6*N+1);
    z[1,2:N+1].=-0.5;
    z[1,N+2:2*N+1].=-1.0;
    z[1,2*N+2:3*N+1].=-0.5;
    z[1,3*N+2:6*N+1].=-z[1,2:3*N+1];
    z[2,2:N+1].=sqrt(3)/2;
    z[2,2*N+2:3*N+1].=-sqrt(3)/2;
    z[2,3*N+2:6*N+1].=-z[2,2:3*N+1];
    VX=0.0;
    VY=0.0;
    x=zeros(6*N);
    y=zeros(6*N);
    for k=1:6*N
        VX=VX+z[1,k];
        VY=VY+z[2,k];
        x[k]=N+VX;
        y[k]=VY;
    end
    x,y
end

#=struct ShiftAndInvert_MUMPS
    m::Mumps
    B
    temp
end

function (M::ShiftAndInvert_MUMPS)(y,x)
    MUMPS.set_job!(M.m,4)
    mul!(M.temp,M.B,x)
    associate_rhs!(M.m,M.temp)
    MUMPS.solve!(M.m)
    MUMPS.get_rhs!(y,M.m)
end

function eigs_MUMPS(A,B;sigma=0,nev::Int64=1,tol::Float64=0.0,restarts::Int64=200,verbose::Bool=false)
    if (!MPI.Initialized())
        MPI.Init();
    end
    icntl = get_icntl(verbose=false)
    m = Mumps{eltype(A)}(mumps_unsymmetric, icntl, default_cntl64)
    associate_matrix!(m,SparseMatrixCSC(A-sigma*B));
    factorize!(m);
    a = ShiftAndInvert_MUMPS(m,B,Vector{eltype(B)}(undef, size(B,1)))
    map_mumps=LinearMap{eltype(A)}(a, size(A,1), ismutating=true)
    if (tol!=0)
        decomp, history  = partialschur(map_mumps, nev=nev, tol=tol, restarts=restarts, which=:LM)
    else
        decomp, history  = partialschur(map_mumps, nev=nev, restarts=restarts, which=:LM)
    end
    if (verbose)
        @show history
    end
    λs_inv, X = partialeigen(decomp);
    MUMPS.finalize!(m);
    λs=(1 ./λs_inv).+sigma
    return λs,X;
end

function eigs_MUMPS(A;sigma=0,nev::Int64=1,tol::Float64=0.0,restarts::Int64=200,verbose::Bool=false)
    eigs_MUMPS(A,eltype(A).(SparseMatrixCSC(I,size(A,1),size(A,2)));sigma=sigma,nev=nev,tol=tol,restarts=restarts,verbose=verbose)
end=#

function eigs_MUMPS(args...;sigma=0,nev::Int64=1,tol::Float64=0.0,restarts::Int64=200,verbose::Bool=false)
     error("MUMPS is not loaded")
end

function eigs_CUDA(args...;sigma=0,nev::Int64=1,tol::Float64=0.0,restarts::Int64=200,verbose::Bool=false,ir_n_steps::Int64=10)
     error("CUDSS is not loaded")
end

struct ShiftAndInvert_LU{TA,TB,TT}
    A_lu::TA
    B::TB
    temp::TT
end

function (M::ShiftAndInvert_LU)(y,x)
    mul!(M.temp, M.B, x)
    ldiv!(y, M.A_lu, M.temp)
end

function eigs_LU(A,B;sigma=0,nev::Int64=1,tol::Float64=0.0,restarts::Int64=200,verbose::Bool=false)
    a = ShiftAndInvert_LU(lu(A-sigma*B),B,Vector{eltype(A)}(undef, size(A,1)))
    map_LU=LinearMap{eltype(A)}(a, size(A,1), ismutating=true)
    if (tol!=0)
        decomp,history  = partialschur(map_LU, nev=nev, tol=tol, restarts=restarts, which=:LM)
    else
        decomp,history  = partialschur(map_LU, nev=nev, restarts=restarts, which=:LM)
    end
    if (verbose)
        @show history
    end
    λs_inv, X = partialeigen(decomp);
    λs=(1 ./λs_inv).+sigma
    return λs,X;
end

function eigs_LU(A;sigma=0,nev::Int64=1,tol::Float64=0.0,restarts::Int64=200,verbose::Bool=false)
    eigs_LU(A,SparseMatrixCSC(I,size(A,1),size(A,2));sigma=sigma,nev=nev,tol=tol,restarts=restarts,verbose=verbose)
end

"""
    add_cylindrical_PML(epsmu::Union{Number,Function},r_pml::Union{Real,realLength},d_pml::Union{Real,realLength},alphaPML::Real)

Function that returns a tensor of permittivity/permeability with a cylindrical PML

- epsmu: permittivity/permeability profile of the fiber. The fiber is assumed to be isotropic when epsmu is a `Number` or a `Function`
- r_pml: distance between the fiber center and the PML beginning
- d_pml: PML thickness
- alphaPML: attenuation factor of the PML
"""
function add_cylindrical_PML(epsmu::Union{Number,Function},r_pml::Union{Real,realLength},d_pml::Union{Real,realLength},alphaPML::Real)
    Gridap.Helpers.@abstractmethod
end


function add_cylindrical_PML(epsmu::Function,r_pml::Union{Real,realLength},d_pml::Union{Real,realLength},alphaPML::Real)
    if (dimension(r_pml)!=dimension(d_pml))
        throw(ArgumentError("r_pml and d_pml must have the same dimension"));
    end
    r(x)=hypot(x[1],x[2]);
    Sr(x)=1.0+im*alphaPML*(r(x)-r_pml)^2/(d_pml)^2;
    Sphi(x)=1.0+im*alphaPML/3.0*(r(x)-r_pml)^3/(d_pml)^2/r(x);
    C(x)=x[1]/r(x);#cos(atan(x[2],x[1]));
    S(x)=x[2]/r(x);#sin(atan(x[2],x[1]));
    if length(epsmu(Point(0,0)*unit(r_pml)))==1
        T1T2(x)=(r(x)<r_pml) ? TensorValue(ComplexF64(1.0),0,0,0,ComplexF64(1.0),0,0,0,ComplexF64(1.0)) : TensorValue(C(x)^2*Sphi(x)/Sr(x)+S(x)^2*Sr(x)/Sphi(x),(Sphi(x)/Sr(x)-Sr(x)/Sphi(x))*C(x)*S(x),0,(Sphi(x)/Sr(x)-Sr(x)/Sphi(x))*C(x)*S(x),S(x)^2*Sphi(x)/Sr(x)+C(x)^2*Sr(x)/Sphi(x),0,0,0,Sr(x)*Sphi(x));
        return x->epsmu(x)*T1T2(x)
    else
        T1(x)=(r(x)<r_pml) ? TensorValue(ComplexF64(1.0),0,0,0,ComplexF64(1.0),0,0,0,ComplexF64(1.0)) : TensorValue(C(x)^2/Sr(x)+S(x)^2/Sphi(x),(1/Sr(x)-1/Sphi(x))*C(x)*S(x),0,(1/Sr(x)-1/Sphi(x))*C(x)*S(x),S(x)^2/Sr(x)+C(x)^2/Sphi(x),0,0,0,1.0);
        T2(x)=(r(x)<r_pml) ? TensorValue(ComplexF64(1.0),0,0,0,ComplexF64(1.0),0,0,0,ComplexF64(1.0)) : TensorValue(C(x)^2*Sphi(x)+S(x)^2*Sr(x),(Sphi(x)-Sr(x))*C(x)*S(x),0,(Sphi(x)-Sr(x))*C(x)*S(x),S(x)^2*Sphi(x)+C(x)^2*Sr(x),0,0,0,Sr(x)*Sphi(x)); #T1/det(T1)
        return x->T1(x)⋅epsmu(x)⋅T2(x)
    end
end

function add_cylindrical_PML(epsmu::Number,r_pml::Union{Real,realLength},d_pml::Union{Real,realLength},alphaPML::Real)
    if (dimension(r_pml)!=dimension(d_pml))
        throw(ArgumentError("r_pml and d_pml must have the same dimension"));
    end
    r(x)=hypot(x[1],x[2]);
    Sr(x)=1.0+im*alphaPML*(r(x)-r_pml)^2/(d_pml)^2;
    Sphi(x)=1.0+im*alphaPML/3.0*(r(x)-r_pml)^3/(d_pml)^2/r(x);
    C(x)=x[1]/r(x);#cos(atan(x[2],x[1]));
    S(x)=x[2]/r(x);#sin(atan(x[2],x[1]));
    if length(epsmu)==1
        T1T2(x)=(r(x)<r_pml) ? TensorValue(ComplexF64(1.0),0,0,0,ComplexF64(1.0),0,0,0,ComplexF64(1.0)) : TensorValue(C(x)^2*Sphi(x)/Sr(x)+S(x)^2*Sr(x)/Sphi(x),(Sphi(x)/Sr(x)-Sr(x)/Sphi(x))*C(x)*S(x),0,(Sphi(x)/Sr(x)-Sr(x)/Sphi(x))*C(x)*S(x),S(x)^2*Sphi(x)/Sr(x)+C(x)^2*Sr(x)/Sphi(x),0,0,0,Sr(x)*Sphi(x));
        return x->epsmu*T1T2(x)
    else
        T1(x)=(r(x)<r_pml) ? TensorValue(ComplexF64(1.0),0,0,0,ComplexF64(1.0),0,0,0,ComplexF64(1.0)) : TensorValue(C(x)^2/Sr(x)+S(x)^2/Sphi(x),(1/Sr(x)-1/Sphi(x))*C(x)*S(x),0,(1/Sr(x)-1/Sphi(x))*C(x)*S(x),S(x)^2/Sr(x)+C(x)^2/Sphi(x),0,0,0,1.0);
        T2(x)=(r(x)<r_pml) ? TensorValue(ComplexF64(1.0),0,0,0,ComplexF64(1.0),0,0,0,ComplexF64(1.0)) : TensorValue(C(x)^2*Sphi(x)+S(x)^2*Sr(x),(Sphi(x)-Sr(x))*C(x)*S(x),0,(Sphi(x)-Sr(x))*C(x)*S(x),S(x)^2*Sphi(x)+C(x)^2*Sr(x),0,0,0,Sr(x)*Sphi(x)); #T1/det(T1)
        return x->T1(x)⋅epsmu⋅T2(x)
    end
end

"""
    add_rectangular_PML(epsmu::Union{Number,Function},x_pml::Union{Real,realLength},dx_pml::Union{Real,realLength},y_pml::Union{Real,realLength},dy_pml::Union{Real,realLength},alphaPML::Real)

Function that returns a tensor of permittivity/permeability with a rectangular PML located at [x_pml,x_pml+dx_pml], [-x_pml-dx_pml,-x_pml], [y_pml,y_pml+dy_pml], [-y_pml-dy_pml,-y_pml]

- epsmu: permittivity/permeability profile of the fiber. The fiber is assumed to be isotropic when epsmu is a `Number` or a `Function`
- x_pml: distance between the fiber center and the PML beginning in the x direction
- dx_pml: PML thickness in the x direction
- y_pml: distance between the fiber center and the PML beginning in the y direction
- dy_pml: PML thickness in the y direction
- alphaPML: attenuation factor of the PML
"""
function add_rectangular_PML(epsmu::Union{Number,Function},x_pml::Union{Real,realLength},dx_pml::Union{Real,realLength},y_pml::Union{Real,realLength},dy_pml::Union{Real,realLength},alphaPML::Real)
    Gridap.Helpers.@abstractmethod
end

function add_rectangular_PML(epsmu::Number,x_pml::Union{Real,realLength},dx_pml::Union{Real,realLength},y_pml::Union{Real,realLength},dy_pml::Union{Real,realLength},alphaPML::Real)
    if !(dimension(x_pml)==dimension(dx_pml)==dimension(y_pml)==dimension(dy_pml))
        throw(ArgumentError("x_pml, dx_pml, y_pml, dy_pml must have the same dimension"));
    end
    Sx=x->1.0+im*alphaPML*(max(abs(x[1])-x_pml,0*unit(x_pml)))^2/(dx_pml)^2;
    Sy=x->1.0+im*alphaPML*(max(abs(x[2])-y_pml,0*unit(x_pml)))^2/(dy_pml)^2;
    if length(epsmu)==1
        S1S2(x)=diagonal_tensor(VectorValue(Sy(x)/Sx(x),Sx(x)/Sy(x),Sx(x)*Sy(x)))
        return x->epsmu*S1S2(x)
    else
        S1(x)=diagonal_tensor(VectorValue(1.0/Sx(x),1.0/Sy(x),1.0));
        S2(x)=diagonal_tensor(VectorValue(Sy(x),Sx(x),Sx(x)*Sy(x)));
        return x->S1(x)⋅epsmu⋅S2(x)
    end
end

function add_rectangular_PML(epsmu::Function,x_pml::Union{Real,realLength},dx_pml::Union{Real,realLength},y_pml::Union{Real,realLength},dy_pml::Union{Real,realLength},alphaPML::Real)
    if !(dimension(x_pml)==dimension(dx_pml)==dimension(y_pml)==dimension(dy_pml))
        throw(ArgumentError("x_pml, dx_pml, y_pml, dy_pml must have the same dimension"));
    end
    Sx=x->1.0+im*alphaPML*(max(abs(x[1])-x_pml,0*unit(x_pml)))^2/(dx_pml)^2;
    Sy=x->1.0+im*alphaPML*(max(abs(x[2])-y_pml,0*unit(x_pml)))^2/(dy_pml)^2;
    if length(epsmu(Point(0,0)*unit(x_pml)))==1
        S1S2(x)=diagonal_tensor(VectorValue(Sy(x)/Sx(x),Sx(x)/Sy(x),Sx(x)*Sy(x)))
        return x->epsmu(x)*S1S2(x)
    else
        S1(x)=diagonal_tensor(VectorValue(1.0/Sx(x),1.0/Sy(x),1.0));
        S2(x)=diagonal_tensor(VectorValue(Sy(x),Sx(x),Sx(x)*Sy(x)));
        return x->S1(x)⋅epsmu(x)⋅S2(x)
    end
end


"""
    add_rectangular_PML(epsmu::Union{Number,Function},xm_pml::Union{Real,realLength},xp_pml::Union{Real,realLength},dx_pml::Union{Real,realLength},ym_pml::Union{Real,realLength},yp_pml::Union{Real,realLength},dy_pml::Union{Real,realLength},alphaPML::Real)

Function that returns a tensor of permittivity/permeability with a rectangular PML located at [xp_pml,xp_pml+dx_pml], [-xm_pml-dx_pml,-xm_pml], [yp_pml,yp_pml+dy_pml], [-ym_pml-dy_pml,-ym_pml]

- epsmu: permittivity/permeability profile of the fiber
- xm_pml: minimum position in the x-direction of the beginning of the PML
- xp_pml: maximum position in the x-direction of the beginning of the PML
- dx_pml: PML thickness in the x direction
- ym_pml: minimum position in the y-direction of the beginning of the PML
- yp_pml: maximum position in the y-direction of the beginning of the PML
- dy_pml: PML thickness in the y direction
- alphaPML: attenuation factor of the PML
"""
function add_rectangular_PML(epsmu::Union{Number,Function},xm_pml::Union{Real,realLength},xp_pml::Union{Real,realLength},dx_pml::Union{Real,realLength},ym_pml::Union{Real,realLength},yp_pml::Union{Real,realLength},dy_pml::Union{Real,realLength},alphaPML::Real)
    Gridap.Helpers.@abstractmethod
end

function add_rectangular_PML(epsmu::Number,xm_pml::Union{Real,realLength},xp_pml::Union{Real,realLength},dx_pml::Union{Real,realLength},ym_pml::Union{Real,realLength},yp_pml::Union{Real,realLength},dy_pml::Union{Real,realLength},alphaPML::Real)
    if !(dimension(xp_pml)==dimension(xm_pml)==dimension(dx_pml)==dimension(yp_pml)==dimension(ym_pml)==dimension(dy_pml))
        throw(ArgumentError("xm_pml, xp_pml, dx_pml, ym_pml, yp_pml, dy_pml must have the same dimension"));
    end
    Sx=x->1.0+im*alphaPML*(max(x[1]-xp_pml,0*unit(xm_pml)))^2/(dx_pml)^2+im*alphaPML*(max(xm_pml-x[1],0*unit(xm_pml)))^2/(dx_pml)^2;
    Sy=x->1.0+im*alphaPML*(max(x[2]-yp_pml,0*unit(xm_pml)))^2/(dy_pml)^2+im*alphaPML*(max(ym_pml-x[2],0*unit(xm_pml)))^2/(dy_pml)^2;
    if length(epsmu)==1
        S1S2(x)=diagonal_tensor(VectorValue(Sy(x)/Sx(x),Sx(x)/Sy(x),Sx(x)*Sy(x)))
        return x->epsmu*S1S2(x)
    else
        S1(x)=diagonal_tensor(VectorValue(1.0/Sx(x),1.0/Sy(x),1.0));
        S2(x)=diagonal_tensor(VectorValue(Sy(x),Sx(x),Sx(x)*Sy(x)));
        return x->S1(x)⋅epsmu⋅S2(x)
    end
end

function add_rectangular_PML(epsmu::Function,xm_pml::Union{Real,realLength},xp_pml::Union{Real,realLength},dx_pml::Union{Real,realLength},ym_pml::Union{Real,realLength},yp_pml::Union{Real,realLength},dy_pml::Union{Real,realLength},alphaPML::Real)
    if !(dimension(xp_pml)==dimension(xm_pml)==dimension(dx_pml)==dimension(yp_pml)==dimension(ym_pml)==dimension(dy_pml))
        throw(ArgumentError("xm_pml, xp_pml, dx_pml, ym_pml, yp_pml, dy_pml must have the same dimension"));
    end
    Sx=x->1.0+im*alphaPML*(max(x[1]-xp_pml,0*unit(xm_pml)))^2/(dx_pml)^2+im*alphaPML*(max(xm_pml-x[1],0*unit(xm_pml)))^2/(dx_pml)^2;
    Sy=x->1.0+im*alphaPML*(max(x[2]-yp_pml,0*unit(xm_pml)))^2/(dy_pml)^2+im*alphaPML*(max(ym_pml-x[2],0*unit(xm_pml)))^2/(dy_pml)^2;
    if length(epsmu(Point(0,0)*unit(xm_pml)))==1
        S1S2(x)=diagonal_tensor(VectorValue(Sy(x)/Sx(x),Sx(x)/Sy(x),Sx(x)*Sy(x)))
        return x->epsmu(x)*S1S2(x)
    else
        S1(x)=diagonal_tensor(VectorValue(1.0/Sx(x),1.0/Sy(x),1.0));
        S2(x)=diagonal_tensor(VectorValue(Sy(x),Sx(x),Sx(x)*Sy(x)));
        return x->S1(x)⋅epsmu(x)⋅S2(x)
    end
end


"""
    add_twist_PML(epsmu::Union{Number,Function},P::Union{Real,realLength},r_pml::Union{Real,realLength},d_pml::Union{Real,realLength},alphaPML::Real)

Function that returns a tensor of permittivity/permeability of a twisted isotropic fiber with a cylindrical PML

- epsmu: function of the tuple (x,y) that describes the permittivity/permeability profile of the fiber
- P: period of the twist
- r_pml: distance between the fiber center and the PML beginning
- d_pml: PML thickness
- alpha: attenuation factor of the PML
"""
function add_twist_PML(epsmu::Union{Number,Function},P::Union{Real,realLength},r_pml::Union{Real,realLength},d_pml::Union{Real,realLength},alphaPML::Real)
    Gridap.Helpers.@abstractmethod
end

function add_twist_PML(epsmu::Function,P::Union{Real,realLength},r_pml::Union{Real,realLength},d_pml::Union{Real,realLength},alphaPML::Real)
    if !(dimension(P)==dimension(r_pml)==dimension(d_pml))
        throw(ArgumentError("P, r_pml, d_pml must have the same dimension"));
    end
    if length(epsmu(Point(0,0)*unit(r_pml)))!=1
        throw(ArgumentError("Twisted PML only works with homogeneous media"));
    end
    alpha=2*pi/P;
    r(x)=hypot(x[1],x[2]);
    Sr(x)=1.0+im*alphaPML*(r(x)-r_pml)^2/(d_pml)^2;
    rp(x)=r(x)+im*alphaPML/3.0*(r(x)-r_pml)^3/(d_pml)^2;
    phi(x)=2*atan(x[2]/(x[1]+r(x)))
    C(x)=cos(phi(x));
    S(x)=sin(phi(x));
    S2(x)=sin(2*phi(x));
    invTxx(x)=(r(x)<r_pml) ? ComplexF64(uconvert(NoUnits,(1.0+alpha^2*x[2]^2)*epsmu(x))) : (C(x)^2/Sr(x)*rp(x)/r(x)+S(x)^2*Sr(x)*r(x)/rp(x)*(1+alpha^2*rp(x)^2))*epsmu(x);
    invTxy(x)=(r(x)<r_pml) ? ComplexF64(uconvert(NoUnits,-alpha^2*x[1]*x[2]*epsmu(x))) : S2(x)*(rp(x)^2-r(x)^2*(1+alpha^2*rp(x)^2)*Sr(x)^2)/(2*r(x)*rp(x)*Sr(x))*epsmu(x)
    invTxz(x)=(r(x)<r_pml) ? ComplexF64(uconvert(NoUnits,-alpha*x[2]*epsmu(x))) : -alpha*rp(x)*Sr(x)*S(x)*epsmu(x)
    invTyy(x)=(r(x)<r_pml) ? ComplexF64(uconvert(NoUnits,(1.0+alpha^2*x[1]^2)*epsmu(x))) : (S(x)^2/Sr(x)*rp(x)/r(x)+C(x)^2*Sr(x)*r(x)/rp(x)*(1+alpha^2*rp(x)^2))*epsmu(x);
    invTyz(x)=(r(x)<r_pml) ? ComplexF64(uconvert(NoUnits,alpha*x[1]*epsmu(x))) : alpha*rp(x)*Sr(x)*C(x)*epsmu(x)
    invTzz(x)=(r(x)<r_pml) ? ComplexF64(uconvert(NoUnits,epsmu(x))) : rp(x)*Sr(x)/r(x)*epsmu(x)
    return x->TensorValue(invTxx(x),invTxy(x),invTxz(x),invTxy(x),invTyy(x),invTyz(x),invTxz(x),invTyz(x),invTzz(x))
end

function add_twist_PML(epsmu::Number,P::Union{Real,realLength},r_pml::Union{Real,realLength},d_pml::Union{Real,realLength},alphaPML::Real)
    if !(dimension(P)==dimension(r_pml)==dimension(d_pml))
        throw(ArgumentError("P, r_pml, d_pml must have the same dimension"));
    end
    if length(epsmu)!=1
        throw(ArgumentError("Twisted PML only works with homogeneous media"));
    end
    alpha=2*pi/P;
    r(x)=hypot(x[1],x[2]);
    Sr(x)=1.0+im*alphaPML*(r(x)-r_pml)^2/(d_pml)^2;
    rp(x)=r(x)+im*alphaPML/3.0*(r(x)-r_pml)^3/(d_pml)^2;
    phi(x)=2*atan(x[2]/(x[1]+r(x)))
    C(x)=cos(phi(x));
    S(x)=sin(phi(x));
    S2(x)=sin(2*phi(x));
    invTxx(x)=(r(x)<r_pml) ? ComplexF64(uconvert(NoUnits,(1.0+alpha^2*x[2]^2)*epsmu)) : (C(x)^2/Sr(x)*rp(x)/r(x)+S(x)^2*Sr(x)*r(x)/rp(x)*(1+alpha^2*rp(x)^2))*epsmu;
    invTxy(x)=(r(x)<r_pml) ? ComplexF64(uconvert(NoUnits,-alpha^2*x[1]*x[2]*epsmu)) : S2(x)*(rp(x)^2-r(x)^2*(1+alpha^2*rp(x)^2)*Sr(x)^2)/(2*r(x)*rp(x)*Sr(x))*epsmu
    invTxz(x)=(r(x)<r_pml) ? ComplexF64(uconvert(NoUnits,-alpha*x[2]*epsmu)) : -alpha*rp(x)*Sr(x)*S(x)*epsmu
    invTyy(x)=(r(x)<r_pml) ? ComplexF64(uconvert(NoUnits,(1.0+alpha^2*x[1]^2)*epsmu)) : (S(x)^2/Sr(x)*rp(x)/r(x)+C(x)^2*Sr(x)*r(x)/rp(x)*(1+alpha^2*rp(x)^2))*epsmu;
    invTyz(x)=(r(x)<r_pml) ? ComplexF64(uconvert(NoUnits,alpha*x[1]*epsmu)) : alpha*rp(x)*Sr(x)*C(x)*epsmu
    invTzz(x)=(r(x)<r_pml) ? ComplexF64(epsmu) : rp(x)*Sr(x)/r(x)*epsmu
    return x->TensorValue(invTxx(x),invTxy(x),invTxz(x),invTxy(x),invTyy(x),invTyz(x),invTxz(x),invTyz(x),invTzz(x))
end

function get_companion(A0,A1,A2)
    n=size(A0,1)
    T=eltype(A0[1]);
    E=spzeros(T,2*n,2*n);
    F=spzeros(T,2*n,2*n);
    E[1:n,1:n]=A2;
    E[n+1:2*n,n+1:2*n]=ones(T).*sparse(I,n,n);
    F[1:n,1:n]=-A1;
    F[1:n,n+1:2*n]=-A0;
    F[n+1:2*n,1:n]=ones(T).*sparse(I,n,n);
    return E,F;
end

"""
    approx_nFSM_PCF(lambda::Real,D::Real,pitch::Real)

Returns an approximate value of nFSM for photonic crystal fiber [Saitoh2005](@cite)  
- lambda: Wavelength
- D: Hole diameter
- pitch: Pitch 

"""
function approx_nFSM_PCF(lambda::Real,D::Real,pitch::Real)
    a=[0.54808 0.71041 0.16904 -1.52736 ; 5.00401 9.73491 1.85765 1.06745 ; -10.43248 47.41496 18.96849 1.93229; 8.22992 -437.50962 -42.4318 3.89]';
    b=[5 1.8 1.7 -0.84 ; 7 7.32 10 1.02 ; 9 22.8 14 13.4]';
    A=zeros(4);
    for i=1:4
        A[i]=a[i,1]+a[i,2]*(D/pitch)^b[i,1]+a[i,3]*(D/pitch)^b[i,2]+a[i,4]*(D/pitch)^b[i,3];
    end
    V=A[1]+A[2]/(1+A[3]*exp(A[4]*lambda/pitch));
    aeff=pitch/sqrt(3)
    return sqrt(1.45^2-(V*lambda/(2*pi*aeff))^2);
end

function approx_nFSM_PCF(lambda::realLength,D::realLength,pitch::realLength)
    return approx_nFSM_PCF(ustrip(u"m",lambda),ustrip(u"m",D),ustrip(u"m",pitch))
end

"""
    approx_neff_PCF(lambda::Real,D::Real,pitch::Real)

Returns an approximate value of the effective index of the fundamental mode of a photonic crystal fiber [Saitoh2005](@cite)  
- lambda: Wavelength
- D: Hole diameter
- pitch: Pitch 

"""
function approx_neff_PCF(lambda::Real,D::Real,pitch::Real)
    nFSM=approx_nFSM_PCF(lambda,D,pitch);
    c=[-0.0973 0.53193 0.24876 5.29801 ; -16.70566 6.70858 2.72423 0.05142 ; 67.13845 52.04855 13.28649 -5.18302 ; -50.25518 -540.66947 -36.80372 2.7641]';
    d=[7 1.49 3.85 -2 ; 9 6.58 10 0.41 ; 10 24.8 15 6]'
    B=zeros(4);
    for i=1:4
        B[i]=c[i,1]+c[i,2]*(D/pitch)^d[i,1]+c[i,3]*(D/pitch)^d[i,2]+c[i,4]*(D/pitch)^d[i,3];
    end
    W=B[1]+B[2]/(1+B[3]*exp(B[4]*lambda/pitch));
    aeff=pitch/sqrt(3)
    return sqrt(nFSM^2+(W*lambda/(2*pi*aeff))^2);
end

function approx_neff_PCF(lambda::realLength,D::realLength,pitch::realLength)
    return approx_neff_PCF(ustrip(u"m",lambda),ustrip(u"m",D),ustrip(u"m",pitch))
end

"""
    nb_args(f::Function)

Returns the possible numbers of arguments of a function

"""
function nb_args(f::Function)
    return unique([methods(f)[i].nargs-1 for i in axes(methods(f),1)])
end

############################# integration #############################
"""
    function_integrate_unitful(f::Function,a::AbstractVector,b::AbstractVector;characteristic_length::AbstractVector=[],kwargs...)

Returns the integral ``\\int_{a_1}^{b_1} \\int_{a_2}^{b_2} ... f(x) dx`` where x is a `Point` of dimension length(a) with unit handling. The integration is performed with a change of variable ```x=tan(t)``` and the `hcubature` function of the `HCubature.jl` package.
- a: Vector of the lower bounds of integration
- b: Vector of the upper bounds of integration
- characteristic_length: Vector of characteristic lengths to make the integration bounds unitless and significant variations around 1. If not provided, the unit of  `a` is used.
- kwargs: keyword arguments of the hcubature function (atol, rtol, initdiv, norm)

"""
function function_integrate_unitful(f::Function,a::AbstractVector,b::AbstractVector;characteristic_length::AbstractVector=[],kwargs...)
    if (length(a)!=length(b))
        throw(ArgumentError("a and b must have the same length"))
    end
    ua=unique(unit.(a))
    ub=unique(unit.(b))
    if ((length(ua)!=1) || (length(ub)!=1) || (dimension.(ua)!=dimension.(ub)))
        throw(ArgumentError("All elements of a and b must have the same unit"))
    end
    ux=ua[1];

    if !isempty(characteristic_length)
        if (length(a)!=length(characteristic_length))
            throw(ArgumentError("a and characteristic_length must have the same length"))
        end
        uc=unique(unit.(characteristic_length))
        if ((length(uc)!=1) || (dimension.(ua)!=dimension.(uc)))
            throw(ArgumentError("All elements of a and characteristic_length must have the same unit"))
        end
    else
        characteristic_length=ones(length(a))*ux;
    end

    x=VectorValue(Tuple(0.0*ux for j=1:length(a)))
    zz=f(x)
    uf=unit(zz[1])
    
    kwargs2=Dict(:norm => x->abs(norm(x)))
    if haskey(kwargs,:atol)
        if dimension(kwargs[:atol])!=dimension(uf*ux^length(a))
            throw(ArgumentError("atol must have the same dimension as the result of integration"))
        end
        kwargs_tmp=Dict(:atol=>ustrip(uf,kwargs[:atol]/prod(characteristic_length)))
        kwargs2=merge(kwargs2,kwargs_tmp)
    end
    if haskey(kwargs,:rtol)
        kwargs_tmp=Dict(:rtol=>kwargs[:rtol])
        kwargs2=merge(kwargs2,kwargs_tmp)
    end
    if haskey(kwargs,:initdiv)
        kwargs_tmp=Dict(:initdiv=>kwargs[:initdiv])
        kwargs2=merge(kwargs2,kwargs_tmp)
    end

    aa=a./characteristic_length;
    bb=b./characteristic_length;
    g(x)=ustrip(uf,f(Point(x.data.*characteristic_length)))
    return uf*prod(characteristic_length)*function_integrate(g,aa,bb;kwargs2...)[1]
end

"""
    function_integrate(f::Function,a::AbstractVector,b::AbstractVector;kwargs...)

Returns the integral ``\\int_{a_1}^{b_1} \\int_{a_2}^{b_2} ... f(x) dx`` where x is a unitless `Point`. The integration is performed with a change of variable ```x=tan(t)``` and the `hcubature` function of the `HCubature.jl` package.
- a: Vector of the lower bounds of integration
- b: Vector of the upper bounds of integration
- kwargs: keyword arguments of the hcubature function (atol, rtol, initdiv, norm)

"""
function function_integrate(f::Function,a::AbstractVector,b::AbstractVector;kwargs...)
    aa=atan.(a)
    bb=atan.(b)
    g2(x)=f(Point((tan.(x.data))))*prod(1.0.+(tan.(x.data)).^2);
    return hcubature(g2,aa,bb;kwargs...)
end
#Peut-être utiliser Integrals.jl à l'avenir mais pas de gain


"""
    compute_kt(KNumber::Int,CellType::Symbol;Irreducible::Bool=true,Pitch::Real=1,MeshType::Symbol=:Internal)

Compute the vectors of a 2D Brillouin zone
- KNumber: Integer related to the number of vectors to compute
- CellType: :Square or :Hexagon
- Irreducible: true-> irreducible Brillouin zone, false->entire Billouin zone
- pitch: Pitch 
- MeshType: Internal (useful to compute the Density Of States) or :Edge (useful to compute the bandgap edge)
"""
function compute_kt(KNumber::Int,CellType::Symbol;Irreducible::Bool=true,Pitch::Union{Real,realLength}=1,MeshType::Symbol=:Internal)
    if (KNumber<1)
        throw(DomainError(KNumber, "KNumber must be at least 1"));
    end
    if (!(CellType in [:Square,:Hexagon]))
        throw(DomainError(CellType, "CellType must be :Square or :Hexagon"));
    end
    if (ustrip(Pitch)<=0)
        throw(DomainError(Pitch, "Pitch must be strictly positive"));
    end
    if (!(MeshType in [:Internal,:Edge]))
        throw(DomainError(MeshType, "MeshType must be :Internal or :Edge"));
    end
    b1=zeros(2)/unit(Pitch)
    b2=zeros(2)/unit(Pitch)
    if (CellType==:Square)
        b1[1]=2*pi/Pitch
        b2[2]=2*pi/Pitch
    else
        b1[1]=2*pi/Pitch
        b1[2]=-2*pi/Pitch/sqrt(3)
        b2[2]=4*pi/Pitch/sqrt(3)
    end
    NumberofKT=1;
    if (KNumber==1)
        return [[0,0]],[1]
    else
        if Irreducible
            if MeshType==:Internal
                NumberofKT=div(KNumber*(KNumber+1),2)
            else
                NumberofKT=3*(KNumber-1)
            end
        else
            if MeshType==:Internal
                if CellType==:Square
                    NumberofKT=KNumber^2
                else
                    NumberofKT=(div(KNumber*(KNumber+1),2)+div((KNumber-2)*(KNumber-1),2))*6-5
                end
            else
                if CellType==:Square
                    NumberofKT=(4*(KNumber-1)-1)*4-3
                else
                    NumberofKT=(4*(KNumber-1)-1)*6-5
                end
            end
        end
        kt=[Vector{typeof(1.0/Pitch)}(undef,2) for _ in 1:NumberofKT]
        weight=Vector{Float64}(undef,NumberofKT)
        counter=1
        if CellType==:Hexagon
            if MeshType==:Internal
                if Irreducible
                    for x=0:(KNumber-1)
                        for y=0:x
                            kt[counter][1]=y*2*pi/3/Pitch/(KNumber-1)
                            kt[counter][2]=x*2*pi/Pitch/sqrt(3)/(KNumber-1)
                            weight[counter]=12
                            if (x+y==0)
                                weight[counter]=1
                            else
                                if (y==x)
                                    weight[counter]=6
                                end
                                if (y==0)
                                    weight[counter]=6
                                end
                                if (x==KNumber-1)
                                    weight[counter]=6
                                    if (y==0)
                                        weight[counter]=3
                                    end
                                    if (y==x)
                                        weight[counter]=2
                                    end
                                end
                            end
                            counter=counter+1
                        end
                    end
                else
                    for x=0:(KNumber-1)
                        for y=0:x
                            kt[counter][1]=y*2*pi/3/Pitch/(KNumber-1)
                            kt[counter][2]=x*2*pi/sqrt(3)/Pitch/(KNumber-1)
                            weight[counter]=6
                            if (x==KNumber-1)
                                weight[counter]=3
                                if (y==KNumber-1)
                                    weight[counter]=2
                                end
                            end
                            counter=counter+1
                        end
                    end
                    for x=1:(KNumber-1)
                        for y=1:(x-1)
                            kt[counter][1]=-y*2*pi/3/Pitch/(KNumber-1)
                            kt[counter][2]=x*2*pi/Pitch/sqrt(3)/(KNumber-1)
                            weight[counter]=6
                            if (x==(KNumber-1))
                                weight[counter]=3
                            end
                            counter=counter+1
                        end
                    end
                    aa=counter
                    for t=1:5
                        for i=2:(aa-1)
                            kt[counter][1]=kt[i][1]*cospi(t/3)-kt[i][2]*sinpi(t/3)
                            kt[counter][2]=kt[i][1]*sinpi(t/3)+kt[i][2]*cospi(t/3)
                            weight[counter]=weight[i]
                            counter=counter+1
                        end
                    end
                end
            else
                if Irreducible
                    for i=0:(KNumber-1)
                        kt[counter][1]=0/Pitch
                        kt[counter][2]=i*2*pi/sqrt(3)/(KNumber-1)/Pitch
                        weight[counter]=6
                        if (i==0)
                            weight[counter]=1
                        end
                        if (i==KNumber-1)
                            weight[counter]=3
                        end
                        counter=counter+1
                    end
                    for i=1:(KNumber-2)
                        kt[counter][1]=i*2*pi/3/(KNumber-1)/Pitch
                        kt[counter][2]=2*pi/sqrt(3)/Pitch
                        weight[counter]=6
                        counter=counter+1
                    end
                    for i=(KNumber-1):-1:1
                        kt[counter][1]=i*2*pi/3/(KNumber-1)/Pitch
                        kt[counter][2]=i*2*pi/sqrt(3)/(KNumber-1)/Pitch
                        weight[counter]=6
                        if (i==(KNumber-1))
                            weight[counter]=2
                        end
                        counter=counter+1
                    end
                else
                    for i=0:KNumber-1
                        kt[counter][1]=0/Pitch
                        kt[counter][2]=i*2*pi/sqrt(3)/(KNumber-1)/Pitch
                        weight[counter]=6
                        if (i==KNumber-1)
                            weight[counter]=3
                        end
                        counter=counter+1
                    end
                    for i=1:KNumber-2
                        kt[counter][1]=i*2*pi/3/(KNumber-1)/Pitch
                        kt[counter][2]=2*pi/sqrt(3)/Pitch
                        weight[counter]=3
                        counter=counter+1
                    end
                    for i=1:KNumber-2
                        kt[counter][1]=-i*2*pi/3/(KNumber-1)/Pitch
                        kt[counter][2]=2*pi/sqrt(3)/Pitch
                        weight[counter]=3
                        counter=counter+1
                    end
                    for i=(KNumber-1):-1:1
                        kt[counter][1]=i*2*pi/3/(KNumber-1)/Pitch
                        kt[counter][2]=i*2*pi/sqrt(3)/(KNumber-1)/Pitch
                        weight[counter]=6
                        if (i==KNumber-1)
                            weight[counter]=2
                        end
                        counter=counter+1
                    end
                    for t=1:5
                        for i=2:(4*(KNumber-1)-1)
                            kt[counter][1]=kt[i][1]*cospi(t/3)-kt[i][2]*sinpi(t/3)
                            kt[counter][2]=kt[i][1]*sinpi(t/3)+kt[i][2]*cospi(t/3)
                            weight[counter]=weight[i]
                            counter=counter+1
                        end
                    end
                end
            end
        else
            if MeshType==:Internal
                if Irreducible
                    for x=0:(KNumber-1)
                        for y=0:x
                            kt[counter][1]=x*pi/Pitch/(KNumber-1)
                            kt[counter][2]=y*pi/Pitch/(KNumber-1)
                            weight[counter]=8
                            if (x+y==0)
                                weight[counter]=1
                            else
                                if (y==0)
                                    weight[counter]=4
                                end
                                if (x==y)
                                    weight[counter]=4
                                end
                                if (x==(KNumber-1))
                                    weight[counter]=4
                                    if (y==0)
                                        weight[counter]=2
                                    end
                                    if (y==(KNumber-1))
                                        weight[counter]=1
                                    end
                                end
                            end
                            counter=counter+1
                        end
                    end
                else
                    for i=0:(KNumber-1)
                        for j=0:(KNumber-1)
                            kt[counter][1]=-pi/Pitch+i*2*pi/Pitch/(KNumber-1)
                            kt[counter][2]=-pi/Pitch+j*2*pi/Pitch/(KNumber-1)
                            weight[counter]=4
                            if (i==0)
                                weight[counter]=2
                            end
                            if (i==(KNumber-1))
                                weight[counter]=2
                            end
                            if ((j==0) || (j==(KNumber-1)))
                                weight[counter]=2
                                if ((i==0) || (i==(KNumber-1)))
                                    weight[counter]=1
                                end
                            end
                            counter=counter+1
                        end
                    end
                end
            else
                if Irreducible
                    for i=0:(KNumber-1)
                        kt[counter][1]=i*pi/(KNumber-1)/Pitch
                        kt[counter][2]=0/Pitch
                        weight[counter]=4
                        if (i==0)
                            weight[counter]=1
                        end
                        if (i==(KNumber-1))
                            weight[counter]=2
                        end
                        counter=counter+1
                    end
                    for i=1:(KNumber-2)
                        kt[counter][1]=i*pi/(KNumber-1)/Pitch
                        kt[counter][2]=i*pi/(KNumber-1)/Pitch
                        weight[counter]=4
                        counter=counter+1
                    end
                    for i=(KNumber-1):-1:1
                        kt[counter][1]=pi/Pitch
                        kt[counter][2]=i*pi/(KNumber-1)/Pitch
                        if (i==(KNumber-1))
                            weight[counter]=1
                        else
                            weight[counter]=4
                        end
                        counter=counter+1
                    end
                else
                    for i=0:(KNumber-1)
                        kt[counter][1]=i*pi/(KNumber-1)/Pitch
                        kt[counter][2]=0/Pitch
                        weight[counter]=4
                        if (i==(KNumber-1))
                            weight[counter]=2
                        end
                        counter=counter+1
                    end
                    for i=1:(KNumber-2)
                        kt[counter][1]=i*pi/(KNumber-1)/Pitch
                        kt[counter][2]=i*pi/(KNumber-1)/Pitch
                        weight[counter]=4
                        counter=counter+1
                    end
                    for i=(KNumber-1):-1:1
                        kt[counter][1]=pi/Pitch
                        kt[counter][2]=i*pi/(KNumber-1)/Pitch
                        if (i==KNumber-1)
                            weight[counter]=1
                        else
                            weight[counter]=2
                        end
                        counter=counter+1
                    end
                    for i=1:(KNumber-2)
                        kt[counter][1]=i*pi/(KNumber-1)/Pitch
                        kt[counter][2]=pi/Pitch
                        weight[counter]=2
                        counter=counter+1
                    end
                    for t=1:3
                        for i=2:(4*(KNumber-1)-1)
                            kt[counter][1]=kt[i][1]*cospi(t/2)-kt[i][2]*sinpi(t/2)
                            kt[counter][2]=kt[i][1]*sinpi(t/2)+kt[i][2]*cospi(t/2)
                            weight[counter]=weight[i]
                            counter=counter+1
                        end
                    end
                end
            end
        end
        weight=weight/sum(weight)
        return kt,weight
    end
end