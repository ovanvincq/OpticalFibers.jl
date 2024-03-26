using LinearAlgebra
using LinearMaps
using MUMPS
using MPI
using SparseArrays
using ArnoldiMethod
using Gridap
using Gridap.Fields
using HCubature

export piecewiseIndex
export meshgrid
export derivative
export ring
export tensor3
export isaNumber
export isaFunction
export inverse
export tensorComponent
export approx_neff_PCF
export approx_nFSM_PCF
export add_cylindrical_PML
export add_rectangular_PML
export add_twist_PML
export nb_args
export integrate1D
export integrate2D


"""
    piecewiseIndex(pos::Real,r::Union{AbstractVector{<:Real},Real},n::AbstractVector{<:Real})

Function that returns the refractive index at the position pos in the case of a step-index fiber profile defined by the vectors r and n.
"""
function piecewiseIndex(pos::Real,r::Union{AbstractVector{<:Real},Real},n::AbstractVector{<:Real})
    if (length(n) != (length(r)+1))
        throw(DimensionMismatch("dim(r) must be equal to dim(n)-1"));
    end
    dn=diff(n);
    result=n[1];
    for j=1:length(dn)
        result=result+dn[j]*(pos>r[j]);
    end
    return result;
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

"""
    derivative(t::Tuple{Vector{<:Real},Vector{<:Number}})

Function that returns the derivative dy/dx of the tuple t=(x,y)
"""
function derivative(t::Tuple{Vector{<:Real},Vector{<:Number}})
    if (length(t[1])!=length(t[2]))
        throw(ArgumentError("x and y must have the same length"));
    end
    xp=(t[1][1:end-1]+t[1][2:end])/2;
    yp=diff(t[2])./diff(t[1]);
    return xp,yp;
end

"""
    derivative(t::Tuple{Vector{<:Real},Vector{<:Number}},order::Int64)

Function that returns the nth order derivative dy/dx of the tuple t=(x,y) 
"""
function derivative(t::Tuple{Vector{<:Real},Vector{<:Number}},order::Int64)
    if (order<0)
        throw(ArgumentError("The order must be positive"));
    end
    temp=t;
    for j=1:order
        temp=derivative(temp);
    end
    return temp;
end

"""
    ring(N::Integer)

Function that returns a tuple of vectors containing the coordinates (x,y) of the holes in the Nth ring of a PCF
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

struct ShiftAndInvert_MUMPS
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
    tensorComponent

Structure that describes a the component of the tensor of permittivity or permeability.

- f::`Union{Function,Number}`
"""
struct tensorComponent
    f::Union{Function,Number}
end

function tensorComponent(x::tensorComponent)
    return tensorComponent(x.f)
end

isaNumber(tc::tensorComponent) = isa(tc.f,Number)
isaFunction(tc::tensorComponent) = isa(tc.f,Function)

Base.show(io::IO, tc::tensorComponent) = (isa(tc.f,Number)) ? print(io,tc.f) : print(io,"fun");

function (tc::tensorComponent)(p::Point)
    if isa(tc.f,Number)
        return tc.f
    else
        return tc.f(p)
    end
end

Base.:*(k::Number, tc::tensorComponent) = (isaNumber(tc)) ? tensorComponent(k*tc.f) : ((k==0) ? tensorComponent(0) : ((k==1) ? tensorComponent(tc.f) : tensorComponent(x->k*tc.f(x))))
Base.:*(tc::tensorComponent,k::Number) = k*tc
Base.:/(k::Number, tc::tensorComponent) = (isaNumber(tc)) ? tensorComponent(k/tc.f) : ((k==0) ? tensorComponent(0) : tensorComponent(x->k/tc.f(x)))
Base.:/(tc::tensorComponent,k::Number) = (isaNumber(tc)) ? tensorComponent(tc.f/k) : tensorComponent(x->tc.f(x)/k)
Base.:+(k::Number, tc::tensorComponent) = (isaNumber(tc)) ? tensorComponent(k+tc.f) : ((k==0) ? tensorComponent(tc.f) : tensorComponent(x->k+tc.f(x)))
Base.:+(tc::tensorComponent,k::Number) = k+tc
Base.:-(k::Number, tc::tensorComponent) = (isaNumber(tc)) ? tensorComponent(k-tc.f) : tensorComponent(x->k-tc.f(x))
Base.:-(tc::tensorComponent,k::Number) = (isaNumber(tc)) ? tensorComponent(tc.f-k) : ((k==0) ? tensorComponent(tc.f) : tensorComponent(x->tc.f(x)-k))

Base.:*(fun::Function, tc::tensorComponent) = (isaNumber(tc)) ? ((tc.f==0) ? tensorComponent(0) : ((tc.f==1) ? tensorComponent(fun) : tensorComponent(x->fun(x)*tc.f))) : tensorComponent(x->fun(x)*tc.f(x))
Base.:*(tc::tensorComponent,fun::Function) = fun*tc
Base.:/(fun::Function, tc::tensorComponent) = (isaNumber(tc)) ? tensorComponent(x->fun(x)/tc.f) : tensorComponent(x->fun(x)/tc.f(x))
Base.:/(tc::tensorComponent,fun::Function) = (isaNumber(tc)) ? tensorComponent(x->tc.f/fun(x)) : tensorComponent(x->tc.f(x)/fun(x))
Base.:+(fun::Function, tc::tensorComponent) = (isaNumber(tc)) ? ((tc.f==0) ? tensorComponent(fun) : tensorComponent(x->fun(x)+tc.f)) : tensorComponent(x->fun(x)+tc.f(x))
Base.:+(tc::tensorComponent,fun::Function) = fun+tc
Base.:-(fun::Function, tc::tensorComponent) = (isaNumber(tc)) ? tensorComponent(x->fun(x)-tc.f) : tensorComponent(x->fun(x)-tc.f(x))
Base.:-(tc::tensorComponent,fun::Function) = (isaNumber(tc)) ? tensorComponent(x->tc.f-fun(x)) : tensorComponent(x->tc.f(x)-fun(x))

Base.:*(tc1::tensorComponent, tc2::tensorComponent) = (isaNumber(tc1)) ? tc1.f*tc2 : tc1.f*tc2
Base.:/(tc1::tensorComponent, tc2::tensorComponent) = (isaNumber(tc1)) ? tc1.f/tc2 : tc1/tc2.f
Base.:+(tc1::tensorComponent, tc2::tensorComponent) = (isaNumber(tc1)) ? tc1.f+tc2 : tc1+tc2.f
Base.:-(tc1::tensorComponent, tc2::tensorComponent) = (isaNumber(tc1)) ? tc1.f-tc2 : tc1-tc2.f
Base.:+(tc::tensorComponent) = tc
Base.:-(tc::tensorComponent) = (isa(tc.f,Number)) ? tensorComponent(-tc.f) : tensorComponent(x->-tc.f(x))

"""
    tensor

Structure that describes a tensor of permittivity or permeability. Each of the 9 components of the a `tensor3` is a function of a tuple (x,y).

- xx :: `tensorComponent`
- yx :: `tensorComponent`
- zx :: `tensorComponent`
- xy :: `tensorComponent`
- yy :: `tensorComponent`
- zy :: `tensorComponent`
- xz :: `tensorComponent`
- yz :: `tensorComponent`
- zz :: `tensorComponent`
"""
struct tensor3
    xx::tensorComponent
    yx::tensorComponent
    zx::tensorComponent
    xy::tensorComponent
    yy::tensorComponent
    zy::tensorComponent
    xz::tensorComponent
    yz::tensorComponent
    zz::tensorComponent
    tensor3(xx::Union{Function,Number,tensorComponent},yx::Union{Function,Number,tensorComponent},zx::Union{Function,Number,tensorComponent},xy::Union{Function,Number,tensorComponent},yy::Union{Function,Number,tensorComponent},zy::Union{Function,Number,tensorComponent},xz::Union{Function,Number,tensorComponent},yz::Union{Function,Number,tensorComponent},zz::Union{Function,Number,tensorComponent}) = new(tensorComponent(xx),tensorComponent(yx),tensorComponent(zx),tensorComponent(xy),tensorComponent(yy),tensorComponent(zy),tensorComponent(xz),tensorComponent(yz),tensorComponent(zz))
    tensor3(xx::Union{Function,Number,tensorComponent},yy::Union{Function,Number,tensorComponent},zz::Union{Function,Number,tensorComponent}) = new(tensorComponent(xx),tensorComponent(0),tensorComponent(0),tensorComponent(0),tensorComponent(yy),tensorComponent(0),tensorComponent(0),tensorComponent(0),tensorComponent(zz))
    tensor3(xx::Union{Function,Number,tensorComponent}) = new(tensorComponent(xx),tensorComponent(0),tensorComponent(0),tensorComponent(0),tensorComponent(xx),tensorComponent(0),tensorComponent(0),tensorComponent(0),tensorComponent(xx))
end

function Base.show(io::IO, t::tensor3)
    println(io,t.xx)
    println(io,t.yx)
    println(io,t.zx)
    println(io,t.xy)
    println(io,t.yy)
    println(io,t.zy)
    println(io,t.xz)
    println(io,t.yz)
    print(io,t.zz)
end

"""
    det(t::tensor3)

Returns the determinant of the tensor t
"""
function LinearAlgebra.:det(t::tensor3)
    return t.xx*t.yy*t.zz+t.xy*t.yz*t.zx+t.xz*t.yx*t.zy-t.xx*t.yz*t.zy-t.xy*t.yx*t.zz-t.xz*t.yy*t.zx
end

"""
    inverse(t::tensor3)

Returns the inverse of the tensor t
"""
function inverse(t::tensor3)
    d=det(t);
    xx=(t.yy*t.zz-t.yz*t.zy)/d
    yx=(t.yz*t.zx-t.yx*t.zz)/d
    zx=(t.yx*t.zy-t.yy*t.zx)/d
    xy=(t.xz*t.zy-t.xy*t.zz)/d
    yy=(t.xx*t.zz-t.xz*t.zx)/d
    zy=(t.xy*t.zx-t.xx*t.zy)/d
    xz=(t.xy*t.yz-t.xz*t.yy)/d
    yz=(t.xz*t.yx-t.xx*t.yz)/d
    zz=(t.xx*t.yy-t.xy*t.yx)/d
    return tensor3(xx,yx,zx,xy,yy,zy,xz,yz,zz);
end

Base.:*(k::Union{Number,Function,tensorComponent}, t::tensor3) = tensor3(k*t.xx,k*t.yx,k*t.zx,k*t.xy,k*t.yy,k*t.zy,k*t.xz,k*t.yz,k*t.zz);
Base.:*(t::tensor3, k::Union{Number,Function,tensorComponent}) = k*t;
Base.:+(t1::tensor3, t2::tensor3) = tensor3(t1.xx+t2.xx,t1.yx+t2.yx,t1.zx+t2.zx,t1.xy+t2.xy,t1.yy+t2.yy,t1.zy+t2.zy,t1.xz+t2.xz,t1.yz+t2.yz,t1.zz+t2.zz);
Base.:+(t::tensor3) = t;
Base.:-(t1::tensor3, t2::tensor3) = tensor3(t1.xx-t2.xx,t1.yx-t2.yx,t1.zx-t2.zx,t1.xy-t2.xy,t1.yy-t2.yy,t1.zy-t2.zy,t1.xz-t2.xz,t1.yz-t2.yz,t1.zz-t2.zz);
Base.:-(t::tensor3) = tensor3(-t.xx,-t.yx,-t.zx,-t.xy,-t.yy,-t.zy,-t.xz,-t.yz,-t.zz);
Base.:/(t::tensor3, k::Union{Number,Function,tensorComponent}) = tensor3(t.xx/k,t.yx/k,t.zx/k,t.xy/k,t.yy/k,t.zy/k,t.xz/k,t.yz/k,t.zz/k);
Base.:/(k::Union{Number,Function,tensorComponent}, t::tensor3) = k*inverse(t);
Base.:*(t1::tensor3, t2::tensor3) = tensor3(t1.xx*t2.xx+t1.xy*t2.yx+t1.xz*t2.zx
,t1.yx*t2.xx+t1.yy*t2.yx+t1.yz*t2.zx
,t1.zx*t2.xx+t1.zy*t2.yx+t1.zz*t2.zx
,t1.xx*t2.xy+t1.xy*t2.yy+t1.xz*t2.zy
,t1.yx*t2.xy+t1.yy*t2.yy+t1.yz*t2.zy
,t1.zx*t2.xy+t1.zy*t2.yy+t1.zz*t2.zy
,t1.xx*t2.xz+t1.xy*t2.yz+t1.xz*t2.zz
,t1.yx*t2.xz+t1.yy*t2.yz+t1.yz*t2.zz
,t1.zx*t2.xz+t1.zy*t2.yz+t1.zz*t2.zz);
Base.:/(t1::tensor3, t2::tensor3) = t1*inverse(t2);


################################

"""
    add_cylindrical_PML(epsmu::Union{Number,Function,tensor3},r_pml::Real,d_pml::Real,alpha::Real)

Function that returns a tensor of permittivity/permeability with a cylindrical PML

- epsmu: permittivity/permeability profile of the fiber. The fiber is assumed to be isotropic when epsmu is a `Number` or a `Function`
- r_pml: distance between the fiber center and the PML beginning
- d_pml: PML thickness
- alphaPML: attenuation factor of the PML
"""
function add_cylindrical_PML(epsmu::Union{Number,Function,tensor3},r_pml::Real,d_pml::Real,alphaPML::Real)
    Gridap.Helpers.@abstractmethod
end

function add_cylindrical_PML(epsmu::Number,r_pml::Real,d_pml::Real,alphaPML::Real)
    function pml_xx(x)
        r=hypot(x[1],x[2]);
        if r<=r_pml
            return ComplexF64(epsmu);
        else
            phi=atan(x[2],x[1]);
            rt=r+im*alphaPML/3.0*(r-r_pml)^3/(d_pml)^2;
            sr=1.0+im*alphaPML*(r-r_pml)^2/(d_pml)^2;
            return epsmu*(rt/(r*sr)*(cos(phi))^2+(r*sr/rt)*(sin(phi))^2);
        end
    end
    function pml_yy(x)
        r=hypot(x[1],x[2]);
        if r<=r_pml
            return ComplexF64(epsmu);
        else
            phi=atan(x[2],x[1]);
            rt=r+im*alphaPML/3.0*(r-r_pml)^3/(d_pml)^2;
            sr=1.0+im*alphaPML*(r-r_pml)^2/(d_pml)^2;
            return epsmu*(rt/(r*sr)*(sin(phi))^2+(r*sr/rt)*(cos(phi))^2);
        end
    end
    function pml_zz(x)
        r=hypot(x[1],x[2]);
        if r<=r_pml
            return ComplexF64(epsmu);
        else
            rt=r+im*alphaPML/3.0*(r-r_pml)^3/(d_pml)^2;
            sr=1.0+im*alphaPML*(r-r_pml)^2/(d_pml)^2;
            return epsmu*(rt/r)*sr;
        end
    end
    function pml_xy(x)
        r=hypot(x[1],x[2]);
        if r<=r_pml
            return ComplexF64(0.0);
        else
            phi=atan(x[2],x[1]);
            rt=r+im*alphaPML/3.0*(r-r_pml)^3/(d_pml)^2;
            sr=1.0+im*alphaPML*(r-r_pml)^2/(d_pml)^2;
            return epsmu*sin(phi)*cos(phi)*(rt/(r*sr)-r*sr/rt);
        end
    end
    return tensor3(pml_xx,pml_xy,0,pml_xy,pml_yy,0,0,0,pml_zz);
end

function add_cylindrical_PML(epsmu::Function,r_pml::Real,d_pml::Real,alphaPML::Real)
    function pml_xx(x)
        r=hypot(x[1],x[2]);
        if r<=r_pml
            return ComplexF64(epsmu(x));
        else
            phi=atan(x[2],x[1]);
            rt=r+im*alphaPML/3.0*(r-r_pml)^3/(d_pml)^2;
            sr=1.0+im*alphaPML*(r-r_pml)^2/(d_pml)^2;
            return epsmu(x)*(rt/(r*sr)*(cos(phi))^2+(r*sr/rt)*(sin(phi))^2);
        end
    end
    function pml_yy(x)
        r=hypot(x[1],x[2]);
        if r<=r_pml
            return ComplexF64(epsmu(x));
        else
            phi=atan(x[2],x[1]);
            rt=r+im*alphaPML/3.0*(r-r_pml)^3/(d_pml)^2;
            sr=1.0+im*alphaPML*(r-r_pml)^2/(d_pml)^2;
            return epsmu(x)*(rt/(r*sr)*(sin(phi))^2+(r*sr/rt)*(cos(phi))^2);
        end
    end
    function pml_zz(x)
        r=hypot(x[1],x[2]);
        if r<=r_pml
            return ComplexF64(epsmu(x));
        else
            rt=r+im*alphaPML/3.0*(r-r_pml)^3/(d_pml)^2;
            sr=1.0+im*alphaPML*(r-r_pml)^2/(d_pml)^2;
            return epsmu(x)*(rt/r)*sr;
        end
    end
    function pml_xy(x)
        r=hypot(x[1],x[2]);
        if r<=r_pml
            return ComplexF64(0.0);
        else
            phi=atan(x[2],x[1]);
            rt=r+im*alphaPML/3.0*(r-r_pml)^3/(d_pml)^2;
            sr=1.0+im*alphaPML*(r-r_pml)^2/(d_pml)^2;
            return epsmu(x)*sin(phi)*cos(phi)*(rt/(r*sr)-r*sr/rt);
        end
    end
    return tensor3(pml_xx,pml_xy,0,pml_xy,pml_yy,0,0,0,pml_zz);
end

function add_cylindrical_PML(epsmu::tensor3,r_pml::Real,d_pml::Real,alphaPML::Real)
    r(x)=hypot(x[1],x[2]);
    Sr(x)=1.0+im*alphaPML*(r(x)-r_pml)^2/(d_pml)^2;
    Sphi(x)=1.0+im*alphaPML/3.0*(r(x)-r_pml)^3/(d_pml)^2/r(x);
    C(x)=cos(atan(x[2],x[1]));
    S(x)=sin(atan(x[2],x[1]));
    Txx(x)=(r(x)<r_pml) ? ComplexF64(1.0) : C(x)^2/Sr(x)+S(x)^2/Sphi(x);
    Txy(x)=(r(x)<r_pml) ? ComplexF64(0.0) : (1/Sr(x)-1/Sphi(x))*C(x)*S(x);
    Tyy(x)=(r(x)<r_pml) ? ComplexF64(1.0) : S(x)^2/Sr(x)+C(x)^2/Sphi(x);
    T1=tensor3(Txx,Txy,0,Txy,Tyy,0,0,0,1);
    Txx2(x)=(r(x)<r_pml) ? ComplexF64(1.0) : C(x)^2*Sphi(x)+S(x)^2*Sr(x);
    Txy2(x)=(r(x)<r_pml) ? ComplexF64(0.0) : (Sphi(x)-Sr(x))*C(x)*S(x);
    Tyy2(x)=(r(x)<r_pml) ? ComplexF64(1.0) : S(x)^2*Sphi(x)+C(x)^2*Sr(x);
    Tzz2(x)=(r(x)<r_pml) ? ComplexF64(1.0) : Sr(x)*Sphi(x)
    T2=tensor3(Txx2,Txy2,0,Txy2,Tyy2,0,0,0,Tzz2);
    return (T1*epsmu*T2)
end

"""
    add_rectangular_PML(epsmu::Union{Number,Function,tensor3},x_pml::Real,dx_pml::Real,y_pml::Real,dy_pml::Real,alpha::Real)

Function that returns a tensor of permittivity/permeability with a rectangular PML located at [x_pml,x_pml+dx_pml], [-x_pml-dx_pml,-x_pml], [y_pml,y_pml+dy_pml], [-y_pml-dy_pml,-y_pml]

- epsmu: permittivity/permeability profile of the fiber. The fiber is assumed to be isotropic when epsmu is a `Number` or a `Function`
- x_pml: distance between the fiber center and the PML beginning in the x direction
- dx_pml: PML thickness in the x direction
- y_pml: distance between the fiber center and the PML beginning in the y direction
- dy_pml: PML thickness in the y direction
- alphaPML: attenuation factor of the PML
"""
function add_rectangular_PML(epsmu::Union{Number,Function,tensor3},x_pml::Real,dx_pml::Real,y_pml::Real,dy_pml::Real,alphaPML::Real)
    Gridap.Helpers.@abstractmethod
end

function add_rectangular_PML(epsmu::Number,x_pml::Real,dx_pml::Real,y_pml::Real,dy_pml::Real,alphaPML::Real)
    Sx=x->1.0+im*alphaPML*(max(abs(x[1])-x_pml,0))^2/(dx_pml)^2;
    Sy=x->1.0+im*alphaPML*(max(abs(x[2])-y_pml,0))^2/(dy_pml)^2;
    return tensor3(x->epsmu*Sy(x)/Sx(x),x->epsmu*Sx(x)/Sy(x),x->epsmu*Sx(x)*Sy(x));
end

function add_rectangular_PML(epsmu::Function,x_pml::Real,dx_pml::Real,y_pml::Real,dy_pml::Real,alphaPML::Real)
    Sx=x->1.0+im*alphaPML*(max(abs(x[1])-x_pml,0))^2/(dx_pml)^2;
    Sy=x->1.0+im*alphaPML*(max(abs(x[2])-y_pml,0))^2/(dy_pml)^2;
    return tensor3(x->epsmu(x)*Sy(x)/Sx(x),x->epsmu(x)*Sx(x)/Sy(x),x->epsmu(x)*Sx(x)*Sy(x));
end

function add_rectangular_PML(epsmu::tensor3,x_pml::Real,dx_pml::Real,y_pml::Real,dy_pml::Real,alphaPML::Real)
    Sx=x->1.0+im*alphaPML*(max(abs(x[1])-x_pml,0))^2/(dx_pml)^2;
    Sy=x->1.0+im*alphaPML*(max(abs(x[2])-y_pml,0))^2/(dy_pml)^2;
    S1=tensor3(x->1/Sx(x),x->1/Sy(x),1);
    S2=tensor3(Sy,Sx,x->Sx(x)*Sy(x))
    return (S1*epsmu*S2)
end


"""
    add_rectangular_PML(epsmu::Union{Number,Function,tensor3},xm_pml::Real,xp_pml::Real,dx_pml::Real,ym_pml::Real,yp_pml::Real,dy_pml::Real,alphaPML::Real)

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
function add_rectangular_PML(epsmu::Union{Number,Function,tensor3},xm_pml::Real,xp_pml::Real,dx_pml::Real,ym_pml::Real,yp_pml::Real,dy_pml::Real,alphaPML::Real)
    Gridap.Helpers.@abstractmethod
end

function add_rectangular_PML(epsmu::Number,xm_pml::Real,xp_pml::Real,dx_pml::Real,ym_pml::Real,yp_pml::Real,dy_pml::Real,alphaPML::Real)
    Sx=x->1.0+im*alphaPML*(max(x[1]-xp_pml,0))^2/(dx_pml)^2+im*alphaPML*(max(xm_pml-x[1],0))^2/(dx_pml)^2;
    Sy=x->1.0+im*alphaPML*(max(x[2]-yp_pml,0))^2/(dy_pml)^2+im*alphaPML*(max(ym_pml-x[2],0))^2/(dy_pml)^2;
    return tensor3(x->epsmu*Sy(x)/Sx(x),x->epsmu*Sx(x)/Sy(x),x->epsmu*Sx(x)*Sy(x));
end

function add_rectangular_PML(epsmu::Function,xm_pml::Real,xp_pml::Real,dx_pml::Real,ym_pml::Real,yp_pml::Real,dy_pml::Real,alphaPML::Real)
    Sx=x->1.0+im*alphaPML*(max(x[1]-xp_pml,0))^2/(dx_pml)^2+im*alphaPML*(max(xm_pml-x[1],0))^2/(dx_pml)^2;
    Sy=x->1.0+im*alphaPML*(max(x[2]-yp_pml,0))^2/(dy_pml)^2+im*alphaPML*(max(ym_pml-x[2],0))^2/(dy_pml)^2;
    return tensor3(x->epsmu(x)*Sy(x)/Sx(x),x->epsmu(x)*Sx(x)/Sy(x),x->epsmu(x)*Sx(x)*Sy(x));
end

function add_rectangular_PML(epsmu::tensor3,xm_pml::Real,xp_pml::Real,dx_pml::Real,ym_pml::Real,yp_pml::Real,dy_pml::Real,alphaPML::Real)
    Sx=x->1.0+im*alphaPML*(max(x[1]-xp_pml,0))^2/(dx_pml)^2+im*alphaPML*(max(xm_pml-x[1],0))^2/(dx_pml)^2;
    Sy=x->1.0+im*alphaPML*(max(x[2]-yp_pml,0))^2/(dy_pml)^2+im*alphaPML*(max(ym_pml-x[2],0))^2/(dy_pml)^2;
    S1=tensor3(x->1/Sx(x),x->1/Sy(x),x->1);
    S2=tensor3(Sy,Sx,x->Sx(x)*Sy(x))
    return (S1*epsmu*S2)
end


"""
    add_twist_PML(epsmu::Union{Function,Number},P::Real,r_pml::Real,d_pml::Real,alphaPML::Real)

Function that returns a tensor of permittivity/permeability of a twisted isotropic fiber with a cylindrical PML

- epsmu: function of the tuple (x,y) that describes the permittivity/permeability profile of the fiber
- P: period of the twist
- r_pml: distance between the fiber center and the PML beginning
- d_pml: PML thickness
- alpha: attenuation factor of the PML
"""
function add_twist_PML(epsmu::Union{Number,Function},P::Real,r_pml::Real,d_pml::Real,alphaPML::Real)
    Gridap.Helpers.@abstractmethod
end

function add_twist_PML(epsmu::Function,P::Real,r_pml::Real,d_pml::Real,alphaPML::Real)
    alpha=2*pi/P;
    r(x)=hypot(x[1],x[2]);
    Sr(x)=1.0+im*alphaPML*(r(x)-r_pml)^2/(d_pml)^2;
    rp(x)=r(x)+im*alphaPML/3.0*(r(x)-r_pml)^3/(d_pml)^2;
    phi(x)=2*atan(x[2]/(x[1]+r(x)))
    C(x)=cos(phi(x));
    S(x)=sin(phi(x));
    S2(x)=sin(2*phi(x));
    invTxx(x)=(r(x)<r_pml) ? ComplexF64((1.0+alpha^2*x[2]^2)*epsmu(x)) : (C(x)^2/Sr(x)*rp(x)/r(x)+S(x)^2*Sr(x)*r(x)/rp(x)*(1+alpha^2*rp(x)^2))*epsmu(x);
    invTxy(x)=(r(x)<r_pml) ? ComplexF64(-alpha^2*x[1]*x[2]*epsmu(x)) : S2(x)*(rp(x)^2-r(x)^2*(1+alpha^2*rp(x)^2)*Sr(x)^2)/(2*r(x)*rp(x)*Sr(x))*epsmu(x)
    invTxz(x)=(r(x)<r_pml) ? ComplexF64(-alpha*x[2]*epsmu(x)) : -alpha*rp(x)*Sr(x)*S(x)*epsmu(x)
    invTyy(x)=(r(x)<r_pml) ? ComplexF64((1.0+alpha^2*x[1]^2)*epsmu(x)) : (S(x)^2/Sr(x)*rp(x)/r(x)+C(x)^2*Sr(x)*r(x)/rp(x)*(1+alpha^2*rp(x)^2))*epsmu(x);
    invTyz(x)=(r(x)<r_pml) ? ComplexF64(alpha*x[1]*epsmu(x)) : alpha*rp(x)*Sr(x)*C(x)*epsmu(x)
    invTzz(x)=(r(x)<r_pml) ? ComplexF64(epsmu(x)) : rp(x)*Sr(x)/r(x)*epsmu(x)
    return tensor3(invTxx,invTxy,invTxz,invTxy,invTyy,invTyz,invTxz,invTyz,invTzz);
end

function add_twist_PML(epsmu::Number,P::Real,r_pml::Real,d_pml::Real,alphaPML::Real)
    alpha=2*pi/P;
    r(x)=hypot(x[1],x[2]);
    Sr(x)=1.0+im*alphaPML*(r(x)-r_pml)^2/(d_pml)^2;
    rp(x)=r(x)+im*alphaPML/3.0*(r(x)-r_pml)^3/(d_pml)^2;
    phi(x)=2*atan(x[2]/(x[1]+r(x)))
    C(x)=cos(phi(x));
    S(x)=sin(phi(x));
    S2(x)=sin(2*phi(x));
    invTxx(x)=(r(x)<r_pml) ? ComplexF64((1.0+alpha^2*x[2]^2)*epsmu) : (C(x)^2/Sr(x)*rp(x)/r(x)+S(x)^2*Sr(x)*r(x)/rp(x)*(1+alpha^2*rp(x)^2))*epsmu;
    invTxy(x)=(r(x)<r_pml) ? ComplexF64(-alpha^2*x[1]*x[2]*epsmu) : S2(x)*(rp(x)^2-r(x)^2*(1+alpha^2*rp(x)^2)*Sr(x)^2)/(2*r(x)*rp(x)*Sr(x))*epsmu
    invTxz(x)=(r(x)<r_pml) ? ComplexF64(-alpha*x[2]*epsmu) : -alpha*rp(x)*Sr(x)*S(x)*epsmu
    invTyy(x)=(r(x)<r_pml) ? ComplexF64((1.0+alpha^2*x[1]^2)*epsmu) : (S(x)^2/Sr(x)*rp(x)/r(x)+C(x)^2*Sr(x)*r(x)/rp(x)*(1+alpha^2*rp(x)^2))*epsmu;
    invTyz(x)=(r(x)<r_pml) ? ComplexF64(alpha*x[1]*epsmu) : alpha*rp(x)*Sr(x)*C(x)*epsmu
    invTzz(x)=(r(x)<r_pml) ? ComplexF64(epsmu) : rp(x)*Sr(x)/r(x)*epsmu
    return tensor3(invTxx,invTxy,invTxz,invTxy,invTyy,invTyz,invTxz,invTyz,invTzz);
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

"""
    nb_args(f::Function)

Returns the possible numbers of arguments of a function

"""
function nb_args(f::Function)
    return unique([methods(f)[i].nargs-1 for i in axes(methods(f),1)])
end

############################# integration #############################
"""
    integrate1D(f::Function;rtol::Real=1E-5;kwargs...)

Integrate a 1D function f(r) from 0 to infinity

The keyword arguments are those of HCubature.jl. In particular,`rtol` and `atol` are the relative and the absolute tolerance on the result
"""
function integrate1D(f::Function;kwargs...)
    f2=x->f(x[1]/(1-x[1]))/(1-x[1])^2
    return hcubature(f2,[0.0],[1.0];kwargs...)
end

"""
    integrate2D(f::Function;kwargs...)

Integrate a 2D function f(r) with r is a the tuple (x,y) for x and y from -infinity to infinity

    The keyword arguments are those of HCubature.jl. In particular,`rtol` and `atol` are the relative and the absolute tolerance on the result
"""
function integrate2D(f::Function;kwargs...)
    f2=x->f(Point(x[1]/(1-x[1]^2),x[2]/(1-x[2]^2)))*(1+x[1]^2)/(1-x[1]^2)^2*(1+x[2]^2)/(1-x[2]^2)^2
    return hcubature(f2,[-1.0,-1.0],[1.0,1.0];kwargs...)
end