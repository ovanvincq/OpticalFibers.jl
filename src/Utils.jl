using LinearAlgebra
using LinearMaps
using MUMPS
using MPI
using SparseArrays
using ArnoldiMethod
using FastGaussQuadrature;

export piecewiseIndex
export meshgrid
export derivative
export ring
export tensor3
export inverse

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

function sdiff1(M::Integer)
    return sparse([ [1.0 zeros(1,M-1)]; diagm(1=>ones(M-1)) - Matrix(I,M,M) ])
end

function Laplacian(Nx::Integer, Ny::Integer, dx::Real, dy::Real)
   Dx = sdiff1(Nx) / dx
   Dy = sdiff1(Ny) / dy
   Ax = Dx' * Dx
   Ay = Dy' * Dy
   return kron(sparse(I,Ny,Ny),Ax) + kron(Ay,sparse(I,Nx,Nx))
end

function indexGaussSampling1D(r::Union{StepRangeLen{<:Real},LinRange{<:Real}},order::Integer,f::Function)
    if (order==1);
        return f.(r);
    end
    dr=step(r);
    permittivity=x->f(abs(x))^2;
    nodes,weight=gausslegendre(order);
    n=sqrt.(vec(sum(transpose(weight).*((permittivity.(transpose(nodes*dr/2).+r))),dims=2))/2);
    return n;
end

function indexGaussSampling1D(rmax::Real,nb::Integer,order::Integer,f::Function)
    r=LinRange(0,rmax,nb);
    n=indexGaussSampling1D(r,order,f);
    return r,n;
end

function indexGaussSampling2D(xmax::Real,ymax::Real,nbx::Integer,nby::Integer,order::Integer,f::Function)
    x=LinRange(-xmax,xmax,nbx);
    y=LinRange(-ymax,ymax,nby);
    dx=step(x);
    dy=step(y);
    permittivity=(x,y)->f(x,y)^2;
    nodes,weight=gausslegendre(order);
    w=weight*weight';
    n=zeros(nbx,nby);
    Threads.@threads for j in eachindex(x)
        for k in eachindex(y)
             n[j,k]=sqrt(sum(w.*permittivity.((nodes*dx/2).+x[j],((nodes*dy/2).+y[k])'))/4);
        end
    end
    return x,y,n;
end

"""
    derivative(t::Tuple{Vector{<:Real},Vector{<:Number}})

Function that returns the derivative of the tuple t=(x,y)
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

Function that returns the nth order derivative of the tuple t=(x,y) 
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

function eigs_MUMPS(A,B;sigma=0,nev::Int64=1,tol::Float64=0.0,restarts::Int64=100)
    if (!MPI.Initialized())
        MPI.Init();
    end
    icntl = get_icntl(verbose=false)
    m = Mumps{eltype(A)}(mumps_unsymmetric, icntl, default_cntl64)
    associate_matrix!(m,SparseMatrixCSC(A-sigma*B));
    factorize!(m);
    a = ShiftAndInvert_MUMPS(m,B,Vector{eltype(B)}(undef, size(B,1)))
    map_mumps=LinearMap{eltype(A)}(a, size(A,1), ismutating=true)
    decomp,  = partialschur(map_mumps, nev=nev, tol=tol, restarts=restarts, which=LM())
    λs_inv, X = partialeigen(decomp);
    MUMPS.finalize!(m);
    λs=(1 ./λs_inv).+sigma
    return λs,X;
end

function eigs_MUMPS(A;sigma=0,nev::Int64=1,tol::Float64=0.0,restarts::Int64=100)
    eigs_MUMPS(A,eltype(A).(SparseMatrixCSC(I,size(A,1),size(A,2)));sigma=sigma,nev=nev,tol=tol,restarts=restarts)
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

function eigs_LU(A,B;sigma=0,nev::Int64=1,tol::Float64=0.0,restarts::Int64=100)
    a = ShiftAndInvert_LU(lu(A-sigma*B),B,Vector{eltype(A)}(undef, size(A,1)))
    map_LU=LinearMap{eltype(A)}(a, size(A,1), ismutating=true)
    decomp,  = partialschur(map_LU, nev=nev, tol=tol, restarts=restarts, which=LM())
    λs_inv, X = partialeigen(decomp);
    λs=(1 ./λs_inv).+sigma
    return λs,X;
end

function eigs_LU(A;sigma=0,nev::Int64=1,tol::Float64=0.0,restarts::Int64=100)
    eigs_LU(A,SparseMatrixCSC(I,size(A,1),size(A,2));sigma=sigma,nev=nev,tol=tol,restarts=restarts)
end

"""
    tensor

Structure that describes a tensor of permittivity or permeability. Each of the 9 components of the a `tensor3` is a function of a tuple (x,y).

- xx :: `Function`
- yx :: `Function`
- zx :: `Function`
- xy :: `Function`
- yy :: `Function`
- zy :: `Function`
- xz :: `Function`
- yz :: `Function`
- zz :: `Function`
"""
struct tensor3
    xx::Function
    yx::Function
    zx::Function
    xy::Function
    yy::Function
    zy::Function
    xz::Function
    yz::Function
    zz::Function
end

"""
    tensor3(xx::Union{Function,Number},yx::Union{Function,Number},zx::Union{Function,Number},xy::Union{Function,Number},yy::Union{Function,Number},zy::Union{Function,Number},xz::Union{Function,Number},yz::Union{Function,Number},zz::Union{Function,Number})

Function that returns a tensor3 defined with functions and constant values. Each function must always return the same type.
"""
function tensor3(xx::Union{Function,Number},yx::Union{Function,Number},zx::Union{Function,Number},xy::Union{Function,Number},yy::Union{Function,Number},zy::Union{Function,Number},xz::Union{Function,Number},yz::Union{Function,Number},zz::Union{Function,Number})
    if (isa(xx,Number))
        xx2=x->xx;
    else
        xx2=xx;
    end
    if (isa(yx,Number))
        yx2=x->yx;
    else
        yx2=yx;
    end
    if (isa(zx,Number))
        zx2=x->zx;
    else
        zx2=zx;
    end
    if (isa(xy,Number))
        xy2=x->xy;
    else
        xy2=xy;
    end
    if (isa(yy,Number))
        yy2=x->yy;
    else
        yy2=yy;
    end
    if (isa(zy,Number))
        zy2=x->zy;
    else
        zy2=zy;
    end
    if (isa(xz,Number))
        xz2=x->xz;
    else
        xz2=xz;
    end
    if (isa(yz,Number))
        yz2=x->yz;
    else
        yz2=yz;
    end
    if (isa(zz,Number))
        zz2=x->zz;
    else
        zz2=zz;
    end
    tensor3(xx2,yx2,zx2,xy2,yy2,zy2,xz2,yz2,zz2)
end

"""
    tensor3(fxx::Union{Function,Number},fyy::Union{Function,Number},fzz::Union{Function,Number})

Function that returns a diagonal tensor3
"""
function tensor3(fxx::Union{Function,Number},fyy::Union{Function,Number},fzz::Union{Function,Number})
    return tensor3(fxx,0,0,0,fyy,0,0,0,fzz)
end

"""
    tensor3(f::Union{Function,Number})

Function that returns a diagonal tensor3 with the same function in the three diagonal terms
"""
function tensor3(f::Union{Function,Number})
    return tensor3(f,0,0,0,f,0,0,0,f)
end

"""
    inverse(t::tensor3)

Returns the inverse of the tensor t
"""
function inverse(t::tensor3)
    det(x)=(t.xx(x)*t.yy(x)*t.zz(x)+t.xy(x)*t.yz(x)*t.zx(x)+t.xz(x)*t.yx(x)*t.zy(x)-t.xx(x)*t.yz(x)*t.zy(x)-t.xy(x)*t.yx(x)*t.zz(x)-t.xz(x)*t.yy(x)*t.zx(x));
    xx(x)=(t.yy(x)*t.zz(x)-t.yz(x)*t.zy(x))/det(x)
    yx(x)=(t.yz(x)*t.zx(x)-t.yx(x)*t.zz(x))/det(x)
    zx(x)=(t.yx(x)*t.zy(x)-t.yy(x)*t.zx(x))/det(x)
    xy(x)=(t.xz(x)*t.zy(x)-t.xy(x)*t.zz(x))/det(x)
    yy(x)=(t.xx(x)*t.zz(x)-t.xz(x)*t.zx(x))/det(x)
    zy(x)=(t.xy(x)*t.zx(x)-t.xx(x)*t.zy(x))/det(x)
    xz(x)=(t.xy(x)*t.yz(x)-t.xz(x)*t.yy(x))/det(x)
    yz(x)=(t.xz(x)*t.yx(x)-t.xx(x)*t.yz(x))/det(x)
    zz(x)=(t.xx(x)*t.yy(x)-t.xy(x)*t.yx(x))/det(x)
    return tensor3(xx,yx,zx,xy,yy,zy,xz,yz,zz);
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

