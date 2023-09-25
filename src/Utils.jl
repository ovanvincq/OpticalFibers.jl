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

function Laplacian(Nx, Ny, dx, dy)
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


function derivative(t::Tuple{Vector{<:Real},Vector{<:Number}})
    if (length(t[1])!=length(t[2]))
        throw(ArgumentError("x and y must have the same length"));
    end
    xp=(t[1][1:end-1]+t[1][2:end])/2;
    yp=diff(t[2])./diff(t[1]);
    return xp,yp;
end

function derivative(t::Tuple{Vector{<:Real},Vector{<:Number}},order::Int64)
    temp=t;
    for j=1:order
        temp=derivative(temp);
    end
    return temp;
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
    eigs_MUMPS(A,Float64.(SparseMatrixCSC(I,size(A,1),size(A,2)));sigma=sigma,nev=nev,tol=tol,restarts=restarts)
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

function tensor3(fxx::Function,fyy::Function,fzz::Function)
    z=x->0;
    return tensor3(fxx,z,z,z,fyy,z,z,z,fzz)
end

function tensor3(f::Function)
    z=x->0;
    return tensor3(f,z,z,z,f,z,z,z,f)
end

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

