
#FD(1E-6,0,1,x->1.48-0.03*(x>2E-6),1000,10E-6)
#lambda=collect(1:0.1:1.5);
#f=[x->Ge_DopedSilica_Fleming(l,0)+0.05*(1-(x>2)) for l in lambda];
#aa=FD.(lambda,0,2,f,1000,10);
##bb=Vector{Vector{ScalarMode1D}}(undef,length(lambda));
##for i=1:length(lambda)
##    bb[i]=FD(lambda[i],0,2,f[i],1000,10);
##end
#Note pour mettre un vecteur en argument sans vectoriser, utiliser Ref
"""
    FD(lambda::Real,l::Integer,mmax::Integer,fonc::Function,nb::Integer,rmax::Real;field::Bool=false,order::Integer=1,solver::Symbol=:Arpack,tol::Float64=0.0)

Returns a vector of `ScalarMode1D` in the case of a cylindrically-symmetric weakly-guiding fiber.  

- lambda: wavelength
- l: azimuthal number
- mmax: maximal number of modes (useful if the fiber is very multimode)
- fonc: function of the radial coordinate r that describes the refractive index profile
- nb: number of nodes for the finite difference method
- rmax: maximal value of r
- field: boolean that indicates if fields must be saved
- order: order of the method to average the refractive index profile for each node (no averaging if order=1)
- solver: can be :Arpack, :LU or :MUMPS
- tol: tolerance for the eigenvalue solver (see documention of Arpack.jl or ArnoldiMethod.jl)
"""
function FD(lambda::Real,l::Integer,mmax::Integer,fonc::Function,nb::Integer,rmax::Real;field::Bool=false,order::Integer=1,solver::Symbol=:Arpack,tol::Float64=0.0)
    if (lambda<=0)
        throw(DomainError(lambda, "Wavelength must be strictly positive"));
    end
    if (l<0)
        throw(DomainError(nu, "The azimuthal number must be positive"));
    end
    if (mmax<1)
        throw(DomainError(mmax, "The number of modes must be strictly positive"));
    end
    if (rmax<=0)
        throw(DomainError(rmax, "The maximum radius must be strictly positive"));
    end
    if (nb<2)
        throw(DomainError(nb, "The number of points must be greater than 2"));
    end
    if (order<1)
        throw(DomainError(order, "The Gauss-legendre order must be greater than 1"));
    end
    if (first(methods(fonc)).nargs!=2)
        throw(DomainError(fonc, "The refractive index function must have 1 argument"));
    end
    if (tol<0.0)
        throw(DomainError(tol, "The tolerance must be positive or null"));
    end
    h=rmax/(nb-1.0);
    if (l>0)
        U=LinRange(h,rmax,nb);
    else
        U=LinRange(h/2,rmax+h/2,nb);
    end
    n=indexGaussSampling1D(U,order,fonc);
    n_max=maximum(n);
    n_cladding=n[end];
    U=U/rmax;
    h=h/rmax;
    k0=2*pi/lambda;
    A=@. (rmax^2*k0^2*n^2)-((l^2)/(U^2))-(2/h^2);
    B=@. 0.5/h/U[1:(end-1)]+1/h^2;
    C=@. -0.5/h/U[2:end]+1/h^2;
    S=Tridiagonal(C,A,B);
    if (solver==:Arpack)
        D,V=eigs(S,nev=mmax,sigma=(k0*n_max*rmax)^2,tol=tol);
    elseif (solver==:LU)
        D,V=eigs_LU(S,nev=mmax,sigma=(k0*n_max*rmax)^2,tol=tol);
    elseif (solver==:MUMPS)
        D,V=eigs_MUMPS(S,nev=mmax,sigma=(k0*n_max*rmax)^2,tol=tol);
    else
        throw(DomainError(solver, "solver must be :Arpack, :LU or :MUMPS"));
    end
    neff=sqrt.(real(D)/k0^2/rmax^2);
    if (field)
        V=real(V[:,neff.>n_cladding]);
        neff=neff[neff.>n_cladding];
        if (l>0)
            r=[0;U*rmax];
            V=[zeros(1,length(neff));V];
        else
            U2=(U+circshift(U,-1))/2;
            r=[0;U2[1:end-1]*rmax];
            V2=(V+circshift(V,-1))/2;
            V=[V[1,:]';V2[1:end-1,:]];
        end
        modes=Vector{ScalarMode1D}(undef,length(neff));
        for i ∈ axes(neff,1)
            name=string("LP ",string(l),",",string(i));
            modes[i]=ScalarMode1D(name,neff[i],lambda,l,r,V[:,i]);
        end
    else
        neff=neff[neff.>n_cladding];
        modes=Vector{ScalarMode1D}(undef,length(neff));
        for i ∈ axes(neff,1)
            name=string("LP ",string(l),",",string(i));
            modes[i]=ScalarMode1D(name,neff[i],lambda,l,[],[]);
        end
    end
    return modes;
end

"""
    FD(lambda::Real,mmax::Integer,fonc::Function,nbx::Integer,nby::Integer,xmax::Real,ymax::Real;field::Bool=false,order::Integer=1,type::Symbol=:Scalar,solver::Symbol=:Arpack,tol::Float64=0.0)

Returns a vector of `ScalarMode2D` if type=:Scalar or a vector of `VectorMode` if type=:Vector.  

- lambda: wavelength
- mmax: maximal number of modes (useful if the fiber is very multimode)
- fonc: function of the cartesian coordinates x/y that describes the refractive index profile
- nbx: number of nodes for the finite difference method in the direction x
- nby: number of nodes for the finite difference method in the direction y
- xmax: maximal value of x
- ymax: maximal value of x
- field: boolean that indicates if fields must be saved
- order: order of the method to average the refractive index profile for each node (no averaging if order=1)
- type: must be :Scalar or :Vector
- solver: can be :Arpack, :LU or :MUMPS
- tol: tolerance for the eigenvalue solver (see documention of Arpack.jl or ArnoldiMethod.jl)
"""
function FD(lambda::Real,mmax::Integer,fonc::Function,nbx::Integer,nby::Integer,xmax::Real,ymax::Real;field::Bool=false,order::Integer=1,type::Symbol=:Scalar,solver::Symbol=:Arpack,tol::Float64=0.0)
    if(type==:Scalar)
        scalarFD(lambda,mmax,fonc,nbx,nby,xmax,ymax,field=field,order=order,solver=solver,tol=tol);
    elseif (type==:Vector)
        vectorFD(lambda,mmax,fonc,nbx,nby,xmax,ymax,field=field,order=order,solver=solver,tol=tol);
    else
        throw(DomainError(type, "type must be :Scalar or :Vector"));
    end
end



#scalarFD(1,5,(X,Y)->1.48-0.03*(sqrt(X^2+Y^2)>2),500,500,5,5,order=5)
function scalarFD(lambda::Real,nbmax::Integer,fonc::Function,nbx::Integer,nby::Integer,xmax::Real,ymax::Real;field::Bool=false,order::Integer=1,solver::Symbol=:Arpack,tol::Float64=0.0)
    if (lambda<=0)
        throw(DomainError(lambda, "Wavelength must be strictly positive"));
    end
    if (nbmax<1)
        throw(DomainError(nbmax, "The number of modes must be strictly positive"));
    end
    if (xmax<=0)
        throw(DomainError(xmax, "The maximum X value must be strictly positive"));
    end
    if (ymax<=0)
        throw(DomainError(ymax, "The maximum Y value must be strictly positive"));
    end
    if (nbx<2)
        throw(DomainError(nbx, "The number of point along x must be greater than 2"));
    end
    if (nby<2)
        throw(DomainError(nbx, "The number of point along y must be greater than 2"));
    end
    if (order<1)
        throw(DomainError(order, "The Gauss-legendre order must be greater than 1"));
    end
    if (first(methods(fonc)).nargs!=3)
        throw(DomainError(fonc, "The refractive index function must have 2 arguments"));
    end
    if (tol<0.0)
        throw(DomainError(tol, "The tolerance must be positive or null"));
    end
    Nx=2*nbx-1;
    Ny=2*nby-1;
    x,y,n=indexGaussSampling2D(xmax,ymax,Nx,Ny,order,fonc);
    dx=step(x);
    dy=step(y);
    nmin=max(maximum(n[1,:]),maximum(n[end,:]),maximum(n[:,1]),maximum(n[:,end]));
    nmax=maximum(n);
    if (nmin>=nmax)
        return Vector{scalarMode2D}[];
    end
    k0=2*pi/lambda;
    xymax=max(xmax,ymax);
    M=Laplacian(Nx,Ny,dx/xymax,dy/xymax);
    M=-M+spdiagm(0=>dropdims(reshape((n).^2,Nx*Ny,1)*(k0*xymax)^2,dims=2));
    if (solver==:Arpack)
        D,V=eigs(M,nev=nbmax,sigma=(k0*nmax*xymax)^2,tol=tol);
    elseif (solver==:LU)
        D,V=eigs_LU(M,nev=nbmax,sigma=(k0*nmax*xymax)^2,tol=tol);
    elseif (solver==:MUMPS)
        D,V=eigs_MUMPS(M,nev=nbmax,sigma=(k0*nmax*xymax)^2,tol=tol);
    else
        throw(DomainError(solver, "solver must be :Arpack, :LU or :MUMPS"));
    end
    neff=sqrt.(real(D)/k0^2/xymax^2);
    if (field)
        V=real(V[:,neff.>nmin]);
        neff=neff[neff.>nmin];
        modes=Vector{ScalarMode2D}(undef,length(neff));
        for i in eachindex(neff)
            name=string("LP - mode ",string(i));
            modes[i]=ScalarMode2D(name,neff[i],lambda,x,y,(reshape(V[:,i],Nx,Ny)));
        end
    else
        neff=neff[neff.>nmin];
        modes=Vector{ScalarMode2D}(undef,length(neff));
        empty_array=Array{Float64}(undef,0,2);
        for i in eachindex(neff)
            name=string("LP - mode ",string(i));
            modes[i]=ScalarMode2D(name,neff[i],lambda,[],[],empty_array);
        end
    end
    return modes;
end

#vectorFD(1,5,(X,Y)->1.48-0.03*(sqrt(X^2+Y^2)>2),500,500,5,5,order=5)
function vectorFD(lambda::Real,nbmax::Integer,fonc::Function,nbx::Integer,nby::Integer,xmax::Real,ymax::Real;field::Bool=false,order::Integer=1,solver::Symbol=:Arpack,tol::Float64=0.0)
    if (xmax<=0)
        throw(DomainError(xmax, "The maximum X value must be strictly positive"));
    end
    if (ymax<=0)
        throw(DomainError(ymax, "The maximum Y value must be strictly positive"));
    end
    if (nbx<2)
        throw(DomainError(nbx, "The number of point along x must be greater than 2"));
    end
    if (nby<2)
        throw(DomainError(nby, "The number of point along y must be greater than 2"));
    end
    if (order<1)
        throw(DomainError(order, "The Gauss-legendre order must be greater than 1"));
    end
    if (first(methods(fonc)).nargs!=3)
        throw(DomainError(fonc, "The refractive index function must have 2 arguments"));
    end
    if (tol<0.0)
        throw(DomainError(tol, "The tolerance must be positive or null"));
    end
    Nx=2*nbx-1;
    Ny=2*nby-1;
    x,y,n=indexGaussSampling2D(xmax,ymax,Nx,Ny,order,fonc);
    dx=x[2]-x[1];
    dy=y[2]-y[1];
    nmin=max(maximum(n[1,:]),maximum(n[end,:]),maximum(n[:,1]),maximum(n[:,end]));
    nmax=maximum(n);
    if (nmin>=nmax)
        return Vector{scalarMode2D}[];
    end
    k0=2*pi/lambda;
    xymax=max(xmax,ymax);
    eps=n.^2;
    epsz_ini=(circshift(eps,(0,1))+circshift(eps,(1,0))+eps+circshift(eps,(1,1)))/4.0;
    invepsz=spdiagm(0=>dropdims(1.0./reshape(epsz_ini,Nx*Ny,1),dims=2));
    epsx_ini=(eps+circshift(eps,(0,1)))/2.0;
    epsy_ini=(eps+circshift(eps,(1,0)))/2.0;
    epsx=spdiagm(0=>dropdims(reshape(epsx_ini,Nx*Ny,1),dims=2));
    epsy=spdiagm(0=>dropdims(reshape(epsy_ini,Nx*Ny,1),dims=2));
    
    Ux=spdiagm(0=>-ones(Nx*Ny),1=>ones(Nx*Ny-1))/dx*xymax;
    Vx=-copy(Transpose(Ux));
    Uy=spdiagm(0=>-ones(Nx*Ny),Nx=>ones(Nx*Ny-Nx))/dy*xymax;
    Vy=-copy(Transpose(Uy));

    k0=k0*xymax;

    Pxx=-Ux*invepsz*Vy*Vx*Uy/k0^2+(k0^2*sparse(I,Nx*Ny,Nx*Ny)+Ux*invepsz*Vx)*(epsx+Vy*Uy/k0^2);
    Pyy=-Uy*invepsz*Vx*Vy*Ux/k0^2+(k0^2*sparse(I,Nx*Ny,Nx*Ny)+Uy*invepsz*Vy)*(epsy+Vx*Ux/k0^2);
    Pxy=Ux*invepsz*Vy*(epsy+Vx*Ux/k0^2)-(k0^2*sparse(I,Nx*Ny,Nx*Ny)+Ux*invepsz*Vx)/k0^2*Vy*Ux;
    Pyx=Uy*invepsz*Vx*(epsx+Vy*Uy/k0^2)-(k0^2*sparse(I,Nx*Ny,Nx*Ny)+Uy*invepsz*Vy)/k0^2*Vx*Uy;
    P=[Pxx Pxy ; Pyx Pyy];
    if (solver==:Arpack)
        D,V=eigs(P,nev=nbmax,sigma=(k0*nmax)^2,tol=tol);
    elseif (solver==:LU)
        D,V=eigs_LU(P,nev=nbmax,sigma=(k0*nmax)^2,tol=tol);
    elseif (solver==:MUMPS)
        D,V=eigs_MUMPS(P,nev=nbmax,sigma=(k0*nmax)^2,tol=tol);
    else
        throw(DomainError(solver, "solver must be :Arpack, :LU or :MUMPS"));
    end
    print(typeof(V));
    neff=sqrt.(real(D)/k0^2);
    if (field)
        V=real(V[:,neff.>nmin]);
        neff=neff[neff.>nmin];
        modes=Vector{VectorMode}(undef,length(neff));
        for i=1:length(neff)
            name=string("Mode ",string(i));
            typeof(V);
            Ex=(reshape(V[1:Nx*Ny,i],(Nx,Ny)));
            Ey=(reshape(V[Nx*Ny+1:2*Nx*Ny,i],(Nx,Ny)));
            Hz=((circshift(Ey,(-1,0))-Ey)/dx-(circshift(Ex,(0,-1))-Ex)/dy)/k0/im;
            Hy=epsx_ini.*Ex/neff[i]+im*(circshift(Hz,(0,1))-Hz)/dy/k0/neff[i];
            Hx=-epsy_ini.*Ey/neff[i]+im*(circshift(Hz,(1,0))-Hz)/dx/k0/neff[i];
            Ez=-((Hy-circshift(Hy,(1,0)))/dx-(Hx-circshift(Hx,(0,1)))/dy)/im/k0./epsz_ini;
            modes[i]=VectorMode(name,neff[i],lambda,x,y,Z0*Ex,Z0*Ey,imag.(Ez)*Z0,real.(Hx),real.(Hy),imag.(Hz));
        end
    else
        neff=neff[neff.>nmin];
        modes=Vector{VectorMode}(undef,length(neff));
        for i=1:length(neff)
            name=string("Mode ",string(i));
            modes[i]=VectorMode(name,neff[i],lambda);
        end
    end
    return modes;
end