#model = GmshDiscreteModel(gmshFile);
"""
    FEM(lambda::Real,mmax::Int64,eps::Function,model::DiscreteModel,approx_neff::Real;order::Int64=2,field::Bool=false,type::Symbol=:Scalar,solver::Symbol=:LU,tol::Float64=0.0)

Returns a vector of `ScalarModeFEM` if type=:Scalar or a vector of `VectorModeFEM` if type=:Vector.  
The fiber is isotropic and described with its relative permittivity epsilon(pos) where pos is the tuple (x,y).

- lambda: wavelength
- mmax: maximal number of modes (useful if the fiber is very multimode)
- eps: function of the tuple (x,y) that describes the relative permittivity profile
- model: DiscreteModel generated with `GridapGmsh.jl`
- approx_neff: effective index around which the modes will be computed
- order: order of the FEM
- field: boolean that indicates if fields must be saved
- type: must be :Scalar or :Vector
- solver: can be :LU or :MUMPS
- tol: tolerance for the eigenvalue solver (see documention of ArnoldiMethod.jl)
"""
function FEM(lambda::Real,mmax::Int64,eps::Function,model::DiscreteModel,approx_neff::Real;order::Int64=2,field::Bool=false,type::Symbol=:Scalar,solver::Symbol=:LU,tol::Float64=0.0)
    if(type==:Scalar)
        FEM1(lambda,mmax,eps,model,approx_neff,order=order,field=field,solver=solver,tol=tol);
    elseif (type==:Vector)
        FEM2(lambda,mmax,eps,model,approx_neff,order=order,field=field,solver=solver,tol=tol);
    else
        throw(DomainError(type, "type must be :Scalar or :Vector"));
    end
end


function FEM1(lambda::Real,mmax::Int64,eps::Function,model::DiscreteModel,approx_neff::Real;order::Int64=2,field::Bool=false,solver::Symbol=:LU,tol::Float64=0.0)
    k2=(2*pi/lambda)^2;
    reffe = ReferenceFE(lagrangian,Float64,order);
    V = TestFESpace(model,reffe,vector_type=Vector{Float64});
    U = V;
    degree = 2*order;
    Ω = Triangulation(model);
    dΩ = Measure(Ω,degree);
    a(u,v) = ∫( -∇(v)⋅∇(u) + k2*(v⋅(eps*u))  )*dΩ;
    b(u,v) = ∫(v⋅u  )dΩ;
    A_task=assemble_matrix(a,U,V);
    B_task=assemble_matrix(b,U,V);
    A=fetch(A_task);
    B=fetch(B_task);
    if (solver==:LU)
        tmp1,tmp2=eigs_LU(A,B,sigma=approx_neff^2*k2,nev=mmax,tol=tol);
    elseif (solver==:MUMPS)
        tmp1,tmp2=eigs_MUMPS(A,B,sigma=approx_neff^2*k2,nev=mmax,tol=tol);
    else
        throw(DomainError(solver, "solver must be :LU or :MUMPS"));
    end
    neff=sqrt.(tmp1/k2);
    sol=Vector{ScalarModeFEM}(undef,0);
    if (field)
        E=[FEFunction(U,tmp2[:,i]) for i in axes(tmp2,2)];
        for i=1:length(neff)
            name=string("Mode ",string(i));
            push!(sol,ScalarModeFEM(name,neff[i],lambda,reffe,U,Ω,dΩ,E[i]));
        end
    else
        for i=1:length(neff)
            name=string("Mode ",string(i));
            push!(sol,ScalarModeFEM(name,neff[i],lambda));
        end
    end
    return sol;
end

function computeDz(Et,Ez,epsi::tensor3)
    ep(x)=VectorValue(epsi.zx(x),epsi.zy(x));
    return eps0*(ep⋅Et+epsi.zz*Ez)
end

function computeDt(Et,Ez,epsi::tensor3)
    eps1(x)=TensorValue(epsi.xx(x),epsi.yx(x),epsi.xy(x),epsi.yy(x));
    eps2(x)=VectorValue(epsi.xz(x),epsi.yz(x));
    return eps0*(eps1⋅Et+eps2*Ez)
end

function computeD(Et,Ez,epsi::tensor3)
    return computeDt(Et,Ez,epsi),computeDz(Et,Ez,epsi);
end

function computeB(Et,Ez,lambda::Real,neff)
    omega=2*pi*c/lambda;
    beta=2*pi*neff/lambda;
    Bz=curl(Et)/(im*omega);
    T1=TensorValue(0,-1,1,0);
    T2=TensorValue(0,im*beta,-im*beta,0);
    Bt=(T2⋅Et+T1⋅(∇(Ez)))*1.0/(im*omega);
    return Bt,Bz;
end

function computeB2(Et,Ez,lambda::Real,neff)
    omega=2*pi*c/lambda;
    beta=2*pi*neff/lambda;
    Bz=curl(Et)/(im*omega);
    T1=TensorValue(0,-im*2.0*pi*neff/lambda,im*2.0*pi*neff/lambda,0);
    T2=TensorValue(0,im*beta,-im*beta,0);
    Bt=(T2⋅Et+T1⋅(∇(Ez)))*1.0/(im*omega);
    return Bt,Bz;
end

function computeH(Bt,Bz,invmu::tensor3)
    invmu1(x)=TensorValue(invmu.xx(x),invmu.yx(x),invmu.xy(x),invmu.yy(x));
    invmu2(x)=VectorValue(invmu.xz(x),invmu.yz(x));
    invmu3(x)=VectorValue(invmu.zx(x),invmu.zy(x));
    return (invmu1⋅Bt+invmu2*Bz)*1/mu0,(invmu3⋅Bt+invmu.zz*Bz)*1/mu0
end

function FEM2(lambda::Real,mmax::Int64,eps::Function,model::DiscreteModel,approx_neff::Real;order::Int64=1,field::Bool=false,solver::Symbol=:LU,tol::Float64=0.0)
    k2=(2*pi/lambda)^2;
    reffe1 = ReferenceFE(nedelec,order)
    reffe2 = ReferenceFE(lagrangian,Float64,2*order)
    Ω = Triangulation(model);
    V1 = TestFESpace(Ω,reffe1,conformity=:Hcurl)
    V2 = TestFESpace(Ω,reffe2,conformity=:H1)
    U1 = V1;
    U2 = V2;
    V = MultiFieldFESpace([V1, V2])
    U = MultiFieldFESpace([U1, U2])
    degree = 4*order;
    dΩ = Measure(Ω,degree);
    #a((u1,u2),(v1,v2))=∫( -curl(v1)⋅curl(u1)*1.0/eps )*dΩ + ∫( (∇(u2))⋅(v1)*1.0/eps )*dΩ+∫( k2*(v1⋅u1+v2*u2) )*dΩ+ ∫( -(∇(u2))⋅(∇(v2))*1.0/eps )*dΩ
    #b((u1,u2),(v1,v2)) = ∫( (v1⋅u1)*1.0/eps )*dΩ+∫( -(∇(v2))⋅(u1)*1.0/eps )*dΩ
    a((Et1,Ez1),(Et2,Ez2))=∫( -curl(Et2)⋅curl(Et1) + ∇(Ez2)⋅(Et1) + k2*(Et2⋅(eps*Et1)+Ez2*(eps*Ez1)) -∇(Ez1)⋅∇(Ez2) )*dΩ
    b((Et1,Ez1),(Et2,Ez2)) = ∫( (Et2⋅Et1)-∇(Ez1)⋅(Et2) )*dΩ
    A_task=assemble_matrix(a,U,V);
    B_task=assemble_matrix(b,U,V);
    A=fetch(A_task);
    B=fetch(B_task);
    if (solver==:LU)
        tmp1,tmp2=eigs_LU(A,B,sigma=approx_neff^2*k2,nev=mmax,tol=tol);
    elseif (solver==:MUMPS)
        tmp1,tmp2=eigs_MUMPS(A,B,sigma=approx_neff^2*k2,nev=mmax,tol=tol);
    else
        throw(DomainError(solver, "solver must be :LU or :MUMPS"));
    end
    neff=sqrt.(tmp1/k2);
    sol=Vector{VectorModeFEM}(undef,0);
    if (field)
        for i=1:length(neff)
            name=string("Mode ",string(i));
            Et,Ez=FEFunction(U,tmp2[:,i])
            ux=VectorValue(1.0,0.0);
            uy=VectorValue(0.0,1.0);
            Ex=Et⋅ux;
            Ey=Et⋅uy;
            
            Bt,Bz=computeB2(Et,Ez,lambda,neff[i]);
            Ht,Hz=computeH(Bt,Bz,tensor3(x->1));
            Hx=Ht⋅ux;
            Hy=Ht⋅uy;
            push!(sol,VectorModeFEM(name,neff[i],lambda,reffe1,reffe2,U1,U2,Ω,dΩ,Ex,Ey,Ez*im*2.0*pi*neff[i]/lambda,Hx,Hy,Hz));
        end
    else
        for i=1:length(neff)
            name=string("Mode ",string(i));
            push!(sol,VectorModeFEM(name,neff[i],lambda));
        end
    end
    return sol
end

"""
    FEM(lambda::Real,mmax::Int64,eps::tensor3,mu::tensor3,model::DiscreteModel,approx_neff::Real;order::Int64=1,field::Bool=false,solver::Symbol=:LU,tol::Float64=0.0)

Returns a vector of `VectorModeFEM`.
The fiber is anisotropic and described with its relative permittivity tensor and its relative permeability tensor.

- lambda: wavelength
- mmax: maximal number of modes (useful if the fiber is very multimode)
- eps: functions of the tuple (x,y) that describes the relative permittivity tensor profile
- mu: functions of the tuple (x,y) that describes the relative permeability tensor profile
- model: DiscreteModel generated with `GridapGmsh.jl`
- approx_neff: effective index around which the modes will be computed
- order: order of the FEM
- field: boolean that indicates if fields must be saved
- solver: can be :LU or :MUMPS
- tol: tolerance for the eigenvalue solver (see documention of ArnoldiMethod.jl)
"""
function FEM(lambda::Real,mmax::Int64,eps::tensor3,mu::tensor3,model::DiscreteModel,approx_neff::Real;order::Int64=1,field::Bool=false,solver::Symbol=:LU,tol::Float64=0.0)
    k=2*pi/lambda;
    reffe1 = ReferenceFE(nedelec,order)
    reffe2 = ReferenceFE(lagrangian,Float64,2*order)
    Ω = Triangulation(model);
    V1 = TestFESpace(Ω,reffe1,conformity=:Hcurl,vector_type=Vector{ComplexF64})
    V2 = TestFESpace(Ω,reffe2,conformity=:H1,vector_type=Vector{ComplexF64}) 
    U1 = V1;
    U2 = V2;
    V = MultiFieldFESpace([V1, V2])
    U = MultiFieldFESpace([U1, U2])
    degree = 4*order;
    dΩ = Measure(Ω,degree);
    invmu=inverse(mu);
    tensor_C(x)=TensorValue(invmu.yy(x),-invmu.yx(x),-invmu.xy(x),invmu.xx(x))
    tensor_B(x)=TensorValue(invmu.yy(x),-invmu.xy(x),-invmu.yx(x),invmu.xx(x))
    vector_B1(x)=VectorValue(-invmu.zy(x),invmu.zx(x))
    vector_B2(x)=VectorValue(invmu.yz(x),-invmu.xz(x))
    
    a((Et1,Ez1),(Et2,Ez2))=∫( -curl(Et2)*(vector_B1⋅(∇(Ez1))+invmu.zz*curl(Et1))-curl(Et1)*(vector_B2⋅(∇(Ez2)))-(∇(Ez1))⋅(tensor_C⋅(∇(Ez2))) + k^2/eps0*(computeDt(Et1,Ez1,eps)⋅Et2+computeDz(Et1,Ez1,eps)*Ez2) )*dΩ
    b((Et1,Ez1),(Et2,Ez2))=∫( Et1⋅(tensor_C⋅(∇(Ez2)))-Et2⋅(tensor_B⋅(∇(Ez1)))+(vector_B1⋅Et1)*curl(Et2)+(vector_B2⋅Et2)*curl(Et1) )*dΩ
    c((Et1,Ez1),(Et2,Ez2))=∫( -(Et1⋅(tensor_C⋅Et2)) )*dΩ

    A_task=assemble_matrix(a,U,V);
    B_task=assemble_matrix(b,U,V);
    C_task=assemble_matrix(c,U,V);
    A=fetch(A_task);
    B=fetch(B_task);
    C=fetch(C_task);
    E,F=get_companion(A,im.*B,C);
    if (solver==:LU)
        tmp1,tmp2=eigs_LU(F,E,sigma=approx_neff*k,nev=mmax,tol=tol);
    elseif (solver==:MUMPS)
        tmp1,tmp2=eigs_MUMPS(F,E,sigma=approx_neff*k,nev=mmax,tol=tol);
    else
        throw(DomainError(solver, "solver must be :LU or :MUMPS"));
    end
    neff=tmp1/k;
    sol=Vector{VectorModeFEM}(undef,0);
    if (field)
        for i=1:length(neff)
            name=string("Mode ",string(i));
            Et,Ez=FEFunction(U,tmp2[:,i])
            ux=VectorValue(1.0,0.0);
            uy=VectorValue(0.0,1.0);
            Ex=Et⋅ux;
            Ey=Et⋅uy;
            Bt,Bz=computeB(Et,Ez,lambda,neff[i]);
            Ht,Hz=computeH(Bt,Bz,invmu);
            Hx=Ht⋅ux;
            Hy=Ht⋅uy;
            push!(sol,VectorModeFEM(name,neff[i],lambda,reffe1,reffe2,U1,U2,Ω,dΩ,Ex,Ey,Ez,Hx,Hy,Hz));
        end
    else
        for i=1:length(neff)
            name=string("Mode ",string(i));
            push!(sol,VectorModeFEM(name,neff[i],lambda));
        end
    end
    return sol
end
