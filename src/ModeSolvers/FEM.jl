"""
    FEM2D(lambda::Real,eps_fonc::Function,model::DiscreteModel;approx_neff::Real=0,neigs::Int64=1,order::Int64=2,field::Bool=false,solver::Symbol=:LU,tol::Float64=0.0,verbose::Bool=false,dPML::Real=0,alphaPML::Real=10,boundary_tag::String="",type::Symbol=:Scalar)

Returns a vector of `Mode{ScalarFieldFEM2D}` if type=:Scalar or a vector of `Mode{VectorFieldFEM2D}` if type=:Vector.  
The fiber is assumed to be isotropic and is described with its relative permittivity.

- lambda: wavelength
- eps_fonc: function of the tuple (x,y) that describes the relative permittivity profile i.e. eps(x)=n(x)^2
- model: DiscreteModel generated with `GridapGmsh.jl`
- neigs: number of modes to compute
- approx_neff: effective index around which the modes will be computed
- order: order of the FEM
- field: boolean that indicates if fields must be saved
- solver: can be :LU or :MUMPS
- tol: tolerance for the eigenvalue solver (see documention of ArnoldiMethod.jl)
- verbose: boolean that enables some outputs
- dPML: Thickness of the PML
- alphaPML: attenuation coefficient of the PML
- boundary_tag: tag of the boundary in model. If "", then the function automatically detects the boundary.
- type: :Scalar or :Vector
"""
function FEM2D(lambda::Real,eps_fonc::Function,model::DiscreteModel;approx_neff::Real=0,neigs::Int64=1,order::Int64=2,field::Bool=false,solver::Symbol=:LU,tol::Float64=0.0,verbose::Bool=false,dPML::Real=0,alphaPML::Real=10,boundary_tag::String="",type::Symbol=:Scalar)
    if (lambda<=0)
        throw(DomainError(lamdba, "The wavelength lambda must be strictly positive"));
    end
    if (first(methods(eps_fonc)).nargs!=2)
        throw(DomainError(eps_fonc, "The permittivity function must have 1 argument"));
    end
    if (num_dims(model)!=2)
        throw(DomainError(model, "The model must be 2D"));
    end
    if (neigs<=0)
        throw(DomainError(neigs, "The number of modes must be strictly positive"));
    end
    if (order<1)
        throw(DomainError(order, "The order must be at least 1"));
    end
    if (!(type in [:Scalar,:Vector]))
        throw(DomainError(solver, "solver must be :Scalar or :Vector"));
    end
    if (!(solver in [:LU,:MUMPS]))
        throw(DomainError(solver, "solver must be :LU or :MUMPS"));
    end
    if (tol<0)
        throw(DomainError(tol, "The tolerance must be positive or null"));
    end
    if (type==:Scalar)
        FEM2D_scalar(lambda,eps_fonc,model;approx_neff=approx_neff,neigs=neigs,order=order,field=field,solver=solver,tol=tol,verbose=verbose,dPML=dPML,alphaPML=alphaPML,boundary_tag=boundary_tag)
    else
        if (dPML<=0)
            FEM2D_vector_guided(lambda,eps_fonc,model;approx_neff=approx_neff,neigs=neigs,order=order,field=field,solver=solver,tol=tol,verbose=verbose,boundary_tag=boundary_tag)
        else
            FEM2D_vector_leaky(lambda,eps_fonc,model;approx_neff=approx_neff,neigs=neigs,order=order,field=field,solver=solver,tol=tol,verbose=verbose,dPML=dPML,alphaPML=alphaPML,boundary_tag=boundary_tag)
        end
    end
end

"""
    FEM1D(lambda::Real,nu::Int64,eps_fonc::Function,model::DiscreteModel;approx_neff::Real=0,neigs::Int64=1,order::Int64=2,field::Bool=false,solver::Symbol=:LU,tol::Real=0.0,verbose::Bool=false,dPML::Real=0,alphaPML::Real=10)

Returns a vector of `Mode{ScalarFieldFEM1D}`.  
The fiber is assumed to be isotropic with a cylindrical symmetry and is described with its relative permittivity.

- lambda: wavelength
- nu: azimuthal number
- eps_fonc: function of r that describes the relative permittivity profile i.e. eps(r)=n(r)^2
- model: DiscreteModel generated with `GridapGmsh.jl`
- neigs: number of modes to compute
- approx_neff: effective index around which the modes will be computed
- order: order of the FEM
- field: boolean that indicates if fields must be saved
- solver: can be :LU or :MUMPS
- tol: tolerance for the eigenvalue solver (see documention of ArnoldiMethod.jl)
- verbose: boolean that enables some outputs
- dPML: Thickness of the PML
- alphaPML: attenuation coefficient of the PML
"""
function FEM1D(lambda::Real,nu::Int64,eps_fonc::Function,model::DiscreteModel;approx_neff::Real=0,neigs::Int64=1,order::Int64=2,field::Bool=false,solver::Symbol=:LU,tol::Real=0.0,verbose::Bool=false,dPML::Real=0,alphaPML::Real=10)
    if (lambda<=0)
        throw(DomainError(lamdba, "The wavelength lambda must be strictly positive"));
    end
    if (nu<0)
        throw(DomainError(nu, "The azymuthal number must be positive or null"));
    end
    if (neigs<=0)
        throw(DomainError(neigs, "The number of modes must be strictly positive"));
    end
    if (first(methods(eps_fonc)).nargs!=2)
        throw(DomainError(eps_fonc, "The permittivity function must have 1 argument"));
    end
    if (num_dims(model)!=1)
        throw(DomainError(model, "The model must be 1D"));
    end
    if (approx_neff<0)
        throw(DomainError(approx_neff, "The approximative effective index must be positive or null"));
    end
    if (order<1)
        throw(DomainError(order, "The order must be at least 1"));
    end
    if (!(solver in [:LU,:MUMPS]))
        throw(DomainError(solver, "solver must be :LU or :MUMPS"));
    end
    if (tol<0)
        throw(DomainError(tol, "The tolerance must be positive or null"));
    end

    boundary=Gridap.Geometry.compute_isboundary_face(model.grid_topology,0);
    if (sum(boundary)!=2)
        throw(DomainError(model, "The number of boundaries is different from 2, check the model"));
    end
    coord=Gridap.ReferenceFEs.get_node_coordinates(model.grid)
    r=[coord[i][1] for i in axes(coord,1)];
    (rmin,posmin)=findmin(r);
    (rmax,posmax)=findmax(r);
    if (verbose)
        println("Find rmin = ",rmin)
        println("Find rmax = ",rmax)
    end
    if (rmin!=0)
        throw(DomainError(model, "The geometry must begin at the origin"));
    end
    if (boundary[posmin]!=1 || boundary[posmax]!=1)
        throw(DomainError(model, "rmin and/or rmax do not correspond to boundaries"));
    end

    labels = get_face_labeling(model)
    ent_min=labels.d_to_dface_to_entity[1][posmin]
    ent_max=labels.d_to_dface_to_entity[1][posmax]
    find_tag=false;
    label_start="";
    label_end="";
    if ((length(findall(==(ent_min),labels.d_to_dface_to_entity[1]))==1) && (length(findall(==(ent_min),labels.d_to_dface_to_entity[2]))==0) && (length(findall(==(ent_max),labels.d_to_dface_to_entity[1]))==1) && (length(findall(==(ent_max),labels.d_to_dface_to_entity[2]))==0))
        list_tag=Gridap.Geometry.get_tag_name(labels)
        tag_entities=get_tag_entities(labels)
        for j in axes(tag_entities,1)
            if tag_entities[j]==[ent_min]
                label_start=list_tag[j]
            end
            if tag_entities[j]==[ent_max]
                label_end=list_tag[j]
            end
        end
        find_tag=(!isempty(label_start) && !isempty(label_end))
    end

    if (find_tag)
        model2=model
        if (verbose)
            println("Boundary tags detected : '",label_start,"' and '",label_end,"'.")
        end
    else
        model2=deepcopy(model)
        labels2 = get_face_labeling(model2)
        list_tag=Gridap.Geometry.get_tag_name(labels2)
        label_start="start";
        while (label_start in list_tag)
            label_start=label_start*"_";
        end
        label_end="end";
        while (label_end in list_tag)
            label_end=label_end*"_";
        end
        tag_start=maximum(get_face_entity(labels2))+1;
        tag_end=tag_start+1;
        Gridap.Geometry.add_tag!(labels2,label_start,[tag_start])
        Gridap.Geometry.add_tag!(labels2,label_end,[tag_end])
        labels2.d_to_dface_to_entity[1][posmin]=tag_start
        labels2.d_to_dface_to_entity[1][posmax]=tag_end
        if (verbose)
            println("Boundary tags '",label_start,"' and '",label_end,"' created.")
        end
    end

    nend=sqrt(eps_fonc(rmax));
    if (verbose)
        println("Cladding refractive index = ",nend)
    end
    if (approx_neff<=0)
        approx_neff=maximum(sqrt.(eps_fonc.(coord)))
        if (verbose)
            println("Set approx_neff to ",approx_neff)
        end
    end

    reffe = ReferenceFE(lagrangian,Float64,order);
    Ω = Triangulation(model2)
    dΩ = Measure(Ω,2*order)
    if (dPML<=0)
        if (nu==0)
            V = TestFESpace(model2,reffe,conformity=:H1,dirichlet_tags=[label_end]);
        else
            V = TestFESpace(model2,reffe,conformity=:H1,dirichlet_tags=[label_start,label_end]);
        end
    else
        if (nu==0)
            V = TestFESpace(model2,reffe,conformity=:H1,dirichlet_tags=[label_end],vector_type=Vector{ComplexF64});
        else
            V = TestFESpace(model2,reffe,conformity=:H1,dirichlet_tags=[label_start,label_end],vector_type=Vector{ComplexF64});
        end
    end
    U=TrialFESpace(V,0);
    k2=(2*pi/lambda)^2
    r2=x->x[1]^2;
    r=x->x[1];
    epsi=x->eps_fonc(x)*x[1]^2*k2-nu^2;
    ux=VectorValue(1.0,);
    if (dPML<=0)
        A_task=Threads.@spawn assemble_matrix((u,v)->∫( -(∇(v)⋅∇(u))*r2 - v*(ux⋅∇(u))*r + (v⋅(epsi*u))  )*dΩ,U,V)
        B_task=Threads.@spawn assemble_matrix((u,v)->∫((v⋅u)*r2  )dΩ,U,V)
    else
        rPML=rmax-dPML;
        alpha=x->alphaPML*(x[1]>rPML)
        s=x->(1.0+im*alpha(x)*(x[1]-rPML)^2/dPML^2);
        A_task=Threads.@spawn assemble_matrix((u,v)->∫( -(∇(v)⋅∇(u))*r2/s - v*(ux⋅∇(u))*r + (v⋅(epsi*u))*s  )*dΩ,U,V)
        B_task=Threads.@spawn assemble_matrix((u,v)->∫((v⋅u)*r2*s  )dΩ,U,V)
    end
    A=fetch(A_task);
    B=fetch(B_task);
    if (verbose)
        println("Matrices of dimension ",size(A), " created.")
    end
    if (solver==:LU)
        tmp1,tmp2=eigs_LU(A,B,sigma=approx_neff^2*k2,nev=neigs,tol=Float64(tol),verbose=verbose);
    else
        tmp1,tmp2=eigs_MUMPS(A,B,sigma=approx_neff^2*k2,nev=neigs,tol=Float64(tol),verbose=verbose);
    end
    neff=sqrt.(tmp1/k2);
    if (verbose)
        println("Found ",length(neff)," eigenvalues.")
        if (dPML<=0)
            println(sum(neff.>nend)," eigenvalues correspond to guided modes (neff>n(cladding)).")
        end
    end
    if (dPML<=0)
        nb_mode=1;
        if (field)
            sol=Vector{Mode{ScalarFieldFEM1D}}(undef,0);
            E=[FEFunction(U,tmp2[:,i]) for i in axes(tmp2,2)];
            for i in axes(neff,1)
                if (neff[i]>nend) #Guided mode only
                    name=string("Mode LP n°",string(nb_mode));
                    push!(sol,Mode(name,neff[i],lambda,ScalarFieldFEM1D(nu,dΩ,E[i])));
                    nb_mode=nb_mode+1;
                end
            end
        else
            sol=Vector{Mode{Nothing}}(undef,0);
            for i in axes(neff,1)
                if (neff[i]>nend)
                    name=string("Mode LP n°",string(nb_mode));
                    push!(sol,Mode(name,neff[i],lambda));
                    nb_mode=nb_mode+1;
                end
            end
        end
    else
        if (field)
            sol=Vector{Mode{ScalarFieldFEM1D}}(undef,0);
            E=[FEFunction(U,tmp2[:,i]) for i in axes(tmp2,2)];
            for i in axes(neff,1)
                name=string("Mode LP n°",string(i));
                push!(sol,Mode(name,neff[i],lambda,ScalarFieldFEM1D(nu,dΩ,E[i])));
            end
        else
            sol=Vector{Mode{Nothing}}(undef,0);
            for i in axes(neff,1)
                name=string("Mode LP n°",string(i));
                push!(sol,Mode(name,neff[i],lambda));
            end
        end
    end
    return sol;
end

function create_model_with_boundary(model::DiscreteModel;verbose::Bool=false)
    #Get boundaries in 0 and 1 dimension
    boundary_1D=Gridap.Geometry.compute_isboundary_face(model.grid_topology,1);
    boundary_0D=Gridap.Geometry.compute_isboundary_face(model.grid_topology,0);
    pos_boundary_1D=findall(boundary_1D);
    pos_boundary_0D=findall(boundary_0D);
    labels = get_face_labeling(model)
    
    #Check if a boundary tag exists
    ent_0D=unique(labels.d_to_dface_to_entity[1][pos_boundary_0D])
    ent_1D=unique(labels.d_to_dface_to_entity[2][pos_boundary_1D])
    ent=sort(unique([ent_0D;ent_1D]))
    label_boundary=""
    if ((length(findall(in(ent_0D),labels.d_to_dface_to_entity[1]))==length(pos_boundary_0D)) && (length(findall(in(ent_0D),labels.d_to_dface_to_entity[2]))==length(pos_boundary_1D)) && (length(findall(in(ent_0D),labels.d_to_dface_to_entity[3]))==0))
        list_tag=Gridap.Geometry.get_tag_name(labels)
        tag_entities=get_tag_entities(labels)
        for j in axes(tag_entities,1)
            if sort(tag_entities[j])==ent
                label_boundary=list_tag[j]
            end
        end
    end
    if (isempty(label_boundary))
        model2=deepcopy(model)
        #add a specific label for boundaries
        labels2 = get_face_labeling(model2)
        #find a new name for boundary tag
        list_tag=Gridap.Geometry.get_tag_name(labels2)
        label_boundary="boundary";
        while (label_boundary in list_tag)
            label_boundary=label_boundary*"_";
        end
        #find the number for the list of the boundary entities
        tag_entity=maximum(get_face_entity(labels2))+1;
        Gridap.Geometry.add_tag!(labels2,label_boundary,[tag_entity])
        labels2.d_to_dface_to_entity[1][pos_boundary_0D].=tag_entity
        labels2.d_to_dface_to_entity[2][pos_boundary_1D].=tag_entity
        if (verbose)
            println("Boundary tag '",label_boundary,"' created.")
        end
        return label_boundary,model2
    else
        if (verbose)
            println("Boundary tag '",label_boundary,"' detected.")
        end
        return label_boundary,model
    end
end

function detect_boundaries(model::DiscreteModel;verbose::Bool=false)
    boundary_0D=Gridap.Geometry.compute_isboundary_face(model.grid_topology,0);
    pos_boundary_0D=findall(boundary_0D)
    coord=Gridap.ReferenceFEs.get_node_coordinates(model.grid)
    cfn=compute_face_nodes(model,0)
    x_boundary=[coord[cfn[pos_boundary_0D[i]][1]][1] for i in axes(pos_boundary_0D,1)]
    y_boundary=[coord[cfn[pos_boundary_0D[i]][1]][2] for i in axes(pos_boundary_0D,1)]
    R_boundary=hypot.(x_boundary,y_boundary)
    R=sum(R_boundary)/length(R_boundary)
    if (sum(isapprox.(R_boundary,R))==length(R_boundary))
        if (verbose)
            println("Circular boundary detected. Radius = ",R)
        end
        boundaries=R;
    else
        x_boundary_min,x_boundary_max=extrema(x_boundary)
        y_boundary_min,y_boundary_max=extrema(y_boundary)
        ok=(isapprox.(x_boundary,x_boundary_max) .| isapprox.(x_boundary,x_boundary_min) .| isapprox.(y_boundary,y_boundary_min) .| isapprox.(y_boundary,y_boundary_max))
        if (sum(ok)==length(x_boundary))
            if (verbose)
                println("Rectangular boundary detected: xmin = ",x_boundary_min," ; xmax = ",x_boundary_max," ; ymin = ",y_boundary_min," ; ymax = ",y_boundary_max);
            end
            boundaries=[x_boundary_min,x_boundary_max,y_boundary_min,y_boundary_max];
        else
            throw(DomainError(model, "When using a PML, the model boundary must be circular or rectangular"));
        end
    end
    return boundaries
end

function FEM2D_scalar(lambda::Real,eps_fonc::Function,model::DiscreteModel;approx_neff::Real=0,neigs::Int64=1,order::Int64=2,field::Bool=false,solver::Symbol=:LU,tol::Real=0.0,verbose::Bool=false,dPML::Real=0,alphaPML::Real=10,boundary_tag::String="")
    if (lambda<=0)
        throw(DomainError(lamdba, "The wavelength lambda must be strictly positive"));
    end
    if (neigs<=0)
        throw(DomainError(neigs, "The number of modes must be strictly positive"));
    end
    if (first(methods(eps_fonc)).nargs!=2)
        throw(DomainError(eps_fonc, "The permittivity function must have 1 argument"));
    end
    if (num_dims(model)!=2)
        throw(DomainError(model, "The model must be 2D"));
    end
    if (approx_neff<0)
        throw(DomainError(approx_neff, "The approximative effective index must be positive or null"));
    end
    if (order<1)
        throw(DomainError(order, "The order must be at least 1"));
    end
    if (!(solver in [:LU,:MUMPS]))
        throw(DomainError(solver, "solver must be :LU or :MUMPS"));
    end
    if (tol<0)
        throw(DomainError(tol, "The tolerance must be positive or null"));
    end

    if (isempty(boundary_tag))
        boundary_tag,model2=create_model_with_boundary(model,verbose=verbose)
    else
        model2=model
    end

    if (dPML>0)
        boundaries=detect_boundaries(model2,verbose=verbose)
        if (length(boundaries)==1)
            circular=true
            R=boundaries
        else
            circular=false
            x_boundary_min=boundaries[1]
            x_boundary_max=boundaries[2]
            y_boundary_min=boundaries[3]
            y_boundary_max=boundaries[4]
        end
    end

    boundary_0D=Gridap.Geometry.compute_isboundary_face(model2.grid_topology,0);
    pos_boundary_0D=findall(boundary_0D)
    coord=Gridap.ReferenceFEs.get_node_coordinates(model2.grid)
    cfn=compute_face_nodes(model2,0)
    nend=maximum(sqrt.(eps_fonc.(coord[first.(cfn[pos_boundary_0D])])));
    if (verbose)
        println("Cladding refractive index = ",nend)
    end
    if (approx_neff<=0)
        approx_neff=maximum(sqrt.(eps_fonc.(coord)))
        if (verbose)
            println("Set approx_neff to ",approx_neff)
        end
    end


    k2=(2*pi/lambda)^2;
    reffe = ReferenceFE(lagrangian,Float64,order);
    if (dPML<=0)
        V = TestFESpace(model2,reffe,conformity=:H1,vector_type=Vector{Float64},dirichlet_tags=[boundary_tag])
    else
        V = TestFESpace(model2,reffe,conformity=:H1,vector_type=Vector{ComplexF64},dirichlet_tags=[boundary_tag])
    end
    U=TrialFESpace(V,0);
    degree = 2*order;
    Ω = Triangulation(model2);
    dΩ = Measure(Ω,degree);
    if (dPML<=0)
        A_task=Threads.@spawn assemble_matrix((u,v)->∫( -∇(v)⋅∇(u) + k2*(v⋅(eps_fonc*u))  )*dΩ,U,V);
    else
        if (circular)
            rPML=R-dPML
            r=x->hypot(x[1],x[2])
            phi=x->atan(x[2],x[1]);
            Sr=x->1.0+im*alphaPML*(r(x)-rPML)^2/(dPML)^2
            Sphi=x->1.0+im*alphaPML/3.0*(r(x)-rPML)^3/(dPML)^2/r(x);
            C=x->cos(phi(x));
            S=x->sin(phi(x));
            iSxx=x->(r(x)<rPML) ? ComplexF64(1.0) : C(x)^2/Sr(x)+S(x)^2/Sphi(x);
            iSxy=x->(r(x)<rPML) ? ComplexF64(0.0) : (1/Sr(x)-1/Sphi(x))*C(x)*S(x);
            iSyy=x->(r(x)<rPML) ? ComplexF64(1.0) : S(x)^2/Sr(x)+C(x)^2/Sphi(x);
            dSr=x->2*im*alphaPML*(r(x)-rPML)/(dPML)^2
            dSphi=x->-im*alphaPML/3.0*(r(x)-rPML)^3/(dPML)^2/r(x)^2+im*alphaPML*(r(x)-rPML)^2/(dPML)^2/r(x)
            diSxx_dphi=x->-2*C(x)*S(x)/Sr(x)+2*C(x)*S(x)/Sphi(x)
            diSxy_dphi=x->(1/Sr(x)-1/Sphi(x))*(C(x)^2-S(x)^2)
            diSyy_dphi=x->2*C(x)*S(x)/Sr(x)-2*C(x)*S(x)/Sphi(x)
            diSxx_dr=x->-C(x)^2*dSr(x)/Sr(x)^2-S(x)^2*dSphi(x)/Sphi(x)^2
            diSxy_dr=x->(-dSr(x)/Sr(x)^2+dSphi(x)/Sphi(x)^2)*C(x)*S(x)
            diSyy_dr=x->-S(x)^2*dSr(x)/Sr(x)^2-C(x)^2*dSphi(x)/Sphi(x)^2
            diSxx_dx=x->(r(x)<rPML) ? ComplexF64(0.0) : C(x)*diSxx_dr(x)-S(x)/r(x)*diSxx_dphi(x)
            diSxy_dx=x->(r(x)<rPML) ? ComplexF64(0.0) : C(x)*diSxy_dr(x)-S(x)/r(x)*diSxy_dphi(x)
            diSxy_dy=x->(r(x)<rPML) ? ComplexF64(0.0) : S(x)*diSxy_dr(x)+C(x)/r(x)*diSxy_dphi(x)
            diSyy_dy=x->(r(x)<rPML) ? ComplexF64(0.0) : S(x)*diSyy_dr(x)+C(x)/r(x)*diSyy_dphi(x)
            ux=VectorValue(1.0,0.0);
            uy=VectorValue(0.0,1.0);
            fx=x->iSxx(x)*diSxx_dx(x)+iSxy(x)*diSyy_dy(x)+iSxx(x)*diSxy_dy(x)+iSxy(x)*diSxy_dx(x)
            fy=x->iSxy(x)*diSxx_dx(x)+iSxy(x)*diSxy_dy(x)+iSyy(x)*diSyy_dy(x)+iSyy(x)*diSxy_dx(x)
            gxx=x->iSxx(x)*iSxx(x)+iSxy(x)*iSxy(x)
            gxy=x->iSxx(x)*iSxy(x)+iSxy(x)*iSyy(x)
            gyx=x->iSxy(x)*iSxx(x)+iSyy(x)*iSxy(x)
            gyy=x->iSxy(x)*iSxy(x)+iSyy(x)*iSyy(x)
            a1=(u,v)->-(∇(u)⋅ux)*(∇(v)⋅ux)*gxx - (∇(u)⋅ux)*v*fx
            a2=(u,v)->-(∇(u)⋅ux)*(∇(v)⋅uy)*gxy 
            a3=(u,v)->-(∇(u)⋅uy)*(∇(v)⋅ux)*gyx - (∇(u)⋅uy)*v*fy
            a4=(u,v)->-(∇(u)⋅uy)*(∇(v)⋅uy)*gyy 
            A_task=Threads.@spawn assemble_matrix((u,v)->∫( a1(u,v) + a2(u,v) + a3(u,v) + a4(u,v) + k2*(v⋅(eps_fonc*u))  )*dΩ,U,V)
        else
            xPML_max=x_boundary_max-dPML
            xPML_min=x_boundary_min+dPML
            yPML_max=y_boundary_max-dPML
            yPML_min=y_boundary_min+dPML
            alphax_max=x->alphaPML*(x[1]>xPML_max)
            alphax_min=x->alphaPML*(x[1]<xPML_min)
            alphay_max=x->alphaPML*(x[2]>yPML_max)
            alphay_min=x->alphaPML*(x[2]<yPML_min)
            sx=x->(1.0+im*alphax_max(x)*(x[1]-xPML_max)^2/dPML^2+im*alphax_min(x)*(x[1]-xPML_min)^2/dPML^2);
            sx_prime=x->(im*alphax_max(x)*2*(x[1]-xPML_max)/dPML^2+im*alphax_min(x)*2*(x[1]-xPML_min)/dPML^2)
            sy=x->(1.0+im*alphay_max(x)*(x[2]-yPML_max)^2/dPML^2+im*alphay_min(x)*(x[2]-yPML_min)^2/dPML^2);
            sy_prime=x->(im*alphay_max(x)*2*(x[2]-yPML_max)/dPML^2+im*alphay_min(x)*2*(x[2]-yPML_min)/dPML^2)
            sx2=x->sx(x)^2
            sy2=x->sy(x)^2
            sx3=x->sx_prime(x)/sx(x)^3
            sy3=x->sy_prime(x)/sy(x)^3
            ux=VectorValue(1.0,0.0);
            uy=VectorValue(0.0,1.0);
            A_task=Threads.@spawn assemble_matrix((u,v)->∫( -(∇(u)⋅ux)*((∇(v)⋅ux)/sx2-v*sx3) -(∇(u)⋅uy)*((∇(v)⋅uy)/sy2-v*sy3)  + k2*(v⋅(eps_fonc*u))  )*dΩ,U,V)
        end
    end
    B_task=Threads.@spawn assemble_matrix((u,v)->∫(v⋅u  )dΩ,U,V);
    A=fetch(A_task);
    B=fetch(B_task);
    if (verbose)
        println("Matrices of dimension ",size(A), " created.")
    end
    if (solver==:LU)
        tmp1,tmp2=eigs_LU(A,B,sigma=approx_neff^2*k2,nev=neigs,tol=Float64(tol),verbose=verbose);
    else
        tmp1,tmp2=eigs_MUMPS(A,B,sigma=approx_neff^2*k2,nev=neigs,tol=Float64(tol),verbose=verbose);
    end
    neff=sqrt.(tmp1/k2);
    if (verbose)
        println("Found ",length(neff)," eigenvalues.")
        if (dPML<=0)
            println(sum(neff.>nend)," eigenvalues correspond to guided modes (neff>n(cladding)).")
        end
    end
    if (dPML<=0)
        nb_mode=1
        if (field)
            sol=Vector{Mode{ScalarFieldFEM2D}}(undef,0);
            E=[FEFunction(U,tmp2[:,i]) for i in axes(tmp2,2)];
            for i in axes(neff,1)
                if (neff[i]>nend) #Guided mode only
                    name=string("Mode LP n°",string(nb_mode));
                    push!(sol,Mode(name,neff[i],lambda,ScalarFieldFEM2D(dΩ,E[i])));
                    nb_mode=nb_mode+1;
                end
            end
        else
            sol=Vector{Mode{Nothing}}(undef,0);
            for i in axes(neff,1)
                if (neff[i]>nend) #Guided mode only
                    name=string("Mode LP n°",string(nb_mode));
                    push!(sol,Mode(name,neff[i],lambda));
                    nb_mode=nb_mode+1;
                end
            end
        end
    else
        if (field)
            sol=Vector{Mode{ScalarFieldFEM2D}}(undef,0);
            E=[FEFunction(U,tmp2[:,i]) for i in axes(tmp2,2)];
            for i in axes(neff,1)
                name=string("Mode LP n°",string(i));
                push!(sol,Mode(name,neff[i],lambda,ScalarFieldFEM2D(dΩ,E[i])));
            end
        else
            sol=Vector{Mode{Nothing}}(undef,0);
            for i in axes(neff,1)
                name=string("Mode LP n°",string(i));
                push!(sol,Mode(name,neff[i],lambda));
            end
        end
    end
    return sol;
end

function computeDz(Et,Ez,epsi::tensor3)
    ep(x)=VectorValue(epsi.zx(x),epsi.zy(x));
    return eps0*(ep⋅Et+epsi.zz.f*Ez)
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
    return (invmu1⋅Bt+invmu2*Bz)*1/mu0,(invmu3⋅Bt+invmu.zz.f*Bz)*1/mu0
end


function FEM2D_vector_guided(lambda::Real,eps_fonc::Function,model::DiscreteModel;approx_neff::Real=0,neigs::Int64=1,order::Int64=2,field::Bool=false,solver::Symbol=:LU,tol::Real=0.0,verbose::Bool=false,boundary_tag::String="")
    if (lambda<=0)
        throw(DomainError(lamdba, "The wavelength lambda must be strictly positive"));
    end
    if (neigs<=0)
        throw(DomainError(nb, "The number of modes must be strictly positive"));
    end
    if (first(methods(eps_fonc)).nargs!=2)
        throw(DomainError(eps_fonc, "The permittivity function must have 1 argument"));
    end
    if (num_dims(model)!=2)
        throw(DomainError(model, "The model must be 2D"));
    end
    if (approx_neff<0)
        throw(DomainError(approx_neff, "The approximative effective index must be positive or null"));
    end
    if (order<1)
        throw(DomainError(order, "The order must be at least 1"));
    end
    if (!(solver in [:LU,:MUMPS]))
        throw(DomainError(solver, "solver must be :LU or :MUMPS"));
    end
    if (tol<0)
        throw(DomainError(tol, "The tolerance must be positive or null"));
    end

    if (isempty(boundary_tag))
        boundary_tag,model2=create_model_with_boundary(model,verbose=verbose)
    else
        model2=model
    end
        
    boundary_0D=Gridap.Geometry.compute_isboundary_face(model2.grid_topology,0);
    pos_boundary_0D=findall(boundary_0D);
    coord=Gridap.ReferenceFEs.get_node_coordinates(model2.grid)
    cfn=compute_face_nodes(model2,0)
    nend=maximum(sqrt.(eps_fonc.(coord[first.(cfn[pos_boundary_0D])])));
    if (verbose)
        println("Cladding refractive index = ",nend)
    end
    if (approx_neff<=0)
        approx_neff=maximum(sqrt.(eps_fonc.(coord)))
        if (verbose)
            println("Set approx_neff to ",approx_neff)
        end
    end   


    k2=(2*pi/lambda)^2;
    reffe1 = ReferenceFE(nedelec,order)
    reffe2 = ReferenceFE(lagrangian,Float64,order+1)
    Ω = Triangulation(model2);
    V1 = TestFESpace(Ω,reffe1,conformity=:Hcurl,dirichlet_tags=[boundary_tag])
    V2 = TestFESpace(Ω,reffe2,conformity=:H1,dirichlet_tags=[boundary_tag])
    U1 = V1;
    U2 = V2;
    V = MultiFieldFESpace([V1, V2])
    U = MultiFieldFESpace([U1, U2])
    degree = 2*order+2;
    dΩ = Measure(Ω,degree);
    a((Et1,Ez1),(Et2,Ez2))=∫( -curl(Et2)⋅curl(Et1) + ∇(Ez2)⋅(Et1) + k2*(Et2⋅(eps_fonc*Et1)+Ez2*(eps_fonc*Ez1)) -∇(Ez1)⋅∇(Ez2) )*dΩ
    b((Et1,Ez1),(Et2,Ez2)) = ∫( (Et2⋅Et1)-∇(Ez1)⋅(Et2) )*dΩ
    A_task=Threads.@spawn assemble_matrix(a,U,V);
    B_task=Threads.@spawn assemble_matrix(b,U,V);
    A=fetch(A_task);
    B=fetch(B_task);
    if (verbose)
        println("Matrices of dimension ",size(A), " created.")
    end
    if (solver==:LU)
        tmp1,tmp2=eigs_LU(A,B,sigma=approx_neff^2*k2,nev=neigs,tol=Float64(tol),verbose=verbose);
    else
        tmp1,tmp2=eigs_MUMPS(A,B,sigma=approx_neff^2*k2,nev=neigs,tol=Float64(tol),verbose=verbose);
    end
    neff=sqrt.(tmp1/k2);
    if (verbose)
        println("Found ",length(neff)," eigenvalues.")
        println(sum(neff.>nend)," eigenvalues correspond to guided modes (neff>n(cladding)).")
    end
    nb_mode=1
    if (field)
        sol=Vector{Mode{VectorFieldFEM2D}}(undef,0);
        for i in axes(neff,1)
            if (real(neff[i])>nend) #Guided mode only
                name=string("Mode ",string(nb_mode));
                Et,Ez=FEFunction(U,tmp2[:,i])
                ux=VectorValue(1.0,0.0);
                uy=VectorValue(0.0,1.0);
                Ex=Et⋅ux;
                Ey=Et⋅uy;
                Bt,Bz=computeB2(Et,Ez,lambda,neff[i]);
                Ht,Hz=computeH(Bt,Bz,tensor3(1));
                Hx=Ht⋅ux;
                Hy=Ht⋅uy;
                push!(sol,Mode(name,neff[i],lambda,VectorFieldFEM2D(dΩ,Ex,Ey,Ez*im*2.0*pi*neff[i]/lambda,Hx,Hy,Hz)));
                nb_mode=nb_mode+1
            end
        end
    else
        sol=Vector{Mode{Nothing}}(undef,0);
        for i in axes(neff,1)
            if (neff[i]>nend) #Guided mode only
                name=string("Mode ",string(nb_mode));
                push!(sol,Mode(name,neff[i],lambda));
                nb_mode=nb_mode+1
            end
        end
    end
    return sol;
end

function FEM2D_vector_leaky(lambda::Real,eps_fonc::Function,model::DiscreteModel;approx_neff::Real=0,neigs::Int64=1,order::Int64=2,field::Bool=false,solver::Symbol=:LU,tol::Real=0.0,verbose::Bool=false,dPML::Real=0,alphaPML::Real=10,boundary_tag::String="")
    if (lambda<=0)
        throw(DomainError(lamdba, "The wavelength lambda must be strictly positive"));
    end
    if (neigs<=0)
        throw(DomainError(nb, "The number of modes must be strictly positive"));
    end
    if (first(methods(eps_fonc)).nargs!=2)
        throw(DomainError(eps_fonc, "The permittivity function must have 1 argument"));
    end
    if (num_dims(model)!=2)
        throw(DomainError(model, "The model must be 2D"));
    end
    if (approx_neff<0)
        throw(DomainError(approx_neff, "The approximative effective index must be positive or null"));
    end
    if (order<1)
        throw(DomainError(order, "The order must be at least 1"));
    end
    if (!(solver in [:LU,:MUMPS]))
        throw(DomainError(solver, "solver must be :LU or :MUMPS"));
    end
    if (tol<0)
        throw(DomainError(tol, "The tolerance must be positive or null"));
    end

    
    if (isempty(boundary_tag))
        boundary_tag,model2=create_model_with_boundary(model,verbose=verbose)
    else
        model2=model
    end
        
    if (dPML>0)
        boundaries=detect_boundaries(model2,verbose=verbose)
        if (length(boundaries)==1)
            circular=true
            R=boundaries
        else
            circular=false
            x_boundary_min=boundaries[1]
            x_boundary_max=boundaries[2]
            y_boundary_min=boundaries[3]
            y_boundary_max=boundaries[4]
        end
    end

    coord=Gridap.ReferenceFEs.get_node_coordinates(model2.grid)
    if (verbose)
        println("Cladding refractive index = ",nend)
    end
    if (approx_neff<=0)
        approx_neff=maximum(sqrt.(eps_fonc.(coord)))
        if (verbose)
            println("Set approx_neff to ",approx_neff)
        end
    end

    if (circular)
        epsilon_tensor=add_cylindrical_PML(eps_fonc,R-dPML,dPML,alphaPML);
        mu_tensor=add_cylindrical_PML(x->1.0,R-dPML,dPML,alphaPML);
    else
        epsilon_tensor=add_rectangular_PML(eps_fonc,x_boundary_min+dPML,x_boundary_max-dPML,dPML,y_boundary_min+dPML,y_boundary_max-dPML,dPML,alphaPML);
        mu_tensor=add_rectangular_PML(x->1.0,x_boundary_min+dPML,x_boundary_max-dPML,dPML,y_boundary_min+dPML,y_boundary_max-dPML,dPML,alphaPML);
    end
    FEM2D_anisotropic(lambda,epsilon_tensor,mu_tensor,model2,approx_neff=approx_neff,neigs=neigs,order=order,field=field,solver=solver,tol=tol,verbose=verbose,boundary_tag=boundary_tag)
end
"""
    FEM2D_anisotropic(lambda::Real,epsilon::tensor3,mu::tensor3,model::DiscreteModel;approx_neff::Real=0,neigs::Int64=1,order::Int64=1,field::Bool=false,solver::Symbol=:LU,tol::Float64=0.0,verbose::Bool=false,boundary_tag::String="")

Returns a vector of `Mode{VectorFieldFEM2D}`.
The fiber is anisotropic and described with its relative permittivity tensor and its relative permeability tensor. The PML is assumed to be already included in the tensors epsilon and mu.

- lambda: wavelength
- epsilon: `tensor3` with functions of (x,y) that describes the relative permittivity tensor profile
- mu: `tensor3` with functions of (x,y) that describes the relative permeability tensor profile
- model: DiscreteModel generated with `GridapGmsh.jl`
- approx_neff: effective index around which the modes will be computed
- neigs: number of modes to calculate
- order: order of the FEM
- field: boolean that indicates if fields must be saved
- solver: can be :LU or :MUMPS
- tol: tolerance for the eigenvalue solver (see documention of ArnoldiMethod.jl)
- verbose: boolean that enables some outputs
- boundary_tag: tag of the boundary in model.

"""
function FEM2D_anisotropic(lambda::Real,epsilon::tensor3,mu::tensor3,model::DiscreteModel;approx_neff::Real=0,neigs::Int64=1,order::Int64=1,field::Bool=false,solver::Symbol=:LU,tol::Float64=0.0,verbose::Bool=false,boundary_tag::String="")
    k=2*pi/lambda;
    reffe1 = ReferenceFE(nedelec,order)
    reffe2 = ReferenceFE(lagrangian,Float64,order+1)
    Ω = Triangulation(model);
    if isempty(boundary_tag)
        V1 = TestFESpace(Ω,reffe1,conformity=:Hcurl,vector_type=Vector{ComplexF64})
        V2 = TestFESpace(Ω,reffe2,conformity=:H1,vector_type=Vector{ComplexF64}) 
    else
        V1 = TestFESpace(Ω,reffe1,conformity=:Hcurl,vector_type=Vector{ComplexF64},dirichlet_tags=[boundary_tag])
        V2 = TestFESpace(Ω,reffe2,conformity=:H1,vector_type=Vector{ComplexF64},dirichlet_tags=[boundary_tag]) 
    end
    U1 = V1;
    U2 = V2;
    V = MultiFieldFESpace([V1, V2])
    U = MultiFieldFESpace([U1, U2])
    degree = 2*order+2;
    dΩ = Measure(Ω,degree);
    invmu=inverse(mu);
    tensor_C(x)=TensorValue(invmu.yy(x),-invmu.yx(x),-invmu.xy(x),invmu.xx(x))
    tensor_B(x)=TensorValue(invmu.yy(x),-invmu.xy(x),-invmu.yx(x),invmu.xx(x))
    vector_B1(x)=VectorValue(-invmu.zy(x),invmu.zx(x))
    vector_B2(x)=VectorValue(invmu.yz(x),-invmu.xz(x))
    
    a((Et1,Ez1),(Et2,Ez2))=∫( -curl(Et2)*(vector_B1⋅(∇(Ez1))+invmu.zz.f*curl(Et1))-curl(Et1)*(vector_B2⋅(∇(Ez2)))-(∇(Ez1))⋅(tensor_C⋅(∇(Ez2))) + k^2/eps0*(computeDt(Et1,Ez1,epsilon)⋅Et2+computeDz(Et1,Ez1,epsilon)*Ez2) )*dΩ
    b((Et1,Ez1),(Et2,Ez2))=∫( Et1⋅(tensor_C⋅(∇(Ez2)))-Et2⋅(tensor_B⋅(∇(Ez1)))+(vector_B1⋅Et1)*curl(Et2)+(vector_B2⋅Et2)*curl(Et1) )*dΩ
    c((Et1,Ez1),(Et2,Ez2))=∫( -(Et1⋅(tensor_C⋅Et2)) )*dΩ

    A_task=Threads.@spawn assemble_matrix(a,U,V);
    B_task=Threads.@spawn assemble_matrix(b,U,V);
    C_task=Threads.@spawn assemble_matrix(c,U,V);
    A=fetch(A_task);
    B=fetch(B_task);
    C=fetch(C_task);
    E,F=get_companion(A,im.*B,C);
    if (verbose)
        println("Matrices of dimension ",size(E), " created.")
    end
    if (solver==:LU)
        tmp1,tmp2=eigs_LU(F,E,sigma=approx_neff*k,nev=neigs,tol=tol,verbose=verbose);
    elseif (solver==:MUMPS)
        tmp1,tmp2=eigs_MUMPS(F,E,sigma=approx_neff*k,nev=neigs,tol=tol,verbose=verbose);
    else
        throw(DomainError(solver, "solver must be :LU or :MUMPS"));
    end
    neff=tmp1/k;
    if (verbose)
        println("Found ",length(neff)," eigenvalues.")
    end
    if (field)
        sol=Vector{Mode{VectorFieldFEM2D}}(undef,0);
        for i in axes(neff,1)
            name=string("Mode n°",string(i));
            Et,Ez=FEFunction(U,tmp2[:,i])
            ux=VectorValue(1.0,0.0);
            uy=VectorValue(0.0,1.0);
            Ex=Et⋅ux;
            Ey=Et⋅uy;
            Bt,Bz=computeB(Et,Ez,lambda,neff[i]);
            Ht,Hz=computeH(Bt,Bz,invmu);
            Hx=Ht⋅ux;
            Hy=Ht⋅uy;
            push!(sol,Mode(name,neff[i],lambda,VectorFieldFEM2D(dΩ,Ex,Ey,Ez,Hx,Hy,Hz)));
        end
    else
        sol=Vector{Mode{Nothing}}(undef,0);
        for i in axes(neff,1)
            name=string("Mode n°",string(i));
            push!(sol,Mode(name,neff[i],lambda));
        end
    end
    return sol
end

"""
    FEM2D_periodic(lambda::Real,eps_fonc::Function,model::DiscreteModel;approx_neff::Real=0,neigs::Int64=1,order::Int64=2,field::Bool=false,solver::Symbol=:LU,tol::Real=0.0,verbose::Bool=false,kt::AbstractVector{<:Real}=[0.0,0.0],type::Symbol=:Scalar)

Returns a vector of `Mode`.
The fiber is isotropic and described with its relative permittivity. The mesh must be periodic.

- lambda: wavelength
- epsilon: function of (x,y) that describes the relative permittivity profile
- model: DiscreteModel generated with `GridapGmsh.jl`
- approx_neff: effective index around which the modes will be computed
- neigs: number of modes to calculate
- order: order of the FEM
- field: boolean that indicates if fields must be saved
- solver: can be :LU or :MUMPS
- tol: tolerance for the eigenvalue solver (see documention of ArnoldiMethod.jl)
- verbose: boolean that enables some outputs
- kt: vector of the 2D Brillouin zone
- type: :Scalar or :Vector

"""
function FEM2D_periodic(lambda::Real,eps_fonc::Function,model::DiscreteModel;approx_neff::Real=0,neigs::Int64=1,order::Int64=3,field::Bool=false,solver::Symbol=:LU,tol::Real=0.0,verbose::Bool=false,kt::AbstractVector{<:Real}=[0.0,0.0],type::Symbol=:Scalar)
    if (lambda<=0)
        throw(DomainError(lamdba, "The wavelength lambda must be strictly positive"));
    end
    if (neigs<=0)
        throw(DomainError(neigs, "The number of modes must be strictly positive"));
    end
    if (first(methods(eps_fonc)).nargs!=2)
        throw(DomainError(eps_fonc, "The permittivity function must have 1 argument"));
    end
    if (num_dims(model)!=2)
        throw(DomainError(model, "The model must be 2D"));
    end
    if (approx_neff<0)
        throw(DomainError(approx_neff, "The approximative effective index must be positive or null"));
    end
    if (order<1)
        throw(DomainError(order, "The order must be at least 1"));
    end
    if (!(solver in [:LU,:MUMPS]))
        throw(DomainError(solver, "solver must be :LU or :MUMPS"));
    end
    if (tol<0)
        throw(DomainError(tol, "The tolerance must be positive or null"));
    end
    if (!(type in [:Scalar,:Vector]))
        throw(DomainError(solver, "solver must be :Scalar or :Vector"));
    end

    coord=Gridap.ReferenceFEs.get_node_coordinates(model.grid)
    if (approx_neff<=0)
        approx_neff=maximum(sqrt.(eps_fonc.(coord)))
        if (verbose)
            println("Set approx_neff to ",approx_neff)
        end
    end
    
    if (type==:Vector)
        k=2*pi/lambda;
        reffe1 = ReferenceFE(nedelec,order)
        reffe2 = ReferenceFE(lagrangian,Float64,order+1)
        Ω = Triangulation(model);
        ikx=VectorValue(im*kt[1],im*kt[2]);
        V1 = TestFESpace(Ω,reffe1,conformity=:Hcurl,vector_type=Vector{ComplexF64})
        V2 = TestFESpace(Ω,reffe2,conformity=:H1,vector_type=Vector{ComplexF64}) 
        U1 = V1;
        U2 = V2;
        V = MultiFieldFESpace([V1, V2])
        U = MultiFieldFESpace([U1, U2])
        degree = 2*order+2;
        dΩ = Measure(Ω,degree);
        invmu=tensor3(1.0)

        a((Et1,Ez1),(Et2,Ez2))=∫( -(curl(Et2)-ikx×Et2)*(curl(Et1)+ikx×Et1)-(∇(Ez1)+ikx*Ez1)⋅(∇(Ez2)-ikx*Ez2) + k^2*(Et1⋅(eps_fonc*Et2)+Ez1⋅(eps_fonc*Ez2)) )*dΩ
        b((Et1,Ez1),(Et2,Ez2))=∫( Et1⋅(∇(Ez2)-ikx*Ez2)-Et2⋅(∇(Ez1)+ikx*Ez1) )*dΩ
        c((Et1,Ez1),(Et2,Ez2))=∫( -(Et1⋅Et2) )*dΩ
    
        A_task=Threads.@spawn assemble_matrix(a,U,V);
        B_task=Threads.@spawn assemble_matrix(b,U,V);
        C_task=Threads.@spawn assemble_matrix(c,U,V);
        A=fetch(A_task);
        B=fetch(B_task);
        C=fetch(C_task);
        E,F=get_companion(A,im.*B,C);
        if (verbose)
            println("Matrices of dimension ",size(E), " created.")
        end
        if (solver==:LU)
            tmp1,tmp2=eigs_LU(F,E,sigma=approx_neff*k,nev=neigs,tol=tol,verbose=verbose);
        elseif (solver==:MUMPS)
            tmp1,tmp2=eigs_MUMPS(F,E,sigma=approx_neff*k,nev=neigs,tol=tol,verbose=verbose);
        else
            throw(DomainError(solver, "solver must be :LU or :MUMPS"));
        end
        neff=tmp1/k;
        if (verbose)
            println("Found ",length(neff)," eigenvalues.")
        end
        if (field)
            sol=Vector{Mode{VectorFieldFEM2D}}(undef,0);
            for i in axes(neff,1)
                name=string("Mode n°",string(i));
                Et,Ez=FEFunction(U,tmp2[:,i])
                ux=VectorValue(1.0,0.0);
                uy=VectorValue(0.0,1.0);
                Ex=Et⋅ux;
                Ey=Et⋅uy;
                Bt,Bz=computeB(Et,Ez,lambda,neff[i]);
                Ht,Hz=computeH(Bt,Bz,invmu);
                Hx=Ht⋅ux;
                Hy=Ht⋅uy;
                push!(sol,Mode(name,neff[i],lambda,VectorFieldFEM2D(dΩ,Ex,Ey,Ez,Hx,Hy,Hz)));
            end
        else
            sol=Vector{Mode{Nothing}}(undef,0);
            for i in axes(neff,1)
                name=string("Mode n°",string(i));
                push!(sol,Mode(name,neff[i],lambda));
            end
        end
    else
        k2=(2*pi/lambda)^2;
        reffe = ReferenceFE(lagrangian,Float64,order);
        V = TestFESpace(model,reffe,conformity=:H1,vector_type=Vector{ComplexF64})

        U=TrialFESpace(V,0);
        degree = 2*order;
        Ω = Triangulation(model);
        dΩ = Measure(Ω,degree);
        ikx=VectorValue(im*kt[1],im*kt[2]);
        A_task=Threads.@spawn assemble_matrix((u,v)->∫( -∇(v)⋅∇(u)+v*(ikx⋅∇(u))-∇(v)⋅(u*ikx)+(v*ikx)⋅(u*ikx) + k2*(v⋅(eps_fonc*u))  )*dΩ,U,V);
        B_task=Threads.@spawn assemble_matrix((u,v)->∫(v⋅u  )dΩ,U,V);
        A=fetch(A_task);
        B=fetch(B_task);
        if (verbose)
            println("Matrices of dimension ",size(A), " created.")
        end
        if (solver==:LU)
            tmp1,tmp2=eigs_LU(A,B,sigma=approx_neff^2*k2,nev=neigs,tol=Float64(tol),verbose=verbose);
        else
            tmp1,tmp2=eigs_MUMPS(A,B,sigma=approx_neff^2*k2,nev=neigs,tol=Float64(tol),verbose=verbose);
        end
        neff=sqrt.(tmp1/k2);
        if (verbose)
            println("Found ",length(neff)," eigenvalues.")
        end
        if (field)
            sol=Vector{Mode{ScalarFieldFEM2D}}(undef,0);
            E=[FEFunction(U,tmp2[:,i]) for i in axes(tmp2,2)];
            for i in axes(neff,1)
                name=string("Mode LP n°",string(i));
                push!(sol,Mode(name,neff[i],lambda,ScalarFieldFEM2D(dΩ,E[i])));
            end
        else
            sol=Vector{Mode{Nothing}}(undef,0);
            for i in axes(neff,1)
                name=string("Mode LP n°",string(i));
                push!(sol,Mode(name,neff[i],lambda));
            end
        end
    end
    return sol;
end

