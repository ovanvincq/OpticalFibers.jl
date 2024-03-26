#Bent fiber
using OpticalFibers
using OpticalFibers.ModeSolvers
using Bessels
using Plots
using GridapGmsh

lambda=1.55

ncore=1.4457
ncladding=1.4378
rho=3.5

delta=(ncore^2-ncladding^2)/(2*ncore^2)

V=2*pi*rho/lambda*sqrt(ncore^2-ncladding^2)

ms0=multi_step_fiber_modes(lambda,0,rho,[ncore,ncladding])
neff=ms0[1].neff

U=2*pi*rho/lambda*sqrt(ncore^2-neff^2)
W=2*pi*rho/lambda*sqrt(neff^2-ncladding^2)

#x bending radius in m
gamma=x->sqrt(pi)/2/(rho*1E-6)*sqrt(rho*1E-6/x)*U^2/V^2/W^(1.5)/(besselk(1,W)^2)*exp(-4/3*x/(rho*1E-6)*W^3/V^2*delta)

model = GmshDiscreteModel("./models/example5.msh");

x=(1:0.01:10)

#without bending
epsilon2D=x->(ncore-(ncore-ncladding)*(hypot(x[1],x[2])>=rho))^2
mFEM2D=FEM2D(lambda,epsilon2D,model,field=true,neigs=1);
neff0=mFEM2D[1].neff

#with bending
RcFEM=1:1:10
gammaFEM=zeros(length(RcFEM))
mFEM=Vector{Mode}(undef,length(RcFEM))
Threads.@threads for i in axes(RcFEM,1)
    epsilon2D_bent=x->epsilon2D(x)*(1+x[1]/RcFEM[i]*1E-3)^2
    mFEM2D_bent=FEM2D(lambda,epsilon2D_bent,model,neigs=1,dPML=3,approx_neff=neff0,field=true);
    mFEM[i]=mFEM2D_bent[1]
    gammaFEM[i]=2*imag(2*pi*mFEM2D_bent[1].neff/(lambda*1E-6))
end
Rc=1:0.01:10
plot(Rc,gamma.(Rc*1E-3),yaxis=:log,xlabel="Rc (mm)",ylabel="γ (m⁻¹)",label="Snyder & Love")
plot!(RcFEM,gammaFEM,label="FEM2D")

using GridapMakie
using GLMakie
using Gridap
fig=GLMakie.Figure(resolution=(750,250))
fig1,ax1=GLMakie.plot(fig[1,1],get_triangulation(mFEM[10]),abs(mFEM[10].field.E),axis=(aspect=DataAspect(),),colormap=:jet);
fig2,ax2=GLMakie.plot(fig[1,2],get_triangulation(mFEM[3]),abs(mFEM[3].field.E),axis=(aspect=DataAspect(),),colormap=:jet);
fig3,ax3=GLMakie.plot(fig[1,3],get_triangulation(mFEM[1]),abs(mFEM[1].field.E),axis=(aspect=DataAspect(),),colormap=:jet);
fig1.title="R = 10 mm"
fig2.title="R = 3 mm"
fig3.title="R = 1 mm"

R=3E3;
eps_anisotrope=tensor3(x->epsilon2D(x)*(1+x[1]/R),x->epsilon2D(x)*(1+x[1]/R),x->epsilon2D(x)/(1+x[1]/R))
mu_anisotrope=tensor3(x->(1+x[1]/R),x->(1+x[1]/R),x->1.0/(1+x[1]/R))
eps_anisotrope_pml=add_cylindrical_PML(eps_anisotrope,32,3,10)
mu_anisotrope_pml=add_cylindrical_PML(mu_anisotrope,32,3,10)

t=FEM2D_general(lambda,eps_anisotrope_pml,mu_anisotrope_pml,model,approx_neff=neff0,neigs=2,field=true,solver=:MUMPS)

#Twisted fiber
eps_twist=OpticalFibers.add_twist_PML(epsilon2D,2E3,32,3,10)
mu_twist=add_twist_PML(1,2E3,32,3,10)
t=FEM2D_anisotropic(lambda,eps_twist,mu_twist,model,approx_neff=neff0,neigs=20,field=true,solver=:MUMPS)

epsilon2D_shift=x->(ncore-(ncore-ncladding)*(hypot(x[1]-10,x[2])>=rho))^2
eps_shift_twist=OpticalFibers.add_twist_PML(epsilon2D_shift,2E3,32,3,10)
model_shift = GmshDiscreteModel("./models/example5-decale.msh");
t=FEM2D_anisotropic(lambda,eps_shift_twist,mu_twist,model_shift,approx_neff=neff0,neigs=20,field=true,solver=:MUMPS)