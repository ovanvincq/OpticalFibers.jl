#Leaky modes in step index fiber with a low index trench
#Circular Core
using OpticalFibers
using OpticalFibers.ModeSolvers
using Gridap
epsilon=x->(1.46-0.05*(x[1]>=4)+0.04*(x[1]>7))^2;
model = CartesianDiscreteModel((0,15),1500)
m0=FEM1D(1.6,0,epsilon,model,field=true,dPML=3,neigs=10)
losses(m0[5])*1e6

using Plots
normalize!.(m0)
r=0:0.01:15
plot(r,real.(computeField(m0[1],r)),xlabel="r (µm)",ylabel="real(E)",label="LP01")
plot!(r,-real.(computeField(m0[5],r)),label="LP02")

m1=FEM1D(1.6,1,epsilon,model,field=true,dPML=3,neigs=15);
pos=argmin(losses.(m1)*1E6);
m1[pos].neff

using GridapGmsh
model = GmshDiscreteModel("./models/Step_index_fiber_pml.msh");
epsilon2D=x->epsilon(hypot(x[1],x[2]));
m=FEM2D(1.6,epsilon2D,model,field=true,neigs=4,approx_neff=real(m1[pos].neff),dPML=3,type=:Vector)

using GridapMakie
using GLMakie
fig,ax,plot_obj=GLMakie.plot(get_triangulation(m[1]),real(m[1].field.Ex),axis=(aspect=DataAspect(),),colormap=:jet)
Colorbar(fig[1,2], plot_obj);

#Elliptical core
model2 = GmshDiscreteModel("./models/Elliptic_fiber_rectangular_pml.msh");
epsilon2=x->(1.46-0.05*(x[1]^2/16+x[2]^2/4>=1)+0.04*(x[1]^2/64+x[2]^2/16>1))^2
m_scalar=FEM2D(1.6,epsilon2,model2,neigs=2,approx_neff=1.44,field=true,solver=:MUMPS,dPML=3)

m_vector=FEM2D(1.6,epsilon2,model2,neigs=4,approx_neff=1.44,field=true,solver=:MUMPS,dPML=3,type=:Vector)