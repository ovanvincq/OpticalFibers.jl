# OpticalFibers.ModeSolvers - Tutorial

```@meta
CurrentModule = OpticalFibers.ModeSolvers
```

## Bimodal Step-index Fiber
This section explains how to modelize a step-index fiber with a core radius $a=2$ µm. The refractive index of the core is $n_{\text{core}}=1.47$ and that of the cladding is $n_{\text{cladding}}=1.45$. 

This fiber is bimodal at $\lambda=1$ µm since the normalized frequency is $V=\frac{2\pi a}{\lambda}\sqrt{n_{\text{core}}^2-n_{\text{cladding}}^2}=3.04$. 

### Scalar Modes
To compute the two scalar modes, we can use the fuction `multi_step_fiber_modes` that returns a vector of modes. The arguments are the wavelength, the azimuthal number, the radius of the core and a vector describing the refractive index. The optional argument `field` indicates that we want to return modes that contain a field.

```@example 1
#using Pkg; nothing # hide
#Pkg.activate("../.."); nothing # hide
using OpticalFibers, OpticalFibers.ModeSolvers
m0=multi_step_fiber_modes(1u"µm",0,2u"µm",[1.47,1.45],field=true);
m01=m0[1]
m1=multi_step_fiber_modes(1u"µm",1,2u"µm",[1.47,1.45],field=true);
m11=m1[1]
nothing; #hide
```

Note that you can also use broadcasting to compute all modes with a single command. The function `Ref` allows to use the same vector of refractive index for all values of the azimuthal number. 
```@example 1
m=multi_step_fiber_modes.(1u"µm",[0,1],2u"µm",Ref([1.47,1.45]),field=true);
m01=m[1][1];
m11=m[2][1];
nothing; #hide
```

The mode profile can be easily plotted:
```@example 1
using Plots
r=(0:0.01:10)*1u"µm";
plot(r,m01.EMField.E(r),label=m01.Name,xlabel="r",ylabel="E",linewidth=2,xlim=[0,10],ylim=[0,:auto])
plot!(r,m11.EMField.E(r),label=m11.Name,linewidth=2)
```

In order to visualize the modes in a 2D plot, we can convert the `Mode{ScalarFieldFunction1D}` to a `Mode{ScalarFieldFunction2D}` but we have to indicate the orientation of the mode by using the `cos` function or the `sin` function. Then, we normalize the modes and plot one of the two LP$_{11}$ modes. Be careful when using `Plots.jl`: unlike this package, the second index of the matrix corresponds to the x-coordinate (as matlab but the opposite of `Makie.jl`).
```@example 1
mm01=convertTo2D(m01)
mm11c=convertTo2D(m11)
mm11s=convertTo2D(m11,90)
mm01=normalize(mm01)
mm11c=normalize(mm11c)
mm11s=normalize(mm11s)
x=(-8:0.125:8)*1u"µm";
contourf(x,x,mm11c.EMField.E(x,x)',levels=100,linewidth=0,aspect_ratio=:equal,title=mm11c.Name,xlabel="x",ylabel="y",cbar=false,ylim=[-8,8],xlim=[-8,8],size=(400,400))
```

The beating between the LP$_{01}$ and the LP$_{11}$ can be observed by plotting the sum of the fields for different values of the distance $z$. The beating length is $\frac{2\pi}{\vert \Delta \beta \vert} = \frac{\lambda}{\vert \Delta n_{eff} \vert}\simeq 107$ µm
```@example 1
L=1u"µm"/(m01.neff-m11.neff)
```

We can take advantage of the possibility to add fields to create the animation.
```@example 1
maxi=maximum(abs2.(mm01.EMField.E(x,x)+mm11c.EMField.E(x,x)));
anim=@animate for j=0:107
    TotalField=getEMField(mm01,j*2u"µm")+getEMField(mm11c,j*2u"µm");
    contourf(x,x,abs2.(TotalField.E(x,x))',levels=100,linewidth=0,title="|E|² - z = $(2*j) µm",xlabel="x",ylabel="y",ylim=[-8,8],xlim=[-8,8],size=(400,400),cbar=false,clim=(0*maxi,maxi))
end;
gif(anim,"anim_field.gif",fps=15)
```

### Vector modes
It is also possible to compute the vector modes of the fiber: LP$_{01}$ mode becomes HE$_{11}$ mode and LP$_{11}$ mode becomes TE$_{01}$, TM$_{01}$ and HE$_{21}$ modes. To compute these modes, we just have to set the argument `type` to `:Vector` (its default value is `:Scalar`).

```@example 1
mv0=multi_step_fiber_modes(1u"µm",0,2u"µm",[1.47,1.45],field=true,type=:Vector)
```

```@example 1
mv1=multi_step_fiber_modes(1u"µm",1,2u"µm",[1.47,1.45],field=true,type=:Vector)
```

```@example 1
mv2=multi_step_fiber_modes(1u"µm",2,2u"µm",[1.47,1.45],field=true,type=:Vector)
```

Then, the Poynting vector of the mode HE$_{21}$ can be computed and plotted.
```@example 1
contourf(x,x,(mv2[1].EMField.Pz(x,x))',linewidth=0,levels=100,xlims=(-4u"µm",4u"µm"),ylims=(-4u"µm",4u"µm"),aspect_ratio=:equal,xlabel="x",ylabel="y",size=(400,400),cbar=false)
xc=x[1:4:end];
X,Y=meshgrid(xc);
quiver!(X,Y,quiver=(real(mv2[1].EMField.Ex(xc,xc)'/20u"V/m"),real(mv2[1].EMField.Ey(xc,xc)'/20u"V/m")),color=:cyan,arrow=arrow(:closed))
```

In order to check the orthoganality of the modes, we can normalize them and compute the overlap integrals:
```@example 1
mv=[mv1;mv0;mv2];
mv=normalize.(mv);
overlap.(mv,transpose(mv))
```

### Scalar modes with FEM1D
We can find the results above using the FEM method. We will start with the 1D method since the fiber has a cylindrical symmetry. The mesh must begin at $r=0$ and we have set its end at 15 µm. The number of elements is set to 1500. We indicate to the the FEM solver that we want to compute 4 eigenvalues but it returns only one since there is only one mode for $\ell=0$.
```@example 1
model = CartesianDiscreteModel((0,15),1500)
RIP=x->1.47-0.02*(x[1]>=2u"µm")
m0FEM=FEM1D(1u"µm",0,RIP,model*u"µm",field=true,neigs=4)
m0FEM=normalize.(m0FEM)
m0=normalize.(m0)
plot(r,m0[1].EMField.E(r),label="Step-index solver",ylabel="E",xlabel="r",linewidth=2,xlim=[0,10],ylim=[0,:auto])
plot!(r,-m0FEM[1].EMField.E(r),label="FEM1D solver",line=:dash,linewidth=2)
```

### Scalar modes with FEM2D
We can also use the FEM2D method. For this, we use a mesh created with GMSH and a function of the tuple $(x,y)$ that describes the refractive index profile. We ask the solver to find 4 modes but it returns only 3 because the fiber can only guides the modes LP01, LP11a and LP11b (a and b refer to the orientation of the mode).
```@example 1
model = GmshDiscreteModel("../../models/example1.msh");
RIP2D=x->1.47-0.02*(hypot(x[1],x[2])>=2u"µm")
mFEM2D=FEM2D(1*u"µm",RIP2D,model*u"µm",field=true,neigs=4)
mFEM2D=normalize.(mFEM2D);
p1=contourf(x,x,mFEM2D[1].EMField.E(x,x)',levels=50,linewidth=0,aspect_ratio=:equal,xlabel="x",ylabel="y",cbar=false,ylim=[-8,8],xlim=[-8,8],title="Mode 1");
p2=contourf(x,x,mFEM2D[2].EMField.E(x,x)',levels=50,linewidth=0,aspect_ratio=:equal,xlabel="x",ylabel="y",cbar=false,ylim=[-8,8],xlim=[-8,8],title="Mode 2");
p3=contourf(x,x,mFEM2D[3].EMField.E(x,x)',levels=50,linewidth=0,aspect_ratio=:equal,xlabel="x",ylabel="y",cbar=false,ylim=[-8,8],xlim=[-8,8],title="Mode 3");
plot(p1, p2, p3, layout=(1,3), legend=false,size=(900,310))
```
Note that there is no preferential orientation for the LP$_{11}$ modes given by the FEM solver. We can check that the modes given by the quasi-analytical solver and the FEM solver are the same by computing the overlap integrals:
```@example 1
abs2.(overlap.(mFEM2D,[mm01 mm11c mm11s]))
```

### Vector modes with FEM2D
To compute the vector, we just have to set `type` to `:Vector`. We ask the solver to find 8 modes but it returns only 6 because the fiber can only guides 2 modes HE$_{11}$, 1 mode becomes TE$_{01}$, 1 mode TM$_{01}$ and 2 modes HE$_{21}$. 
```@example 1
mvFEM2D=FEM2D(1*u"µm",RIP2D,model*u"µm",field=true,neigs=8,type=:Vector)
mvFEM2D=normalize.(mvFEM2D)
```
We can verify that the modes given by the quasi-analytical solver and the FEM solver are the same by plotting the fields.
```@example 1
p1=contourf(x,x,real(mvFEM2D[4].EMField.Ex(x,x)'),levels=50,linewidth=0,aspect_ratio=:equal,title="Ex",xlabel="x",ylabel="y",cbar=false,ylim=[-8,8],xlim=[-8,8]);
p2=contourf(x,x,real(mvFEM2D[4].EMField.Ey(x,x)'),levels=50,linewidth=0,aspect_ratio=:equal,title="Ey",xlabel="x",ylabel="y",cbar=false,ylim=[-8,8],xlim=[-8,8]);
p3=contourf(x,x,imag(mvFEM2D[4].EMField.Ez(x,x)'),levels=50,linewidth=0,aspect_ratio=:equal,title="Ez",xlabel="x",ylabel="y",cbar=false,ylim=[-8,8],xlim=[-8,8]);
p4=contourf(x,x,real(mvFEM2D[4].EMField.Hx(x,x)'),levels=50,linewidth=0,aspect_ratio=:equal,title="Hx",xlabel="x",ylabel="y",cbar=false,ylim=[-8,8],xlim=[-8,8]);
p5=contourf(x,x,real(mvFEM2D[4].EMField.Hy(x,x)'),levels=50,linewidth=0,aspect_ratio=:equal,title="Hy",xlabel="x",ylabel="y",cbar=false,ylim=[-8,8],xlim=[-8,8]);
p6=contourf(x,x,imag(mvFEM2D[4].EMField.Hz(x,x)'),levels=50,linewidth=0,aspect_ratio=:equal,title="Hz",xlabel="x",ylabel="y",cbar=false,ylim=[-8,8],xlim=[-8,8]);
plot(p1, p2, p3, p4, p5, p6, layout=(2,3), legend=false,size=(900,600))
```
```@example 1
p1=contourf(x,x,real(mv[4].EMField.Ex(x,x)'),levels=50,linewidth=0,aspect_ratio=:equal,title="Ex",xlabel="x",ylabel="y",cbar=false,ylim=[-8,8],xlim=[-8,8]);
p2=contourf(x,x,real(mv[4].EMField.Ey(x,x)'),levels=50,linewidth=0,aspect_ratio=:equal,title="Ey",xlabel="x",ylabel="y",cbar=false,ylim=[-8,8],xlim=[-8,8]);
p3=contourf(x,x,imag(mv[4].EMField.Ez(x,x)'),levels=50,linewidth=0,aspect_ratio=:equal,title="Ez",xlabel="x",ylabel="y",cbar=false,ylim=[-8,8],xlim=[-8,8]);
p4=contourf(x,x,real(mv[4].EMField.Hx(x,x)'),levels=50,linewidth=0,aspect_ratio=:equal,title="Hx",xlabel="x",ylabel="y",cbar=false,ylim=[-8,8],xlim=[-8,8]);
p5=contourf(x,x,real(mv[4].EMField.Hy(x,x)'),levels=50,linewidth=0,aspect_ratio=:equal,title="Hy",xlabel="x",ylabel="y",cbar=false,ylim=[-8,8],xlim=[-8,8]);
p6=contourf(x,x,imag(mv[4].EMField.Hz(x,x)'),levels=50,linewidth=0,aspect_ratio=:equal,title="Hz",xlabel="x",ylabel="y",cbar=false,ylim=[-8,8],xlim=[-8,8]);
plot(p1, p2, p3, p4, p5, p6, layout=(2,3), legend=false,size=(900,600))
```

## Gradient index fiber and dispersion
In this tutorial, a germanium-doped parabolic gradient index fiber will be studied. The maximum Ge concentration is 20% and the core radius is 3.5 µm. The silica dispersion will be taken into account.

First, a vector of `Function` is created to modelized the dispersive refractive index profile between 1 and 1.5 µm:
```@example 2
using Plots, OpticalFibers, OpticalFibers.PhysicalData, OpticalFibers.ModeSolvers
lambda=(1:0.01:1.5)u"µm";
RIP=[x->n_Ge_Doped_Silica(l,0)+(n_Ge_Doped_Silica(l,0.2)-n_Ge_Doped_Silica(l,0))*(x[1]<=3.5u"µm")*(1-x[1]^2/(3.5u"µm")^2) for l in lambda];
r=(0:0.1:5)u"µm";
plot(r,RIP[1].(r),label="λ = 1 µm",xlabel="r",ylabel="Refractive index",linewidth=2,xlim=[0,5])
plot!(r,RIP[end].(r),label="λ = 1.5 µm",linewidth=2)
```

The modal content is then computed (note that we must use Ref(model) because length(model) is not defined in Gridap):
```@example 2
model = CartesianDiscreteModel((0,20),2000)
m=FEM1D.(lambda,[0,1,2]',RIP,Ref(model*u"µm"),neigs=2);
#mode LP01 always exists
neff01=[real(m[j,1][1].neff) for j in 1:length(lambda)];
N02=sum((length.(m[:,1])).>=2);
neff02=[real(m[j,1][2].neff) for j in 1:N02];
N11=sum((length.(m[:,2])).>=1);
neff11=[real(m[j,2][1].neff) for j in 1:N11];
N21=sum((length.(m[:,3])).>=1);
neff21=[real(m[j,3][1].neff) for j in 1:N21];
plot(lambda,n_Ge_Doped_Silica.(lambda,0),label="Silica",xlabel="λ",ylabel="Effective index",color=:black,linewidth=2,xlim=[1,1.5]);
plot!(lambda,n_Ge_Doped_Silica.(lambda,0.2),label="Ge-doped Silica (20%)",color=:black,line=:dash,linewidth=2);
plot!([lambda,lambda[1:N11],lambda[1:N02],lambda[1:N21]],[neff01,neff11,neff02,neff21],label=["LP01" "LP11" "LP02" "LP21"],linewidth=[2 2 2 2])
```

The second-order dispersion is defined by $\beta_2=\frac{\partial^2 \beta}{\partial \omega^2}$.
```@example 2
beta01=neff01*2*pi./lambda;
omega=2*pi*OpticalFibers.PhysicalData.c_Unitful./lambda
omega2,beta2=derivative((omega,beta01),2);
lambda2=2*pi*OpticalFibers.PhysicalData.c_Unitful./omega2;
plot(lambda2,uconvert.(u"ps^2/km",beta2),xlabel="λ",ylabel="β₂",label="LP01",linewidth=2,xlim=[1,1.5])
```

To compute the effective area and the non-linear coefficient of the fundamental mode, the fields must be calculated:
```@example 2
m=FEM1D.(lambda,0,RIP,Ref(model*u"µm"),field=true);
m=[m[i][1] for i in 1:length(lambda)];
A=Aeff.(m);
gamma=nonLinearCoefficient.(m,2.53E-20u"m^2/W")
plot(lambda,A,label="Aeff",ylabel="Effective area",xlabel="λ",color=:blue,leg=:topright,linewidth=2,xlim=[1,1.5]);
plot!(twinx(),lambda,uconvert.(u"W^-1*km^-1",gamma),label="γ",ylabel="Non-linear coefficient",color=:red,leg=:topleft,linewidth=2,xlim=[1,1.5])
```

## Leaky modes in step index fiber with a low index trench

### Circular core
In this example, we study a fiber with the parameters below:
- core with a radius of 4 µm and a refractive index of 1.46
- a trench located between 4 µm and 7 µm with a refractive index of 1.41
- a cladding with a refractive index of 1.45

In order to compute the leaky scalar modes at $\lambda=1.6$ µm, a PML must be added when using a FEM solver. In this example, the PML is located between 12 µm and 15 µm (thickness of 3 µm).
```@example 5
using Plots, OpticalFibers, OpticalFibers.ModeSolvers
RIP=x->1.46-0.05*(x[1]>=4u"µm")+0.04*(x[1]>7u"µm");
model = CartesianDiscreteModel((0,15),1500)
m0=FEM1D(1.6u"µm",0,RIP,model*u"µm",field=true,dPML=3u"µm",neigs=10)
```
The LP$_{01}$ mode is the mode 1. Since its effective index is greater than the refractive index of the cladding, this mode is guided and its losses are not significant.

The LP$_{02}$ mode is the mode 5, its losses in dB/km can be calculated:
```@example 5
m0[5].losses
```

```@example 5
r=(0:0.01:15)u"µm"
m0=normalize.(m0)
plot(r,real(m0[1].EMField.E(r)),xlabel="r",ylabel="real(E)",label="LP01",linewidth=2,xlim=[0,15])
plot!(r,-real(m0[5].EMField.E(r)),label="LP02",linewidth=2)
```

We can also compute the effective index of the LP$_{11}$ mode:
```@example 5
m1=FEM1D(1.6u"µm",1,RIP,model*u"µm",field=true,dPML=3u"µm",neigs=15);
pos=argmin(getproperty.(m1,:losses))
m1[pos].neff
```

To compute the four leaky vector modes that correspond to the LP$_{11}$ mode, we can also add a PML to the FEM2D solver:
```@example 5
model = GmshDiscreteModel("../../models/Step_index_fiber_pml.msh");
RIP2D=x->RIP(hypot(x[1],x[2]));
m=FEM2D(1.6u"µm",RIP2D,model*u"µm",field=true,neigs=4,approx_neff=real(m1[pos].neff),dPML=3u"µm",type=:Vector)
```

The field can be plotted using GridapMakie:
```@example 5
using GridapMakie, GLMakie
fig,ax,plot_obj=GLMakie.plot(real(m[1].EMField.Ex.value),axis=(aspect=DataAspect(),),colormap=:jet)
Colorbar(fig[1,2], plot_obj);
save("FEM_PML_Ex.png",fig); nothing #hide
```
![Real(Ex) for the first mode](FEM_PML_Ex.png)

### Elliptical core 
In this example, the fiber is similar to the previous one but the core and the trench are elliptic:
- core with a radius of 4 µm in the x direction and 2 µm in the y direction and a refractive index of 1.46
- a trench with a refractive index of 1.41 and located between the core and an ellipse with radii 8 and 4 µm in x and y direction respectively
- a cladding with a refractive index of 1.45

![Elliptic Fiber RIP](./assets/Elliptic_profile.png)

A rectangular PML is added for $12<\vert x \vert < 15$ µm and $8<\vert x \vert < 11$ µm to compute the leaky modes.

We can compute the modes in the scalar approximation or the vector modes. The results are similar since the refractive index steps are very low.

```@example 5
model2 = GmshDiscreteModel("../../models/Elliptic_fiber_rectangular_pml.msh");
RIP2=x->1.46-0.05*(x[1]^2/16u"µm^2"+x[2]^2/4u"µm^2">=1)+0.04*(x[1]^2/64u"µm^2"+x[2]^2/16u"µm^2">1)
using MUMPS
m_scalar=FEM2D(1.6u"µm",RIP2,model2*u"µm",neigs=2,approx_neff=1.44,field=true,solver=:MUMPS,dPML=3u"µm")
```
```@example 5
m_vector=FEM2D(1.6u"µm",RIP2,model2*u"µm",neigs=4,approx_neff=1.44,field=true,solver=:MUMPS,dPML=3u"µm",type=:Vector)
```

The first mode can be saved in a file that can be opened with ParaView:
```@example 5
writevtk("mode1",m_vector[1]);
```

![Elliptic Fiber Ex](./assets/Elliptic_Ex.png)

## Bent fiber
In this example, we will compute the bending losses of a step-index fiber at a wavelength of 1.55 µm. The core of the fiber has a refractive index of 1.4457 and a radius of 3.5 µm. The refractive index of the cladding is 1.4378.

First, we can use the step-index solver to calculate the effective index. Then, the approximate formula given in the Snyder & Love [SnyderLove:1983](@cite) (Formula 23-23, page 481) can be used to compute the attenuation coefficient $\gamma=2*Im(\beta)$:

```@example 6
using Bessels, OpticalFibers, OpticalFibers.ModeSolvers
lambda=1.55u"µm"
ncore=1.4457
ncladding=1.4378
rho=3.5u"µm"
delta=(ncore^2-ncladding^2)/(2*ncore^2)
V=2*pi*rho/lambda*sqrt(ncore^2-ncladding^2)
ms0=multi_step_fiber_modes(lambda,0,rho,[ncore,ncladding])
neff=ms0[1].neff
U=2*pi*rho/lambda*sqrt(ncore^2-neff^2)
W=2*pi*rho/lambda*sqrt(neff^2-ncladding^2)
gamma=x->sqrt(pi)/2/rho*sqrt(rho/x)*U^2/V^2/W^(1.5)/(besselk(1,W)^2)*exp(-4/3*x/rho*W^3/V^2*delta)
neff
```
To compute losses with FEM in the scalar approximation and neglecting the elasto-optics effect, we must add a PML and multiply the refractive index by $(1+x/Rc)$ with $Rc$ the bending radius. Note that the cladding must be large enough to contain the turning point (point at which the refractive index equals the effective index): in the case studied here, the cladding outer radius must be at least 25 µm for a bent radius of 10 mm. Then, we construct a model with a PML located between 32 and 35 µm.

We first compute the effective index with the FEM2D function to check that we obtain the same value as with the step-index solver.
```@example 6
model = GmshDiscreteModel("../../models/example5.msh");
RIP2D=x->ncore-(ncore-ncladding)*(hypot(x[1],x[2])>=rho)
mFEM2D=FEM2D(lambda,RIP2D,model*u"µm",field=true,neigs=1);
neff0=mFEM2D[1].neff
```

Then we calculate modes and the attenuation coefficent for Rc between 1 and 10 mm. Note that we use threads to accelerate the computation (this is not possible when using the MUMPS solver because it uses its own threads system and this causes an error).
```@example 6
RcFEM=(1:1:10)u"mm"
gammaFEM=zeros(length(RcFEM))*u"m^-1"
mFEM=Vector{Mode}(undef,length(RcFEM))
Threads.@threads for i in axes(RcFEM,1)
    RIP2D_bent=x->RIP2D(x)*(1+x[1]/RcFEM[i])
    mFEM2D_bent=FEM2D(lambda,RIP2D_bent,model*u"µm",neigs=1,dPML=3u"µm",approx_neff=neff0,field=true);
    mFEM[i]=mFEM2D_bent[1]
    gammaFEM[i]=mFEM2D_bent[1].alpha;
end
using Plots
Rc=(1:0.01:10)u"mm"
plot(Rc,uconvert.(u"m^-1",gamma.(Rc)),yaxis=:log,xlabel="Rc",ylabel="γ",label="Snyder & Love",linewidth=2,xlim=[0,10])
plot!(RcFEM,uconvert.(u"m^-1",gammaFEM),label="FEM2D",linewidth=2)
```
The value given by FEM2D is slighlty different from that predicted by the analytical formula and the difference increases as the bent radius decreases because the mode increasingly distorted.
```@example 6
using GridapMakie
using GLMakie
fig=GLMakie.Figure(size=(750,250))
fig1,ax1=GLMakie.plot(fig[1,1],abs(mFEM[10].EMField.E.value),axis=(aspect=DataAspect(),),colormap=:jet);
fig2,ax2=GLMakie.plot(fig[1,2],abs(mFEM[3].EMField.E.value),axis=(aspect=DataAspect(),),colormap=:jet);
fig3,ax3=GLMakie.plot(fig[1,3],abs(mFEM[1].EMField.E.value),axis=(aspect=DataAspect(),),colormap=:jet);
fig1.title="Rc = 10 mm"
fig2.title="Rc = 3 mm"
fig3.title="Rc = 1 mm"
save("Bent_fiber.png",fig); nothing #hide
```
![Electric field in the bent fiber](Bent_fiber.png)

In the case of vector modes, we can also multiply the refractive index by $(1+x/Rc)$ but, more rigourously, we have to multiply the permittivity and permeability tensors by $\left(\begin{matrix} (1+x/Rc) & 0 & 0\\
0 & (1+x/Rc) & 0 \\ 0 & 0 & (1+x/Rc)^{-1} \end{matrix}\right)$. Below, we compute the modes of the fiber with a bent radius of 3 mm.

```@example 6
R=3u"mm";
epsilon2D=x->RIP2D(x)^2
eps_anisotrope=x->diagonal_tensor(VectorValue(epsilon2D(x)*(1+x[1]/R),epsilon2D(x)*(1+x[1]/R),epsilon2D(x)/(1+x[1]/R)))
mu_anisotrope=x->diagonal_tensor(VectorValue(1+x[1]/R,1+x[1]/R,1.0/(1+x[1]/R)))
eps_anisotrope_pml=add_cylindrical_PML(eps_anisotrope,32u"µm",3u"µm",10)
mu_anisotrope_pml=add_cylindrical_PML(mu_anisotrope,32u"µm",3u"µm",10)
t=FEM2D_anisotropic(lambda,eps_anisotrope_pml,mu_anisotrope_pml,model*u"µm",approx_neff=neff0,neigs=2,field=true,solver=:MUMPS)
```

## Twisted fiber
The method used to compute the modes of an isotropic fiber is described in a paper written by Nicolet et al. [Nicolet2007](@cite). In the following example, we use the same fiber as in the bent fiber example but instead of being bent, the fiber is twisted with a period of $P=2$ mm.
```@example 6
eps_twist=OpticalFibers.add_twist_PML(epsilon2D,2u"mm",32u"µm",3u"µm",10)
mu_twist=add_twist_PML(1,2u"mm",32u"µm",3u"µm",10)
using MUMPS
t=FEM2D_anisotropic(lambda,eps_twist,mu_twist,model*u"µm",approx_neff=neff0,neigs=20,field=true,solver=:MUMPS)
```
Modes 2 and 20 are the HE$_{11}$ modes. The losses are not significant since the twist has no effect on a circular centered core. However, the effective indices are no more degenerate: as explained by Napiorkowski et al. [Napiorkowski:2014](@cite), in the helicoidal coordinates system, the effective index is neff$\pm\frac{\nu \lambda}{P}$ with $\nu$ the azimuthal number ($\nu=1$ for the HE$_{11}$ mode).

If the core is shifted by 10 µm, losses due to the twist appear:
```@example 6
epsilon2D_shift=x->(ncore-(ncore-ncladding)*(hypot(x[1]-10u"µm",x[2])>=rho))^2
eps_shift_twist=OpticalFibers.add_twist_PML(epsilon2D_shift,2u"mm",32u"µm",3u"µm",10)
model_shift = GmshDiscreteModel("../../models/example5-decale.msh");
t=FEM2D_anisotropic(lambda,eps_shift_twist,mu_twist,model_shift*u"µm",approx_neff=neff0,neigs=20,field=true,solver=:MUMPS)
```

## Photonic crystal fiber (PCF)
In a PCF, the modes are not guided modes but leaky modes so that the computation requires a PML. The fiber is constituted of three rings of air hole (n=1) inserted in silica (n=1.45). The pitch is 2 µm, the hole diameter is 1.5 µm and the PML begins at 8 µm from the fiber center and its thickness is 2 µm. 

![PCF mesh](./assets/PCF.png)

First, the mesh is loaded:
```@example 4
using OpticalFibers, OpticalFibers.ModeSolvers, GridapMakie, GLMakie
model = GmshDiscreteModel("../../models/PCF.msh");
```
Then we define the refractive index profile function:
```@example 4
Pitch=2u"µm"
r_hole=0.75u"µm"
pos=ring.(1:3)
xc=vcat(first.(pos)...)*Pitch
yc=vcat(last.(pos)...)*Pitch
RIP_PCF(x) = (any(@. hypot(x[1]-xc,x[2]-yc)<r_hole)) ? 1.0 : 1.45
```
Then we can compute four modes whose effective indices are close to the approximate value calculated for the fundamental mode at the wavelength of 1.3 µm:
```@example 4
neff_approx=approx_neff_PCF(1.3u"µm",1.5u"µm",2u"µm");
using MUMPS
m=FEM2D(1.3u"µm",RIP_PCF,model*u"µm",neigs=4,approx_neff=neff_approx,field=true,solver=:MUMPS,type=:Vector,dPML=2u"µm")
```
The last two modes are fundamental modes. We can compute and plot the z-component of the Poynting vector of the last mode:
```@example 4
Pz=m[end].EMField.Pz
fig,ax,plot_obj=GLMakie.plot(Pz.value,axis=(aspect=DataAspect(),),colormap=:jet)
ax.xlabel="x (µm)";
ax.ylabel="y (µm)";
ax.title="neff = $(m[end].neff)";
Colorbar(fig[1,2], plot_obj);
save("FEM2_Pz.png",fig); nothing #hide
```
![Pz for FM computed with FEM](FEM2_Pz.png)

## Photonic bandgap fiber
Before designing a PBG fiber, one has to compute the PBG of the infinite microstructured media that constitutes the cladding. In this example the cladding is a hexagonal lattice of circular rods with a pitch of 10 µm. The rods has a diameter of 3 µm and a refractive index of 1.47 while the cladding background media has a refractive index of 1.45. To compute the PBG, the mesh must be periodic. Since the cell is highly symmetric, the computation of the bands can be restricted to the highest symmetry points of the irreducible Brillouin zone (Γ, M and K).
```@example 7
using Plots, OpticalFibers, OpticalFibers.ModeSolvers
kt,weight=compute_kt(2,:Hexagon,Irreducible=true,MeshType=:Edge,Pitch=10u"µm")
RIP2D=x->1.47-0.02*(hypot(x[1],x[2])>=3u"µm")
model = GmshDiscreteModel("../../models/Cell_Hexagon2.msh")
lambda=(0.8:0.025:1.6)u"µm"
mm=Matrix{Vector{Mode}}(undef,length(lambda),length(kt))
Threads.@threads for i in 1:length(kt)
    mm[:,i]=FEM2D_periodic.(lambda,RIP2D,Ref(model*u"µm"),neigs=30,field=true,kt=kt[i]);
end
neff=zeros(length(lambda),length(kt),30)
P=plot()
a=palette([:red,:green,:blue],length(kt))
plot!(P,lambda,1.45*ones(length(lambda)),color=:black,label="",xlim=[0.8,1.6],ylim=[1.43,1.47],xlabel="Wavelength",ylabel="Effective index")
for i=1:30
    neff[:,:,i]=[real(m[i].neff) for m in mm]
    for k=1:length(kt)
        plot!(P,lambda,neff[:,k,i],color=a[k],label="")
    end
end
P.series_list[2].plotattributes[:label]="Γ"
P.series_list[3].plotattributes[:label]="M"
P.series_list[4].plotattributes[:label]="K"
P
```

To compute the fundamental mode, it is necessary to compute an approximative value of the effective index based on the band diagram.
```@example 7
using Interpolations
neff_approx=(max.(neff[5:18,2,7],neff[5:18,1,7])+min.(neff[5:18,2,6],1.45))/2
interp=LinearInterpolation(lambda[5:18],neff_approx)
model_PBG = GmshDiscreteModel("../../models/PBG.msh")
Pitch=10u"µm"
r_hole=3u"µm"
pos=ring.(1:4)
xc=vcat(first.(pos)...)*Pitch
yc=vcat(last.(pos)...)*Pitch
RIP_PBG(x) = (any(@. hypot(x[1]-xc,x[2]-yc)<r_hole)) ? 1.47 : 1.45
lambda_PBG=(0.925:0.005:1.115)u"µm"
mode_PBG=Vector{Vector{Mode}}(undef,length(lambda_PBG))
for i=1:length(lambda_PBG)
    mode_PBG[i]=FEM2D(lambda_PBG[i],RIP_PBG,model_PBG*u"µm",neigs=1,approx_neff=interp(lambda_PBG[i]),field=true,solver=:MUMPS,dPML=4u"µm")
end
mode_PBG=vcat(mode_PBG...)
neff_PBG=getproperty.(mode_PBG,:neff)
Aeff_PBG=Aeff.(mode_PBG)
losses_PBG= getproperty.(mode_PBG,:losses)
omega=2*pi*OpticalFibers.PhysicalData.c_Unitful./lambda_PBG
beta=real(getproperty.(mode_PBG,:beta))
omega2,beta2=derivative((omega,beta),2)
lambda2=2*pi*OpticalFibers.PhysicalData.c_Unitful./omega2
p1=plot(lambda_PBG,ustrip.(losses_PBG),yaxis=:log,xlabel="λ",ylabel="Losses (dB/km)",linewidth=2,xlim=[0.925,1.115])
p2=plot(lambda_PBG,Aeff_PBG,xlabel="λ",ylabel="Aeff",linewidth=2,xlim=[0.925,1.115])
p3=plot(lambda2,uconvert.(u"ps^2/km",beta2),xlabel="λ",ylabel="β₂",linewidth=2,xlim=[0.925,1.115])
plot(p1, p2, p3, layout=(1,3), legend=false,size=(650,400))
```
