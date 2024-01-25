# OpticalFibers.ModeSolvers - Tutorial

```@meta
CurrentModule = OpticalFibers.ModeSolvers
```

## Bimodal Step-index Fiber
This section explains how to modelize a step-index fiber with a core radius $a=2$ µm. The refractive index of the core is $n_{\text{core}}=1.47$ and that of the cladding is $n_{\text{cladding}}=1.45$. 

This fiber is bimodal at $\lambda=1$ µm since the normalized frequency is $V=\frac{2\pi a}{\lambda}\sqrt{n_{\text{core}}^2-n_{\text{cladding}}^2}=3.04$. 

To compute the two modes, we can use the fuction `multi_step_fiber_modes` that returns a vector of modes:

```@example 1
#using Pkg; nothing # hide
#Pkg.activate("../.."); nothing # hide
using OpticalFibers
using OpticalFibers.ModeSolvers
m0=multi_step_fiber_modes(1,0,2,[1.47,1.45],maxPosition=10);
m01=m0[1]
m1=multi_step_fiber_modes(1,1,2,[1.47,1.45],maxPosition=10);
m11=m1[1]
nothing; #hide
```

Note that you can also use broadcasting to compute all modes with a single command:
```@example 1
m=multi_step_fiber_modes.(1,[0,1],2,Ref([1.47,1.45]),maxPosition=10);
m01=m[1][1];
m11=m[2][1];
nothing; #hide
```

It is possible to plot the mode profile:
```@example 1
using Plots
plot(m01.r,m01.E,label=m01.Name)
plot!(m11.r,m11.E,label=m11.Name)
```

In order to visualize the modes in a 2D plot, one should use the conversion from `ScalarMode1D` to `ScalarMode2D`. Be careful when using `Plots.jl`: unlike this package, the second index of the matrix corresponds to the x-coordinate (as matlab but the opposite of `Makie.jl`).
```@example 1
mm01=ScalarMode2D(m01);
mm11s=ScalarMode2D(m11,sincos='s');
mm11c=ScalarMode2D(m11,sincos='c');
contourf(mm11c.x,mm11c.y,mm11c.E',levels=100,linewidth=0)
```

To observe the beating between the LP01 and the LP11, the modes must be first normalized. Then, the addition of the two fields at the distance z ∈ [0,214] µm is required. The beating length is $\frac{2\pi}{\vert \Delta \beta \vert} = \frac{\lambda}{\vert \Delta n_{eff} \vert}= 107$ µm
```@example 1
normalize!(mm01)
normalize!(mm11c)
normalize!(mm11s)
L=1/(m01.neff-m11.neff)
```

```@example 1
anim=@animate for j=0:214
    TotalField=ScalarField(mm01,j)+ScalarField(mm11c,j);
    contourf(TotalField.x,TotalField.y,abs2.(TotalField.E'),levels=100,linewidth=0)
    title!("z = $j µm");
end;
gif(anim,"anim_field.gif",fps=15)
```

It is also possible to compute the vector modes of the fiber: LP$_{01}$ mode becomes HE$_{11}$ mode and LP$_{11}$ mode becomes TE$_{01}$, TM$_{01}$ and HE$_{21}$ modes:

```@example 1
mv0=multi_step_fiber_modes(1,0,2,[1.47,1.45],maxPosition=10,type=:Vector)
```

```@example 1
mv1=multi_step_fiber_modes(1,1,2,[1.47,1.45],maxPosition=10,type=:Vector)
```

```@example 1
mv2=multi_step_fiber_modes(1,2,2,[1.47,1.45],maxPosition=10,type=:Vector)
```

Then, the Poynting Vector of the mode HE$_{21}$ can be computed and plotted:
```@example 1
Px,Py,Pz=PoyntingVector(mv2[1]);
contourf(mv2[1].x,mv2[1].y,Pz',linewidth=0,levels=100,xlims=(-4,4),ylims=(-4,4))
X,Y=meshgrid(mv2[1].x,mv2[1].y);
quiver!(X[5:5:end,5:5:end],Y[5:5:end,5:5:end],quiver=(mv2[1].Ex[5:5:end,5:5:end]'/20,mv2[1].Ey[5:5:end,5:5:end]'/20),color=:cyan,arrow=arrow(:closed))
```

In order to check the orthoganality of the modes, we can normalize them and compute the overlap integrals:
```@example 1
mv=[mv0;mv1;mv2];
normalize!.(mv);
overlap.(mv,transpose(mv))
```

## Gradient index fiber
In this tutorial, a germanium-doped parabolic gradient index fiber will be studied. The maximum Ge concentration is 20% and the core radius is 3.5 µm. The silica dispersion will be taken into account.

First, a vector of `Function` is created to modelized the dispersive refractive index profile between 1 and 1.5 µm:
```@example 2
using OpticalFibers
using OpticalFibers.PhysicalData
using Plots
lambda=1:0.01:1.5;
f=[r->n_Ge_Doped_Silica(l*1e-6,0)+(n_Ge_Doped_Silica(l*1e-6,0.2)-n_Ge_Doped_Silica(l*1e-6,0))*(r<=3.5)*(1-r^2/3.5^2) for l in lambda];
r=0:0.1:5;
plot(r,f[1].(r),label="λ = 1 µm",xlabel="r (µm)",ylabel="Refractive index");
plot!(r,f[end].(r),label="λ = 1.5 µm")
```

The modal content is then computed:
```@example 2
using OpticalFibers.ModeSolvers
m=FD.(lambda,[0,1,2]',2,f,1000,10);
#mode LP01 always exists
neff01=[m[j,1][1].neff for j in 1:length(lambda)];
N02=sum((length.(m[:,1])).>=2);
neff02=[m[j,1][2].neff for j in 1:N02];
N11=sum((length.(m[:,2])).>=1);
neff11=[m[j,2][1].neff for j in 1:N11];
N21=sum((length.(m[:,3])).>=1);
neff21=[m[j,3][1].neff for j in 1:N21];
plot(lambda,n_Ge_Doped_Silica.(lambda*1E-6,0),label="Silica",xlabel="λ (µm)",ylabel="Effective index",color=:black);
plot!(lambda,n_Ge_Doped_Silica.(lambda*1E-6,0.2),label="Ge-doped Silica (20%)",color=:black,line=:dash);
plot!([lambda,lambda[1:N11],lambda[1:N02],lambda[1:N21]],[neff01,neff11,neff02,neff21],label=["LP01" "LP11" "LP02" "LP21"])
```

The second-order dispersion is defined by $\beta_2=\frac{\partial^2 \beta}{\partial \omega^2}$.
```@example 2
beta01=neff01*2*pi./lambda*1E6;
omega=2*pi*OpticalFibers.PhysicalData.c./lambda*1E6;
omega2,beta2=derivative((omega,beta01),2);
lambda2=2*pi*OpticalFibers.PhysicalData.c./omega2*1E6;
plot(lambda2,beta2*1E26,xlabel="λ (µm)",ylabel="β₂ (10⁻²⁶ s²/m)")
```

To compute the effective area and the non-linear coefficient, the fields must be calculated:
```@example 2
m=FD.(lambda,0,1,f,1000,10,field=true);
m=[m[i][1] for i in 1:length(lambda)];
A=Aeff.(m);
gamma=nonLinearCoefficient.(m,2.53E-20)*1E21;
plot(lambda,A,label="Aeff",ylabel="Effective area (µm²)",xlabel="λ (µm)",color=:blue,leg=:topright);
plot!(twinx(),lambda,gamma,label="γ",ylabel="Non-linear coefficient ((W.km)⁻¹)",color=:red,leg=:topleft)
```

## FEM: Step-index fiber
As an example, we will compute the five first vector modes in a step-index fiber with a 1 µm-core radius, a refractive index of $\sqrt{3}$ in the core and 1 in the cladding. The mesh was generated with GMSH.
```@example 3
using OpticalFibers
using OpticalFibers.ModeSolvers
using Gridap
using GridapGmsh
using GridapMakie
using GLMakie
model = GmshDiscreteModel("../../models/Step_index_fiber.msh");
permittivity=x->1+2*(x[1]^2+x[2]^2<=1);
m=FEM(1,5,permittivity,model,sqrt(3),order=2,solver=:MUMPS,field=true,type=:Vector);
fig,ax,plot_obj=wireframe(m[1].Ω, color=:black, linewidth=2,axis=(aspect=DataAspect(),))
ax.xlabel="x (µm)";
ax.ylabel="y (µm)";
save("FEM1_mesh.png",fig); nothing #hide
```
![Pz for FM computed with FEM](FEM1_mesh.png)

Now we can plot the z-component of the Poynting Vector with the package GridapMakie.jl (it is also possible to save the field in a vtk file with the funcion `writevtk` implemented in Gridap.jl and open it with ParaView):
```@example 3
Px,Py,Pz=PoyntingVector(m[1]);
fig,ax,plot_obj=GLMakie.plot(m[1].Ω,Pz,axis=(aspect=DataAspect(),),colormap=:jet)
ax.xlabel="x (µm)";
ax.ylabel="y (µm)";
ax.title="neff = $(real(m[1].neff))";
Colorbar(fig[1,2], plot_obj);
save("FEM1_Pz.png",fig); nothing #hide
```
![Pz for FM computed with FEM](FEM1_Pz.png)

## Leaky modes in step index fiber with a low index trench
In this example, we study a fiber with the parameters below:
- core with a radius of 4 µm and a refractive index of 1.46
- a trench located between 4 µm and 7 µm with a refractive index of 1.41
- a cladding with a refractive index of 1.45

In order to compute the leaky scalar modes at $\lambda=1.6$ µm, a PML must be added when using the 1-dimensional FD solver. The PML will be located between 12 µm and 15 µm.
```@example 5
using OpticalFibers
using OpticalFibers.ModeSolvers
f=r->1.46-0.05*(r<=7)*(r>=4)-0.01*(r>7);
m0=FD(1.6,0,10,f,1000,15,field=true,rPML=12)
```
The LP01 mode is the mode 1. Since its effective index is greater than the refractive index of the cladding, this mode is guided and its losses are not significant.

The LP02 mode is the mode 7, its losses in dB/km can be calculated:
```@example 5
losses(m0[7])*1e6
```

```@example 5
using Plots
plot(m0[1].r,-real.(m0[1].E),xlabel="r (µm)",ylabel="real(E)",label="LP01")
plot!(m0[7].r,-real.(m0[7].E),label="LP02")
```

We can also compute the effective index of the LP11 mode:
```@example 5
m1=FD(1.6,1,10,f,1000,15,field=true,rPML=12)
pos=argmin(losses.(m1)*1E6);
m1[pos].neff
```

To compute the leaky vector modes that correspond to the LP11 mode, we can add a PML for the FEM solver:
```@example 5
using Gridap
using GridapGmsh
model = GmshDiscreteModel("../../models/Step_index_fiber_pml.msh");
permittivity=x->f(hypot(x[1],x[2]))^2;
epsilon=add_cylindrical_PML(permittivity,12,3,10);
mu=add_cylindrical_PML(x->1.0,12,3,10);
m=FEM(1.6,4,epsilon,mu,model,real(m1[pos].neff),field=true,order=2,solver=:MUMPS)
```

The field can be plotted using GridapMakie:
```@example 5
using GridapMakie
using GLMakie
fig,ax,plot_obj=GLMakie.plot(m[end].Ω,real(m[1].Ex),axis=(aspect=DataAspect(),),colormap=:jet)
Colorbar(fig[1,2], plot_obj);
save("FEM_PML_Ex.png",fig); nothing #hide
```
![Real(Ex) for the first mode](FEM_PML_Ex.png)

## Elliptical core fiber with a low index trench
In this example, the fiber is similar to the previous one but the core and the trench are elliptic:
- core with a radius of 4 µm in the x direction and 2 µm in the y direction and a refractive index of 1.46
- a trench with a refractive index of 1.41 and located between the core and an ellipse with radii 8 and 4 µm in x and y direction respectively
- a cladding with a refractive index of 1.45

![Elliptic Fiber RIP](./assets/Elliptic_profile.png)

A rectangular PML is added for $12<\vert x \vert < 15$ µm and $8<\vert x \vert < 11$ µm to compute the leaky modes.

The first mode is saved in a file that can be opened with ParaView.

```@example 6
using OpticalFibers
using OpticalFibers.ModeSolvers
f=x->1.46-0.05*(x[1]^2/16+x[2]^2/4>=1)+0.04*(x[1]^2/64+x[2]^2/16>1);
epsilon=add_rectangular_PML(x->f(x)^2,12,3,8,3,10);
mu=add_rectangular_PML(x->1,12,3,8,3,10);
m=FEM(1.6,4,epsilon,mu,model,1.44,field=true,order=2,solver=:MUMPS);
writevtk("mode1.vtu",m[1]); nothing #hide
```

![Elliptic Fiber Ex](./assets/Elliptic_Ex.png)

## FEM: Photonic Crystal Fiber
In a PCF, the modes are not guided modes but leaky modes so that the computation requires a PML. The fiber is constituted of three rings of air hole (n=1) inserted in silica (n=1.45). The pitch is 2 µm, the hole diameter is 1.5 µm and the PML begins at 8 µm from the fiber center and its thickness is 2 µm. 

![PCF mesh](./assets/PCF.png)

First, the mesh is loaded:
```@example 4
using OpticalFibers
using OpticalFibers.ModeSolvers
using Gridap
using GridapGmsh
using GridapMakie
using GLMakie
model = GmshDiscreteModel("../../models/PCF.msh");
```
Then we define the permittivity function and add a PML to obtain the final permittivity and permeability tensors:
```@example 4
function eps_PCF(x)
    Pitch=2;
    r_hole=0.75;
    x1,y1=ring(1);
    x2,y2=ring(2);
    x3,y3=ring(3);
    xc=[x1;x2;x3]*Pitch;
    yc=[y1;y2;y3]*Pitch;
    for i=1:length(xc)
        rc=hypot(x[1]-xc[i],x[2]-yc[i]);
        if rc<r_hole
            return 1.0;
        end
    end
    return 1.45^2;
end
epsilon=add_cylindrical_PML(eps_PCF,8,2,10);
mu=add_cylindrical_PML(x->1.0,8,2,10);
```
Then we can compute ten modes whose effective indices are close to the approximate value calculated for the fundamental mode at the wavelength of 1.3 µm:
```@example 4
neff_approx=approx_neff_PCF(1.3,1.5,2);
m=FEM(1.3,10,epsilon,mu,model,neff_approx,field=true,order=2,solver=:MUMPS)
```
The last two modes are fundamental modes. We can compute and plot the z-component of the Poynting vector of the last mode:
```@example 4
Px,Py,Pz=PoyntingVector(m[end]);
fig,ax,plot_obj=GLMakie.plot(m[end].Ω,Pz,axis=(aspect=DataAspect(),),colormap=:jet)
ax.xlabel="x (µm)";
ax.ylabel="y (µm)";
ax.title="neff = $(m[end].neff)";
Colorbar(fig[1,2], plot_obj);
save("FEM2_Pz.png",fig); nothing #hide
```
![Pz for FM computed with FEM](FEM2_Pz.png)
