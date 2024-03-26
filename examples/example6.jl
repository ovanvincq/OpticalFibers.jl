#Photonic Crystal Fiber (PCF)
using OpticalFibers
using OpticalFibers.ModeSolvers
using Gridap
using GridapGmsh
using GridapMakie
using GLMakie
model = GmshDiscreteModel("./models/PCF.msh");

function eps_PCF(x)
    Pitch=2;
    r_hole=0.75;
    x1,y1=ring(1);
    x2,y2=ring(2);
    x3,y3=ring(3);
    xc=[x1;x2;x3]*Pitch;
    yc=[y1;y2;y3]*Pitch;
    for i in axes(xc,1)
        rc=hypot(x[1]-xc[i],x[2]-yc[i]);
        if rc<r_hole
            return 1.0;
        end
    end
    return 1.45^2;
end

neff_approx=approx_neff_PCF(1.3,1.5,2);
m=FEM2D(1.3,eps_PCF,model,neigs=4,approx_neff=neff_approx,field=true,solver=:MUMPS,type=:Vector,dPML=2)

Px,Py,Pz=PoyntingVector(m[end]);
fig,ax,plot_obj=GLMakie.plot(get_triangulation(Pz),Pz,axis=(aspect=DataAspect(),),colormap=:jet)
ax.xlabel="x (µm)";
ax.ylabel="y (µm)";
ax.title="neff = $(m[end].neff)";
Colorbar(fig[1,2], plot_obj);
