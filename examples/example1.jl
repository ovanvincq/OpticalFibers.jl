#Bimodal Step-Index Fiber

#Scalar Modes
using OpticalFibers
using OpticalFibers.ModeSolvers
m0=multi_step_fiber_modes(1,0,2,[1.47,1.45],field=true);
m01=m0[1]
m1=multi_step_fiber_modes(1,1,2,[1.47,1.45],field=true);
m11=m1[1]

m=multi_step_fiber_modes.(1,[0,1],2,Ref([1.47,1.45]),field=true);
m01=m[1][1];
m11=m[2][1];

using Plots
r=0:0.01:10;
plot(r,computeField(m01,r),label=m01.Name)
plot!(r,computeField(m11,r),label=m11.Name)

mm01=Mode{ScalarFieldFunction2D}(m01)
mm11c=Mode{ScalarFieldFunction2D}(m11)
mm11s=Mode{ScalarFieldFunction2D}(m11,90)
normalize!(mm01)
normalize!(mm11c)
normalize!(mm11s)
x=-8:0.125:8;
contourf(x,x,computeField(mm11c,x,x')',levels=100,linewidth=0)
L=1/(m01.neff-m11.neff)

anim=@animate for j=0:214
    TotalField=getField(mm01,j)+getField(mm11c,j);
    contourf(x,x,abs2.(computeField(TotalField,x,x'))',levels=100,linewidth=0,title="z = $j µm")
    #title!("z = $j µm");
end;
gif(anim,"anim_field.gif",fps=15)

#Vector Modes
mv0=multi_step_fiber_modes(1,0,2,[1.47,1.45],field=true,type=:Vector)
mv1=multi_step_fiber_modes(1,1,2,[1.47,1.45],field=true,type=:Vector)
mv2=multi_step_fiber_modes(1,2,2,[1.47,1.45],field=true,type=:Vector)

Px,Py,Pz=PoyntingVector(mv2[1]);
contourf(x,x,Pz.(tuple.(x,x'))',linewidth=0,levels=100,xlims=(-4,4),ylims=(-4,4))
X,Y=meshgrid(x[1:4:end]);
quiver!(X,Y,quiver=(mv2[1].field.Ex.(tuple.(x[1:4:end],x[1:4:end]'))'/20,mv2[1].field.Ey.(tuple.(x[1:4:end],x[1:4:end]'))'/20),color=:cyan,arrow=arrow(:closed))

mv=[mv1;mv0;mv2];
normalize!.(mv,atol=1E-5);
overlap.(mv,transpose(mv),atol=1E-15,rtol=1E-5)

#Scalar Modes with FEM1D
using Gridap
model = CartesianDiscreteModel((0,15),1500)
epsilon=x->(1.47-0.02*(x[1]>=2))^2
m0FEM=FEM1D(1,0,epsilon,model,field=true,neigs=4)
normalize!(m0FEM[1])
normalize!(m0[1])
r=0:0.01:10
plot(r,computeField(m0[1],r),label="Step-index solver")
plot!(r,-real(computeField(m0FEM[1],r)),label="FEM1D solver",line=:dash)

#Scalar Modes with FEM2D
using GridapGmsh
model = GmshDiscreteModel("./models/example1.msh");
epsilon2D=x->(1.47-0.02*(hypot(x[1],x[2])>=2))^2
mFEM2D=FEM2D(1,epsilon2D,model,field=true,neigs=4)
normalize!.(mFEM2D)
p1=contourf(x,x,real(computeField(mFEM2D[1],x,x'))',levels=50,linewidth=0,aspect_ratio=:equal);
p2=contourf(x,x,real(computeField(mFEM2D[2],x,x'))',levels=50,linewidth=0,aspect_ratio=:equal);
p3=contourf(x,x,real(computeField(mFEM2D[3],x,x'))',levels=50,linewidth=0,aspect_ratio=:equal);
plot(p1, p2, p3, layout=(1,3), legend=false,size=(900,300))
abs2.(overlap.(mFEM2D,[mm01 mm11c mm11s]))

#Vector Modes with FEM2D
mvFEM2D=FEM2D(1,epsilon2D,model,field=true,neigs=8,type=:Vector)
normalize!.(mvFEM2D)
abs2.(overlap.(mvFEM2D,transpose(mv)))
p1=contourf(x,x,real(computeField(mvFEM2D[4],x,x',:Ex))',levels=50,linewidth=0,aspect_ratio=:equal,title="Ex");
p2=contourf(x,x,real(computeField(mvFEM2D[4],x,x',:Ey))',levels=50,linewidth=0,aspect_ratio=:equal,title="Ey");
p3=contourf(x,x,imag(computeField(mvFEM2D[4],x,x',:Ez))',levels=50,linewidth=0,aspect_ratio=:equal,title="Ez");
p4=contourf(x,x,real(computeField(mvFEM2D[4],x,x',:Hx))',levels=50,linewidth=0,aspect_ratio=:equal,title="Hx");
p5=contourf(x,x,real(computeField(mvFEM2D[4],x,x',:Hy))',levels=50,linewidth=0,aspect_ratio=:equal,title="Hy");
p6=contourf(x,x,imag(computeField(mvFEM2D[4],x,x',:Hz))',levels=50,linewidth=0,aspect_ratio=:equal,title="Hz");
plot(p1, p2, p3, p4, p5, p6, layout=(2,3), legend=false,size=(900,600))

p1=contourf(x,x,real(computeField(mv[4],x,x',:Ex))',levels=50,linewidth=0,aspect_ratio=:equal,title="Ex");
p2=contourf(x,x,real(computeField(mv[4],x,x',:Ey))',levels=50,linewidth=0,aspect_ratio=:equal,title="Ey");
p3=contourf(x,x,imag(computeField(mv[4],x,x',:Ez))',levels=50,linewidth=0,aspect_ratio=:equal,title="Ez");
p4=contourf(x,x,real(computeField(mv[4],x,x',:Hx))',levels=50,linewidth=0,aspect_ratio=:equal,title="Hx");
p5=contourf(x,x,real(computeField(mv[4],x,x',:Hy))',levels=50,linewidth=0,aspect_ratio=:equal,title="Hy");
p6=contourf(x,x,imag(computeField(mv[4],x,x',:Hz))',levels=50,linewidth=0,aspect_ratio=:equal,title="Hz");
plot(p1, p2, p3, p4, p5, p6, layout=(2,3), legend=false,size=(900,600))
