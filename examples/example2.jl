#Gradient Index Fiber and Dispersion
using OpticalFibers
using OpticalFibers.PhysicalData
using Plots
lambda=1:0.01:1.5;
epsilon=[x->(n_Ge_Doped_Silica(l*1e-6,0)+(n_Ge_Doped_Silica(l*1e-6,0.2)-n_Ge_Doped_Silica(l*1e-6,0))*(x[1]<=3.5)*(1-x[1]^2/3.5^2))^2 for l in lambda];
r=0:0.1:5;
plot(r,sqrt.(epsilon[1].(r)),label="λ = 1 µm",xlabel="r (µm)",ylabel="Refractive index");
plot!(r,sqrt.(epsilon[end].(r)),label="λ = 1.5 µm")

using OpticalFibers.ModeSolvers
using Gridap
model = CartesianDiscreteModel((0,20),2000)
m=FEM1D.(lambda,[0,1,2]',epsilon,Ref(model),neigs=2);
#mode LP01 always exists
neff01=[real(m[j,1][1].neff) for j in 1:length(lambda)];
N02=sum((length.(m[:,1])).>=2);
neff02=[real(m[j,1][2].neff) for j in 1:N02];
N11=sum((length.(m[:,2])).>=1);
neff11=[real(m[j,2][1].neff) for j in 1:N11];
N21=sum((length.(m[:,3])).>=1);
neff21=[real(m[j,3][1].neff) for j in 1:N21];
plot(lambda,n_Ge_Doped_Silica.(lambda*1E-6,0),label="Silica",xlabel="λ (µm)",ylabel="Effective index",color=:black);
plot!(lambda,n_Ge_Doped_Silica.(lambda*1E-6,0.2),label="Ge-doped Silica (20%)",color=:black,line=:dash);
plot!([lambda,lambda[1:N11],lambda[1:N02],lambda[1:N21]],[neff01,neff11,neff02,neff21],label=["LP01" "LP11" "LP02" "LP21"])

beta01=neff01*2*pi./lambda*1E6;
omega=2*pi*OpticalFibers.PhysicalData.c./lambda*1E6;
omega2,beta2=derivative((omega,beta01),2);
lambda2=2*pi*OpticalFibers.PhysicalData.c./omega2*1E6;
plot(lambda2,beta2*1E26,xlabel="λ (µm)",ylabel="β₂ (10⁻²⁶ s²/m)")

m=FEM1D.(lambda,0,epsilon,Ref(model),field=true);
m=[m[i][1] for i in 1:length(lambda)];
A=Aeff.(m);
gamma=nonLinearCoefficient.(m,2.53E-20)*1E21;
plot(lambda,A,label="Aeff",ylabel="Effective area (µm²)",xlabel="λ (µm)",color=:blue,leg=:topright);
plot!(twinx(),lambda,gamma,label="γ",ylabel="Non-linear coefficient ((W.km)⁻¹)",color=:red,leg=:topleft)