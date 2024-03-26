using FFTW
using LinearAlgebra
using ProgressMeter

struct GNLSE_workspace
    dz::Float64
    precision::Float64
    precision_type::Symbol
    gini::Vector{ComplexF64}
    iba::Vector{ComplexF64}
    tmp::Vector{ComplexF64}
    k1::Vector{ComplexF64}
    k2::Vector{ComplexF64}
    k3::Vector{ComplexF64}
    k4::Vector{ComplexF64}
    Ai::Vector{ComplexF64}
    g::Vector{ComplexF64}
    D::Vector{ComplexF64}
    tmp2::Vector{ComplexF64}
    tmp3::Vector{ComplexF64}
    tmp_abs::Vector{ComplexF64}
    pfft::FFTW.cFFTWPlan{ComplexF64, -1, false, 1, Tuple{Int64}};
    fR::Float64
    hR::Vector{ComplexF64}
end  

#=@kwdef mutable struct test2
    a :: Float64 =1
    test2(a) = a < 0 ? error("out of order") : new(a)
end 

Core.setfield!(f::test2,name::Symbol,x::Real) = x<0 ? error("erreur") : setfield!(f,name,x);=#

function GNLSE_workspace(p::GNLSE_param,m::GNLSE_mode,r::GNLSE_Raman)
    dz=Float64(p.dz);
    precision=Float64(p.precision);
    precision_type=p.precision_type;
    gini=fftshift(im*get_γ(m,p));
    iba=fftshift(im*get_β(m,p)-get_α(m,p)/2.0)/2.0;
    tmp=Vector{ComplexF64}(undef,p.N);
    k1=Vector{ComplexF64}(undef,p.N);
    k2=Vector{ComplexF64}(undef,p.N);
    k3=Vector{ComplexF64}(undef,p.N);
    k4=Vector{ComplexF64}(undef,p.N);
    Ai=Vector{ComplexF64}(undef,p.N);
    g=Vector{ComplexF64}(undef,p.N);
    D=Vector{ComplexF64}(undef,p.N);
    tmp2=Vector{ComplexF64}(undef,p.N);
    if (r.fR==0)
        tmp3=Vector{ComplexF64}(undef,0);
        tmp_abs=Vector{ComplexF64}(undef,0);
    else
        tmp3=Vector{ComplexF64}(undef,p.N);
        tmp_abs=Vector{ComplexF64}(undef,p.N);
    end
    pfft=plan_fft(tmp);
    fR=r.fR;
    hR=get_hR(r,p);
    return GNLSE_workspace(dz,precision,precision_type,gini,iba,tmp,k1,k2,k3,k4,Ai,g,D,tmp2,tmp3,tmp_abs,pfft,fR,hR)
end

function step_without_Raman!(vv::Vector{ComplexF64},w::GNLSE_workspace)
        mul_perso!(w.Ai,w.D,vv);
    #Compute k1
        mul!(w.tmp2,w.pfft,vv)
        mul3!(w.tmp2);
        mul!(w.tmp,inv(w.pfft),w.tmp2);
    #Compute k2
        computeK1!(w.k1,w.tmp,w.D,w.g,w.Ai);
        mul!(w.tmp2,w.pfft,w.tmp)
        mul3!(w.tmp2);
        mul!(w.tmp,inv(w.pfft),w.tmp2);
    #Compute k3
        computeK2!(w.k2,w.tmp,w.g,w.Ai);
        mul!(w.tmp2,w.pfft,w.tmp)
        mul3!(w.tmp2);
        mul!(w.tmp,inv(w.pfft),w.tmp2);
    #compute k4
        computeK3!(w.k3,w.tmp,w.D,w.g,w.Ai);
        mul!(w.tmp2,w.pfft,w.tmp)
        mul3!(w.tmp2);
        mul!(w.tmp,inv(w.pfft),w.tmp2);
        kernel_final!(w.Ai,w.k1,w.k2,w.k3,w.g,w.D,w.tmp,vv)
end

function step_with_Raman!(vv::Vector{ComplexF64},w::GNLSE_workspace)
    mul_perso!(w.Ai,w.D,vv);
#Compute k1
    mul!(w.tmp,w.pfft,vv)
    compute_abs2!(w.tmp,w.tmp_abs);
    mul!(w.tmp2,inv(w.pfft),w.tmp_abs);
    mul_perso!(w.tmp2,w.hR,w.tmp2);
    mul!(w.tmp3,w.pfft,w.tmp2);
    mulhr2!(w.tmp,w.tmp3,w.tmp_abs,w.fR);
    mul!(w.tmp2,inv(w.pfft),w.tmp3);
#Compute k2
    computeK1!(w.k1,w.tmp2,w.D,w.g,w.Ai);
    mul!(w.tmp,w.pfft,w.tmp2);
    compute_abs2!(w.tmp,w.tmp_abs);
    mul!(w.tmp2,inv(w.pfft),w.tmp_abs);
    mul_perso!(w.tmp2,w.hR,w.tmp2);
    mul!(w.tmp3,w.pfft,w.tmp2);
    mulhr2!(w.tmp,w.tmp3,w.tmp_abs,w.fR);
    mul!(w.tmp2,inv(w.pfft),w.tmp3);
#Compute k3
    computeK2!(w.k2,w.tmp2,w.g,w.Ai);
    mul!(w.tmp,w.pfft,w.tmp2);
    compute_abs2!(w.tmp,w.tmp_abs);
    mul!(w.tmp2,inv(w.pfft),w.tmp_abs);
    mul_perso!(w.tmp2,w.hR,w.tmp2);
    mul!(w.tmp3,w.pfft,w.tmp2);
    mulhr2!(w.tmp,w.tmp3,w.tmp_abs,w.fR);
    mul!(w.tmp2,inv(w.pfft),w.tmp3);
#compute k4
    computeK3!(w.k3,w.tmp2,w.D,w.g,w.Ai);
    mul!(w.tmp,w.pfft,w.tmp2);
    compute_abs2!(w.tmp,w.tmp_abs);
    mul!(w.tmp2,inv(w.pfft),w.tmp_abs);
    mul_perso!(w.tmp2,w.hR,w.tmp2);
    mul!(w.tmp3,w.pfft,w.tmp2);
    mulhr2!(w.tmp,w.tmp3,w.tmp_abs,w.fR);
    mul!(w.tmp2,inv(w.pfft),w.tmp3);
    kernel_final!(w.Ai,w.k1,w.k2,w.k3,w.g,w.D,w.tmp2,vv);
end
 
function computeK1!(k1::Vector{ComplexF64},tmp::Vector{ComplexF64},D::Vector{ComplexF64},gamma::Vector{ComplexF64},Ai::Vector{ComplexF64})
    Threads.@threads for i in eachindex(k1)
        @inbounds k1[i] = D[i]*tmp[i]*gamma[i]/2.0;
        @inbounds tmp[i]=Ai[i]+k1[i];
    end
end

function computeK2!(k2::Vector{ComplexF64},tmp::Vector{ComplexF64},gamma::Vector{ComplexF64},Ai::Vector{ComplexF64})
    Threads.@threads for i in eachindex(k2)
        @inbounds k2[i] = tmp[i]*gamma[i];
        @inbounds tmp[i]=Ai[i]+k2[i]/2.0;
    end
end

function computeK3!(k3::Vector{ComplexF64},tmp::Vector{ComplexF64},D::Vector{ComplexF64},gamma::Vector{ComplexF64},Ai::Vector{ComplexF64})
    Threads.@threads for i in eachindex(k3)
        @inbounds k3[i] = tmp[i]*gamma[i];
        @inbounds tmp[i]=(Ai[i]+k3[i])*D[i];
    end
end

function mul3!(tmp::Vector{ComplexF64})
    Threads.@threads for i in eachindex(tmp)
        @inbounds tmp[i] = abs2(tmp[i])*tmp[i];
    end
end

function mul_perso!(Ai::Vector{ComplexF64},D::Vector{ComplexF64},vv::Vector{ComplexF64})
    Threads.@threads for i in eachindex(Ai)
        @inbounds Ai[i] = D[i]*vv[i];
    end
end

function kernel_final!(Ai::Vector{ComplexF64},k1::Vector{ComplexF64},k2::Vector{ComplexF64},k3::Vector{ComplexF64},gamma::Vector{ComplexF64},D::Vector{ComplexF64},tmp::Vector{ComplexF64},vv::Vector{ComplexF64})
    Threads.@threads for i in eachindex(Ai)
        @inbounds vv[i] = D[i]*((k1[i]+k2[i]+k3[i])/3.0+Ai[i])+gamma[i]*tmp[i]/6.0;
    end
end

function compute_abs2!(tmp::Vector{ComplexF64},tmp_abs::Vector{ComplexF64})
    Threads.@threads for i in eachindex(tmp)
        @inbounds tmp_abs[i] = abs2(tmp[i]);
    end
end

function mulhr2!(tmp::Vector{ComplexF64},tmp3::Vector{ComplexF64},tmp_abs::Vector{ComplexF64},fR::Float64)
    Threads.@threads for i in eachindex(tmp)
        @inbounds tmp3[i] = tmp[i]*real((1.0-fR)*tmp_abs[i]+fR*tmp3[i]);
    end
end

function propagation(vv::Vector{ComplexF64},L::Real,p::GNLSE_param,m::GNLSE_mode,r::GNLSE_Raman=GNLSE_Raman(0,(x->1,:Omega)),ProgressBar::Bool=false)
    ok,msg=verify(p);
    if (ok)
        w=GNLSE_workspace(p,m,r)
        vvf=copy(vv);
        if (w.precision==0)
            propagation_constant_step!(w,vvf,Float64(L),ProgressBar);
        end
        return vvf;
    else
        error("Error in GNLSE_param: "*msg)
    end
end

function compute_Dg!(w,dz)
    Threads.@threads for i in eachindex(w.D)
        @inbounds w.g[i] = w.gini[i]*dz;
        @inbounds w.D[i] = exp(w.iba[i]*dz);
    end
end

function propagation_constant_step!(w::GNLSE_workspace,vvf::Vector{ComplexF64},L::Float64,ProgressBar::Bool=false)
    nb_step=floor(Int,L/w.dz);
    dz_final=L-nb_step*w.dz;
    compute_Dg!(w,w.dz);
    if (ProgressBar)
        if (w.fR==0)
            p = Progress(nb_step)
            for i=1:nb_step
                step_without_Raman!(vvf,w);
                next!(p)
            end
            if (dz_final>0)
                computeDg!(w,dz_final);
                step_without_Raman!(vvf,w);
            end
            finish!(p)
        else
            p = Progress(nb_step)
            for i=1:nb_step
                step_with_Raman!(vvf,w);
                next!(p)
            end
            if (dz_final>0)
                computeDg!(w,dz_final);
                step_with_Raman!(vvf,w);
            end
            finish!(p)
        end
    else
        if (w.fR==0)
            for i=1:nb_step
                step_without_Raman!(vvf,w);
            end
            if (dz_final>0)
                computeDg!(w,dz_final);
                step_without_Raman!(vvf,w);
            end
        else
            for i=1:nb_step
                step_with_Raman!(vvf,w);
            end
            if (dz_final>0)
                computeDg!(w,dz_final);
                step_with_Raman!(vvf,w);
            end
        end
    end
end