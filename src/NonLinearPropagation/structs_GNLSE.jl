@kwdef mutable struct GNLSE_param
    N::Int #Number of points, must be a power of 2
    τmax::Real #Maximum time
    λ0::Real #Central Wavelength
    dz::Real=1E-3 #Step of the R4KIP method to solve the GNLSE
    precision::Real=0 #Precision of the method. If 0, constant step.
    precision_type::Symbol=:Global #type of error (useless if precision=0)
    MGNLSE::Bool=false #True for using the dispersion of the effective area (Modified GNLSE)
end

Base.show(io::IO, p::GNLSE_param) = print(io,"N = ",p.N,"\nτₘₐₓ = ",p.τmax,"\ndτ = ",get_dτ(p),"\nτ = ",get_τmin(p)," : ",get_τmax(p),"\nλ₀ = ",p.λ0,"\nω₀ = ",get_ω0(p),"\ndΩ = ",get_dΩ(p),"\nΩ = ",get_Ωmin(p)," : ",get_Ωmax(p),"\nω = ",get_ωmin(p)," : ",get_ωmax(p),"\nλ = ",get_λmin(p)," : ",get_λmax(p),"\ndz = ",p.dz,"\nprecision = ",p.precision,"\nprecision type = ",p.precision_type,"\nMGNLSE = ",p.MGNLSE);

@kwdef mutable struct GNLSE_mode
    α::Union{Real,Vector{<:Real},Tuple{Function,Symbol}}=0 #Losses (m⁻¹)
    β::Union{Real,Vector{<:Real},Tuple{Function,Symbol}} #β-β₀-β₁Ω (m⁻¹)
    γ::Union{Real,Vector{<:Real},Tuple{Function,Symbol}} #Nonlinear coefficient (W⁻¹m⁻¹)
    Aeff::Union{Real,Vector{<:Real},Tuple{Function,Symbol}}=1 #Effective Area (m²), usefull if MGLNSE=true in parameters
end

@kwdef mutable struct GNLSE_Raman
    fR::Real
    hR::Tuple{Function,Symbol}
end

function get_ω0(p::GNLSE_param)
    return 2*pi*c/p.λ0;
end

function get_dτ(p::GNLSE_param)
    return 2*p.τmax/p.N;
end

function get_τ(p::GNLSE_param)
    return get_dτ(p)*collect(-div(p.N,2):1:(div(p.N,2)-1));
end

function get_τmin(p::GNLSE_param)
    return -get_dτ(p)*div(p.N,2);
end

function get_τmax(p::GNLSE_param)
    return get_dτ(p)*(div(p.N,2)-1);
end

function get_dΩ(p::GNLSE_param)
    return pi/p.τmax;
end

function get_Ω(p::GNLSE_param)
    return get_dΩ(p)*collect(-div(p.N,2):1:(div(p.N,2)-1));
end

function get_Ωmin(p::GNLSE_param)
    return -get_dΩ(p)*div(p.N,2);
end

function get_Ωmax(p::GNLSE_param)
    return get_dΩ(p)*(div(p.N,2)-1);
end

function get_ω(p::GNLSE_param)
    return get_Ω(p).+get_ω0(p);
end

function get_ωmin(p::GNLSE_param)
    return get_Ωmin(p).+get_ω0(p);
end

function get_ωmax(p::GNLSE_param)
    return get_Ωmax(p).+get_ω0(p);
end

function get_λ(p::GNLSE_param)
    return 2*pi*c./get_ω(p);
end

function get_λmin(p::GNLSE_param)
    return 2*pi*c/get_ωmax(p);
end

function get_λmax(p::GNLSE_param)
    return 2*pi*c/get_ωmin(p);
end

function verify(p::GNLSE_param)
    ok=true;
    msg="ok";
    if p.N<32
        ok=false;
        msg="N must be greater than 32";
        if (p.N & (p.N-1))!=0
            ok=false;
            msg="N must be a power of 2";
        end
    end
    if p.τmax<=0
        ok=false;
        msg="τmax must be greater than 0";
    end
    if p.λ0<=0
        ok=false;
        msg="λ0 must be greater than 0";
    end
    if p.dz<=0
        ok=false;
        msg="dz must be greater than 0";
    end
    if p.precision<0
        ok=false;
        msg="precision must be postive or null";
    end
    return ok,msg;
end

function get_α(m::GNLSE_mode,p::GNLSE_param)
    if isa(m.α,Union{Real,Vector{<:Real}})
        if (length(m.α)==1)
            alpha=m.α*ones(p.N);
        else
            Omega=get_Ω(p);
            alpha=zeros(p.N);
            for i=1:length(m.α)
                @. alpha=alpha+m.α[i]*Omega^(i-1)/factorial(i-1);
            end
        end
    else
        if m.α[2]==:lambda
            lambda=get_λ(p);
            alpha=m.α[1].(lambda);
        elseif m.α[2]==:Omega
            Omega=get_Ω(p);
            alpha=m.α[1].(Omega);
        elseif m.α[2]==:omega
            omega=get_ω(p);
            alpha=m.α[1].(omega);
        else
            alpha=[];
        end 
    end
    return alpha;
end

function get_γ(m::GNLSE_mode,p::GNLSE_param)
    if isa(m.γ,Union{Real,Vector{<:Real}})
        if (length(m.γ)==1)
            gamma=m.γ*ones(p.N);
        else
            Omega=get_Ω(p);
            gamma=zeros(p.N);
            for i=1:length(m.γ)
                @. gamma=gamma+m.γ[i]*Omega^(i-1)/factorial(i-1);
            end
        end
    else
        if m.γ[2]==:lambda
            lambda=get_λ(p);
            gamma=m.γ[1].(lambda);
        elseif m.γ[2]==:Omega
            Omega=get_Ω(p);
            gamma=m.γ[1].(Omega);
        elseif m.γ[2]==:omega
            omega=get_ω(p);
            gamma=m.γ[1].(omega);
        else
            gamma=[];
        end 
    end
    return gamma;
end

function get_Aeff(m::GNLSE_mode,p::GNLSE_param)
    if isa(m.Aeff,Union{Real,Vector{<:Real}})
        if (length(m.Aeff)==1)
            Aeff=m.Aeff*ones(p.N);
        else
            Omega=get_Ω(p);
            Aeff=zeros(p.N);
            for i=1:length(m.Aeff)
                @. Aeff=Aeff+m.Aeff[i]*Omega^(i-1)/factorial(i-1);
            end
        end
    else
        if m.Aeff[2]==:lambda
            lambda=get_λ(p);
            Aeff=m.Aeff[1].(lambda);
        elseif m.Aeff[2]==:Omega
            Omega=get_Ω(p);
            Aeff=m.Aeff[1].(Omega);
        elseif m.Aeff[2]==:omega
            omega=get_ω(p);
            Aeff=m.Aeff[1].(omega);
        else
            Aeff=[];
        end 
    end
    return Aeff;
end

function get_β(m::GNLSE_mode,p::GNLSE_param)
    if isa(m.β,Union{Real,Vector{<:Real}})
        Omega=get_Ω(p);
        beta=zeros(p.N);
        for i=1:length(m.β)
            @. beta=beta+m.β[i]*Omega^(i+1)/factorial(i+1);
        end
    else
        if m.β[2]==:lambda
            lambda=get_λ(p);
            beta=m.β[1].(lambda);
        elseif m.β[2]==:Omega
            Omega=get_Ω(p);
            beta=m.β[1].(Omega);
        elseif m.β[2]==:omega
            omega=get_ω(p);
            beta=m.β[1].(omega);
        else
            beta=[];
        end 
    end
    return beta;
end

function get_hR(r::GNLSE_Raman,p::GNLSE_param)
    if (r.fR==0)
        return Vector{ComplexF64}(undef,0);
    else
        if (r.hR[2]==:Omega)
            hR=(r.hR[1].(get_Ω(p)))/real(r.hR[1](0));
        else
            tau=get_τ(p).+p.τmax;
            hRw=ifft(r.hR[1].(tau));
            hR=hRw/real(hRw[1]);
        end
        return hR;
    end
end

function hr(t::Float64,type::Integer)
    if type==1
        c=299792458.0;
        integrale=0.375E-11;
        CP=[5625,10000,23125,36250,46300,49700,61150,69167,79367,83550,93000,108000,121500];
        Ai=[1,11.4,36.67,67.67,74,4.5,6.8,4.6,4.2,4.5,2.7,3.1,3.0];
        GFWHM=[5210,11042,17500,16250,13533,2450,4150,15500,5950,6430,15000,9100,16000];
        LFWHM=[1737,3881,5833,5417,4511,817,1383,5167,1983,2143,5000,3033,5333];
        output=0
        for i=1:13
            output=output+Ai[i]*exp.(-pi*c*LFWHM[i]*t).*exp.(-(pi*c*GFWHM[i]*t/2.0).^2).*sin.(2.0*pi*c*CP[i]*t);
        end
        output=output/integrale.*(t.>0);
    else
        tau1=12.2E-15;
        tau2=32.0E-15;
        output = (tau1^2+tau2^2)/(tau1*tau2^2)*exp.(-t/tau2).*sin.(t/tau1).*(t.>0);
    end
end
            



