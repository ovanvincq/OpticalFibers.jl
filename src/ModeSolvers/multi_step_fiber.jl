function besselj_prime(nu,x)
    return 0.5*(besselj(nu-1,x)-besselj(nu+1,x));
end

function besseli_prime(nu,x)
    return 0.5*(besseli(nu-1,x)+besseli(nu+1,x));
end

function besselk_prime(nu,x)
    return -0.5*(besselk(nu-1,x)+besselk(nu+1,x));
end

function bessely_prime(nu,x)
    return 0.5*(bessely(nu-1,x)-bessely(nu+1,x));
end

function phi1(index::Float64,radius::Float64,neff::Vector{Float64},r::Float64,k0::Float64,nu::Int64)
    U=sqrt.(abs.((neff.^2).-(index^2)))*k0*radius;
    U1=U[neff.<index];
    U2=U[neff.==index];
    U3=U[neff.>index];
    y1=besselj.(nu,U1*r/radius)./besselj.(nu,U1);
    y2=(r/radius)^nu*ones(size(U2));
    y3=besseli.(nu,U3*r/radius)./besseli.(nu,U3);
    return [y1;y2;y3]
end

function phi1(index::Float64,radius::Float64,neff::Float64,r::Union{Float64,Vector{Float64}},k0::Float64,nu::Int64)
    U=sqrt(abs((neff^2)-(index^2)))*k0*radius;
    if (neff<index)
        return besselj.(nu,U.*r./radius)./besselj(nu,U);
    end
    if (neff==index)
        return (r./radius).^nu;
    end
    return besseli.(nu,U.*r./radius)./besseli(nu,U);
end

function phi1_prime(index::Float64,radius::Float64,neff::Vector{Float64},r::Float64,k0::Float64,nu::Int64)
    U=sqrt.(abs.((neff.^2).-(index^2)))*k0*radius;
    U1=U[neff.<index];
    U2=U[neff.==index];
    U3=U[neff.>index];
    y1=besselj_prime.(nu,U1*r/radius)./besselj.(nu,U1).*U1/radius;
    y2=(r/radius)^(nu-1)*nu/radius*ones(size(U2));
    y3=besseli_prime.(nu,U3*r/radius)./besseli.(nu,U3).*U3/radius;
    return [y1;y2;y3]
end

function phi1_prime(index::Float64,radius::Float64,neff::Float64,r::Union{Float64,Vector{Float64}},k0::Float64,nu::Int64)
    U=sqrt(abs((neff^2)-(index^2)))*k0*radius;
    if (neff<index)
        return besselj_prime.(nu,U.*r./radius)./besselj(nu,U).*U./radius;
    end
    if (neff==index)
        return (r./radius).^(nu-1).*nu./radius;
    end
    return besseli_prime.(nu,U.*r./radius)./besseli(nu,U).*U./radius;
end

function phi2(index::Float64,radius::Float64,neff::Vector{Float64},r::Float64,k0::Float64,nu::Int64)
    U=sqrt.(abs.((neff.^2).-(index^2)))*k0*radius;
    U1=U[neff.<index];
    U2=U[neff.==index];
    U3=U[neff.>index];
    y1=bessely.(nu,U1*r/radius)./bessely.(nu,U1);
    y2=(radius/r)^nu*ones(size(U2));
    y3=besselk.(nu,U3*r/radius)./besselk.(nu,U3);
    return [y1;y2;y3]
end

function phi2(index::Float64,radius::Float64,neff::Float64,r::Union{Float64,Vector{Float64}},k0::Float64,nu::Int64)
    U=sqrt(abs((neff^2)-(index^2)))*k0*radius;
    if (neff<index)
        return bessely.(nu,U.*r./radius)./bessely(nu,U);
    end
    if (neff==index)
        return (radius./r).^nu;
    end
    return besselk.(nu,U.*r./radius)./besselk(nu,U);
end

function phi2_prime(index::Float64,radius::Float64,neff::Vector{Float64},r::Float64,k0::Float64,nu::Int64)
    U=sqrt.(abs.((neff.^2).-(index^2)))*k0*radius;
    U1=U[neff.<index];
    U2=U[neff.==index];
    U3=U[neff.>index];
    y1=bessely_prime.(nu,U1*r/radius)./bessely.(nu,U1).*U1/radius;
    y2=-(r/radius)^(-nu-1)*nu/radius*ones(size(U2));
    y3=besselk_prime.(nu,U3*r/radius)./besselk.(nu,U3).*U3/radius;
    return [y1;y2;y3]
end

function phi2_prime(index::Float64,radius::Float64,neff::Float64,r::Union{Float64,Vector{Float64}},k0::Float64,nu::Int64)
    U=sqrt(abs((neff^2)-(index^2)))*k0*radius;
    if (neff<index)
        return bessely_prime.(nu,U.*r./radius)./bessely(nu,U).*U./radius;
    end
    if (neff==index)
        return -(r./radius).^(-nu-1).*nu./radius;
    end
    return besselk_prime.(nu,U.*r./radius)./besselk(nu,U).*U./radius;
end

function EHz1(index::Float64,radius::Float64,neff::Vector{Float64},r::Float64,k0::Float64,nu::Int64)
#idem phi1
    U=sqrt.(abs.((neff.^2).-(index^2)))*k0*radius;
    U1=U[neff.<index];
    U2=U[neff.==index];
    U3=U[neff.>index];
    y1=besselj.(nu,U1*r/radius)./besselj.(nu,U1);
    y2=(r/radius)^nu*ones(size(U2));
    y3=besseli.(nu,U3*r/radius)./besseli.(nu,U3);
    return [y1;y2;y3]
end

function EHz1(index::Float64,radius::Float64,neff::Float64,r::Union{Float64,Vector{Float64}},k0::Float64,nu::Int64)
    U=sqrt(abs((neff^2)-(index^2)))*k0*radius;
    if (neff<index)
        return besselj.(nu,U.*r./radius)./besselj(nu,U);
    end
    if (neff==index)
        return (r./radius).^nu;
    end
    return besseli.(nu,U.*r./radius)./besseli(nu,U);
end

function EHz2(index::Float64,radius::Float64,neff::Vector{Float64},r::Float64,k0::Float64,nu::Int64)
    #idem phi2
    U=sqrt.(abs.((neff.^2).-(index^2)))*k0*radius;
    U1=U[neff.<index];
    U2=U[neff.==index];
    U3=U[neff.>index];
    y1=bessely.(nu,U1*r/radius)./bessely.(nu,U1);
    y2=(radius/r)^nu*ones(size(U2));
    y3=besselk.(nu,U3*r/radius)./besselk.(nu,U3);
    return [y1;y2;y3]
end

function EHz2(index::Float64,radius::Float64,neff::Float64,r::Union{Float64,Vector{Float64}},k0::Float64,nu::Int64)
    U=sqrt(abs((neff^2)-(index^2)))*k0*radius;
    if (neff<index)
        return bessely.(nu,U.*r./radius)./bessely(nu,U);
    end
    if (neff==index)
        return (radius./r).^nu;
    end
    return besselk.(nu,U.*r./radius)./besselk(nu,U);
end

function Ephi1(index::Float64,radius::Float64,neff::Vector{Float64},r::Float64,k0::Float64,nu::Int64)
    facteur=1.0./((k0^2*index^2).-(k0^2*neff.^2)).*neff*k0;
    U=sqrt.(abs.((neff.^2).-(index^2)))*k0*radius;
    U1=U[neff.<index];
    U2=U[neff.==index];
    U3=U[neff.>index];
    f1=facteur[neff.<index];
    f2=facteur[neff.==index];
    f3=facteur[neff.>index];
    y1=f1.*besselj.(nu,U1*r/radius)./besselj.(nu,U1)*nu/r;
    if (nu==0)
        y2=zeros(size(U2));
    else
        y2=f2.*(r/radius)^nu*nu/r;
    end
    y3=f3.*besseli.(nu,U3*r/radius)./besseli.(nu,U3)*nu/r;
    return [y1;y2;y3]
end

function Ephi1(index::Float64,radius::Float64,neff::Float64,r::Vector{Float64},k0::Float64,nu::Int64)
    facteur=1.0/(k0^2*index^2-k0^2*neff^2)*neff*k0;
    U=sqrt(abs((neff^2)-(index^2)))*k0*radius;
    r1=r[r.!=0];
    r0=r[r.==0];
    result=zeros(size(r));
    if (neff<index)
        result1=facteur*besselj.(nu,U*r1/radius)/besselj(nu,U)*nu./r1;
        if (nu==1)
            result0=facteur*U/2.0/radius/besselj(nu,U)*ones(size(r0));
        else
            result0=zeros(size(r0));
        end
    elseif (neff==index)
        result1=facteur*(r1/radius).^nu*nu./r1;
        if (nu==1)
            result0=facteur/radius*ones(size(r0));
        else
            result0=zeros(size(r0));
        end
    else
        result1=facteur*besseli.(nu,U*r1/radius)/besseli(nu,U)*nu./r1;
        if (nu==1)
            result0=facteur*U/2.0/radius/besseli(nu,U)*ones(size(r0));
        else
            result0=zeros(size(r0));
        end
    end
    result[r.!=0]=result1;
    result[r.==0]=result0;
    return result;
end

function Ephi1(index::Float64,radius::Float64,neff::Float64,r::Float64,k0::Float64,nu::Int64)
    facteur=1.0/(k0^2*index^2-k0^2*neff^2)*neff*k0;
    U=sqrt(abs((neff^2)-(index^2)))*k0*radius;
    if (r==0)
        if (neff<index)
            if (nu==1)
                result=facteur*U/2.0/radius/besselj(nu,U);
            else
                result=0.0;
            end
        elseif (neff==index)
            if (nu==1)
                result=facteur/radius;
            else
                result=0.0;;
            end
        else
            if (nu==1)
                result=facteur*U/2.0/radius/besseli(nu,U);
            else
                result=0.0;
            end
        end
    else
        if (neff<index)
            result=facteur*besselj.(nu,U*r/radius)/besselj(nu,U)*nu/r;
        elseif (neff==index)
            result=(r/radius)^nu;
        else
            result=facteur*besseli.(nu,U*r/radius)/besseli(nu,U)*nu/r;
        end
    end
    return result;
end

function Ephi2(index::Float64,radius::Float64,neff::Vector{Float64},r::Float64,k0::Float64,nu::Int64)
    facteur=1.0./((k0^2*index^2).-(k0^2*neff.^2)).*neff*k0/r;
    U=sqrt.(abs.((neff.^2).-(index^2)))*k0*radius;
    U1=U[neff.<index];
    U2=U[neff.==index];
    U3=U[neff.>index];
    f1=facteur[neff.<index];
    f2=facteur[neff.==index];
    f3=facteur[neff.>index];
    y1=f1.*bessely.(nu,U1*r/radius)./bessely.(nu,U1)*nu;
    if (nu==0)
        y2=zeros(size(U2));
    else
        y2=f2*(r/radius)^(-nu)*nu;
    end
    y3=f3.*besselk.(nu,U3*r/radius)./besselk.(nu,U3)*nu;
    return [y1;y2;y3]
end

function Ephi2(index::Float64,radius::Float64,neff::Float64,r::Union{Float64,Vector{Float64}},k0::Float64,nu::Int64)
    facteur=1.0/(k0^2*index^2-k0^2*neff^2)*neff*k0;
    U=sqrt(abs((neff^2)-(index^2)))*k0*radius;
    if (neff<index)
        result=facteur*bessely.(nu,U*r/radius)/bessely(nu,U)*nu./r;
    elseif (neff==index)
        if (nu==0)
            result=zeros(size(r));
        else
            result=facteur*(r/radius).^(-nu)*nu./r;
        end
    else
        result=facteur*besselk.(nu,U*r/radius)/besselk(nu,U)*nu./r;
    end
    return result;
end

function Ephi3(index::Float64,radius::Float64,neff::Vector{Float64},r::Float64,k0::Float64,nu::Int64)
    facteur=sqrt(mu0/eps0)./((k0^2*index^2).-(k0^2*neff.^2))*k0;
    U=sqrt.(abs.((neff.^2).-(index^2)))*k0*radius;
    U1=U[neff.<index];
    U2=U[neff.==index];
    U3=U[neff.>index];
    f1=facteur[neff.<index];
    f2=facteur[neff.==index];
    f3=facteur[neff.>index];
    y1=-f1.*besselj_prime.(nu,U1*r/radius)./besselj.(nu,U1).*U1/radius;
    if (nu==0)
        y2=sqrt(mu0/eps0)*k0*r/2.0*ones(size(U2));
    else
        y2=-f2*(r/radius)^(nu-1)*nu/radius;
    end
    y3=-f3.*besseli_prime.(nu,U3*r/radius)./besseli.(nu,U3).*U3/radius;
    return [y1;y2;y3]
end

function Ephi3(index::Float64,radius::Float64,neff::Float64,r::Union{Float64,Vector{Float64}},k0::Float64,nu::Int64)
    facteur=sqrt(mu0/eps0)/(k0^2*index^2-k0^2*neff^2)*k0;
    U=sqrt(abs((neff^2)-(index^2)))*k0*radius;
    if (neff<index)
        result=-facteur*besselj_prime.(nu,U*r/radius)/besselj(nu,U)*U/radius;
    elseif (neff==index)
        if (nu==0)
            result=sqrt(mu0/eps0)*k0*r/2.0;
        else
            result=-facteur*(r/radius).^(nu-1)*nu/radius;
        end
    else
        result=-facteur*besseli_prime.(nu,U*r/radius)/besseli(nu,U)*U/radius;
    end
    return result;
end

function Ephi4(index::Float64,radius::Float64,neff::Vector{Float64},r::Float64,k0::Float64,nu::Int64)
    facteur=sqrt(mu0/eps0)./((k0^2*index^2).-(k0^2*neff.^2))*k0;
    U=sqrt.(abs.((neff.^2).-(index^2)))*k0*radius;
    U1=U[neff.<index];
    U2=U[neff.==index];
    U3=U[neff.>index];
    f1=facteur[neff.<index];
    f2=facteur[neff.==index];
    f3=facteur[neff.>index];
    y1=-f1.*bessely_prime.(nu,U1*r/radius)./bessely.(nu,U1).*U1/radius;
    if (nu==0)
        y2=-f2/r;
    else
        y2=f2*(r/radius).^(-nu-1)*nu/radius;
    end
    y3=-f3.*besselk_prime.(nu,U3*r/radius)./besselk.(nu,U3).*U3/radius;
    return [y1;y2;y3]
end

function Ephi4(index::Float64,radius::Float64,neff::Float64,r::Union{Float64,Vector{Float64}},k0::Float64,nu::Int64)
    facteur=sqrt(mu0/eps0)/(k0^2*index^2-k0^2*neff^2)*k0;
    U=sqrt(abs((neff^2)-(index^2)))*k0*radius;
    if (neff<index)
        result=-facteur*bessely_prime.(nu,U*r/radius)/bessely(nu,U)*U/radius;
    elseif (neff==index)
        if (nu==0)
            result=-facteur./r;
        else
            result=-facteur*(r/radius).^(-nu-1)*nu/radius;
        end
    else
        result=-facteur*besselk_prime.(nu,U*r/radius)/besselk(nu,U)*U/radius;
    end
    return result;
end

function Hphi1(index::Float64,radius::Float64,neff::Vector{Float64},r::Float64,k0::Float64,nu::Int64)
    facteur=sqrt(eps0/mu0)./((k0^2*index^2).-(k0^2*neff.^2))*k0*index^2;
    U=sqrt.(abs.((neff.^2).-(index^2)))*k0*radius;
    U1=U[neff.<index];
    U2=U[neff.==index];
    U3=U[neff.>index];
    f1=facteur[neff.<index];
    f2=facteur[neff.==index];
    f3=facteur[neff.>index];
    y1=f1.*besselj_prime.(nu,U1*r/radius)./besselj.(nu,U1).*U1/radius;
    if (nu==0)
        y2=-sqrt(eps0/mu0)*k0*r/2.0*index^2*ones(size(U2));
    else
        y2=f2.*(r/radius)^(nu-1)*nu/radius;
    end
    y3=f3.*besseli_prime.(nu,U3*r/radius)./besseli.(nu,U3).*U3/radius;
    return [y1;y2;y3]
end

function Hphi1(index::Float64,radius::Float64,neff::Float64,r::Union{Float64,Vector{Float64}},k0::Float64,nu::Int64)
    facteur=sqrt(eps0/mu0)/(k0^2*index^2-k0^2*neff^2)*k0*index^2;
    U=sqrt(abs((neff^2)-(index^2)))*k0*radius;
    if (neff<index)
        result=facteur*besselj_prime.(nu,U*r/radius)/besselj(nu,U)*U/radius;
    elseif (neff==index)
        if (nu==0)
            result=-sqrt(eps0/mu0)*k0*r/2.0*index^2;
        else
            result=facteur*(r/radius).^(nu-1)*nu/radius;
        end
    else
        result=facteur*besseli_prime.(nu,U*r/radius)/besseli(nu,U)*U/radius;
    end
    return result;
end

function Hphi2(index::Float64,radius::Float64,neff::Vector{Float64},r::Float64,k0::Float64,nu::Int64)
    facteur=sqrt(eps0/mu0)./((k0^2*index^2).-(k0^2*neff.^2))*k0*index^2;
    U=sqrt.(abs.((neff.^2).-(index^2)))*k0*radius;
    U1=U[neff.<index];
    U2=U[neff.==index];
    U3=U[neff.>index];
    f1=facteur[neff.<index];
    f2=facteur[neff.==index];
    f3=facteur[neff.>index];
    y1=f1.*bessely_prime.(nu,U1*r/radius)./bessely.(nu,U1).*U1/radius;
    if (nu==0)
        y2=f2/r;
    else
        y2=-f2.*(r/radius)^(-nu-1)*nu/radius;
    end
    y3=f3.*besselk_prime.(nu,U3*r/radius)./besselk.(nu,U3).*U3/radius;
    return [y1;y2;y3]
end

function Hphi2(index::Float64,radius::Float64,neff::Float64,r::Union{Float64,Vector{Float64}},k0::Float64,nu::Int64)
    facteur=sqrt(eps0/mu0)/(k0^2*index^2-k0^2*neff^2)*k0*index^2;
    U=sqrt(abs((neff^2)-(index^2)))*k0*radius;
    if (neff<index)
        result=facteur*bessely_prime.(nu,U*r/radius)/bessely(nu,U)*U/radius;
    elseif (neff==index)
        if (nu==0)
            result=facteur./r;
        else
            result=-facteur*(r/radius).^(-nu-1)*nu/radius;
        end
    else
        result=facteur*besselk_prime.(nu,U*r/radius)/besselk(nu,U)*U/radius;
    end
    return result;
end

function Hphi3(index::Float64,radius::Float64,neff::Vector{Float64},r::Float64,k0::Float64,nu::Int64)
    facteur=1.0./((k0^2*index^2).-(k0^2*neff.^2))*k0.*neff;
    U=sqrt.(abs.((neff.^2).-(index^2)))*k0*radius;
    U1=U[neff.<index];
    U2=U[neff.==index];
    U3=U[neff.>index];
    f1=facteur[neff.<index];
    f2=facteur[neff.==index];
    f3=facteur[neff.>index];
    y1=-f1.*besselj.(nu,U1*r/radius)./besselj.(nu,U1)*nu/r;
    if (nu==0)
        y2=zeros(size(U2));
    else
        y2=f2.*(r/radius)^(nu)*nu/r;
    end
    y3=-f3.*besseli.(nu,U3*r/radius)./besseli.(nu,U3)*nu/r;
    return [y1;y2;y3]
end

function Hphi3(index::Float64,radius::Float64,neff::Float64,r::Vector{Float64},k0::Float64,nu::Int64)
    facteur=1.0/(k0^2*index^2-k0^2*neff^2)*k0*neff;
    U=sqrt(abs((neff^2)-(index^2)))*k0*radius;
    r1=r[r.!=0];
    r0=r[r.==0];
    result=zeros(size(r));
    if (neff<index)
        result1=-facteur*besselj.(nu,U*r1/radius)/besselj(nu,U)*nu./r1;
        if (nu==1)
            result0=-facteur*U/2.0/radius/besselj(nu,U)*ones(size(r0));
        else
            result0=zeros(size(r0));
        end
    elseif (neff==index)
        if (nu==0)
            result1=zeros(size(r1));
        else
            result1=facteur*(r1/radius).^(nu)*nu./r1;
        end
        if (nu==1)
            result0=facteur/radius*ones(size(r0));
        else
            result0=zeros(size(r0));
        end
    else
        result1=-facteur*besseli.(nu,U*r1/radius)/besseli(nu,U)*nu./r1;
        if (nu==1)
            result0=-facteur*U/2.0/radius/besseli(nu,U)*ones(size(r0));
        else
            result0=zeros(size(r0));
        end
    end
    result[r.!=0]=result1;
    result[r.==0]=result0;
    return result;
end

function Hphi3(index::Float64,radius::Float64,neff::Float64,r::Float64,k0::Float64,nu::Int64)
    facteur=1.0/(k0^2*index^2-k0^2*neff^2)*k0*neff;
    U=sqrt(abs((neff^2)-(index^2)))*k0*radius;
    if r==0
        if (neff<index)
            if (nu==1)
                result=-facteur*U/2.0/radius/besselj(nu,U);
            else
                result=0.0;
            end
        elseif (neff==index)
            if (nu==1)
                result=facteur/radius;
            else
                result=0.0;
            end
        else
            if (nu==1)
                result=-facteur*U/2.0/radius/besseli(nu,U);
            else
                result=0.0;
            end
        end
    else
        if (neff<index)
            result=-facteur*besselj.(nu,U*r/radius)/besselj(nu,U)*nu/r;
        elseif (neff==index)
            if (nu==0)
                result=0.0;
            else
                result=facteur*(r/radius)^(nu)*nu/r;
            end
        else
            result=facteur*besseli.(nu,U*r/radius)/besseli(nu,U)*nu/r;
        end
    end
    return result;
end

function Hphi4(index::Float64,radius::Float64,neff::Vector{Float64},r::Float64,k0::Float64,nu::Int64)
    facteur=1.0./((k0^2*index^2).-(k0^2*neff.^2))*k0.*neff/r;
    U=sqrt.(abs.((neff.^2).-(index^2)))*k0*radius;
    U1=U[neff.<index];
    U2=U[neff.==index];
    U3=U[neff.>index];
    f1=facteur[neff.<index];
    f2=facteur[neff.==index];
    f3=facteur[neff.>index];
    y1=-f1.*bessely.(nu,U1*r/radius)./bessely.(nu,U1)*nu;
    if (nu==0)
        y2=zeros(size(U2));
    else
        y2=f2.*(r/radius)^(-nu)*nu;
    end
    y3=-f3.*besselk.(nu,U3*r/radius)./besselk.(nu,U3)*nu;
    return [y1;y2;y3]
end

function Hphi4(index::Float64,radius::Float64,neff::Float64,r::Union{Float64,Vector{Float64}},k0::Float64,nu::Int64)
    facteur=1.0/(k0^2*index^2-k0^2*neff^2)*k0*neff;
    U=sqrt(abs((neff^2)-(index^2)))*k0*radius;
    if (neff<index)
        result=-facteur*bessely.(nu,U*r/radius)/bessely(nu,U)*nu./r;
    elseif (neff==index)
        if (nu==0)
            result=zeros(size(r));
        else
            result=facteur*(r/radius).^(-nu)*nu./r;
        end
    else
        result=-facteur*besselk.(nu,U*r/radius)/besselk(nu,U)*nu./r;
    end
    return result;
end

function Er1(index::Float64,radius::Float64,neff::Float64,r::Union{Float64,Vector{Float64}},k0::Float64,nu::Int64)
    facteur=1.0/(k0^2*index^2-k0^2*neff^2)*k0*neff;
    U=sqrt(abs((neff^2)-(index^2)))*k0*radius;
    if (neff<index)
        result=facteur*besselj_prime.(nu,U*r/radius)/besselj(nu,U)*U/radius;
    elseif (neff==index)
        if (nu==0)
            result=-k0*r/2.0*neff;
        else
            result=facteur*(r/radius).^(nu-1)*nu/radius;
        end
    else
        result=facteur*besseli_prime.(nu,U*r/radius)/besseli(nu,U)*U/radius;
    end
    return result;
end

function Er2(index::Float64,radius::Float64,neff::Float64,r::Union{Float64,Vector{Float64}},k0::Float64,nu::Int64)
    facteur=1.0/(k0^2*index^2-k0^2*neff^2)*k0*neff;
    U=sqrt(abs((neff^2)-(index^2)))*k0*radius;
    if (neff<index)
        result=facteur*bessely_prime.(nu,U*r/radius)/bessely(nu,U)*U/radius;
    elseif (neff==index)
        if (nu==0)
            result=facteur./r;
        else
            result=-facteur*(r/radius).^(-nu-1)*nu/radius;
        end
    else
        result=facteur*besselk_prime.(nu,U*r/radius)/besselk(nu,U)*U/radius;
    end
    return result;
end

function Er3(index::Float64,radius::Float64,neff::Float64,r::Union{Float64,Vector{Float64}},k0::Float64,nu::Int64)
    facteur=sqrt(mu0/eps0)/(k0^2*index^2-k0^2*neff^2)*k0;
    U=sqrt(abs((neff^2)-(index^2)))*k0*radius;
    r1=r[r.!=0];
    r0=r[r.==0];
    result=zeros(size(r));
    result1=-nu*facteur./r1.*EHz1(index,radius,neff,r1,k0,nu);
    if (nu==1)
        if neff<index
            result0=-facteur*U/2.0/radius/besselj(1,U)*ones(size(r0));
        elseif neff==index
            result0=-facteur/radius*ones(size(r0))
        else
            result0=-facteur*U/2.0/radius/besseli(1,U)*ones(size(r0));
        end
    else
        result0=zeros(size(r0));
    end
    result[r.!=0]=result1;
    result[r.==0]=result0;
    return result;
end

function Er4(index::Float64,radius::Float64,neff::Float64,r::Union{Float64,Vector{Float64}},k0::Float64,nu::Int64)
    facteur=sqrt(mu0/eps0)/(k0^2*index^2-k0^2*neff^2)*k0;
    return -facteur*nu./r.*EHz2(index,radius,neff,r,k0,nu);
end

function Hr1(index::Float64,radius::Float64,neff::Float64,r::Union{Float64,Vector{Float64}},k0::Float64,nu::Int64)
    facteur=sqrt(eps0/mu0)/(k0^2*index^2-k0^2*neff^2)*k0*index^2;
    U=sqrt(abs((neff^2)-(index^2)))*k0*radius;
    r1=r[r.!=0];
    r0=r[r.==0];
    result=zeros(size(r));
    result1=-nu*facteur./r1.*EHz1(index,radius,neff,r1,k0,nu);
    if (nu==1)
        if neff<index
            result0=-facteur*U/2.0/radius/besselj(1,U)*ones(size(r0));
        elseif neff==index
            result0=-facteur/radius*ones(size(r0))
        else
            result0=-facteur*U/2.0/radius/besseli(1,U)*ones(size(r0));
        end
    else
        result0=zeros(size(r0));
    end
    result[r.!=0]=result1;
    result[r.==0]=result0;
    return result;
end

function Hr2(index::Float64,radius::Float64,neff::Float64,r::Union{Float64,Vector{Float64}},k0::Float64,nu::Int64)
    facteur=sqrt(eps0/mu0)/(k0^2*index^2-k0^2*neff^2)*k0*index^2;
    return -facteur*nu./r.*EHz2(index,radius,neff,r,k0,nu);
end

function Hr3(index::Float64,radius::Float64,neff::Float64,r::Union{Float64,Vector{Float64}},k0::Float64,nu::Int64)
    facteur=1.0/(k0^2*index^2-k0^2*neff^2)*k0*neff;
    U=sqrt(abs((neff^2)-(index^2)))*k0*radius;
    if (neff<index)
        result=facteur*besselj_prime.(nu,U*r/radius)/besselj(nu,U)*U/radius;
    elseif (neff==index)
        if (nu==0)
            result=-k0*r/2.0*neff;
        else
            result=facteur*(r/radius).^(nu-1)*nu/radius;
        end
    else
        result=facteur*besseli_prime.(nu,U*r/radius)/besseli(nu,U)*U/radius;
    end
    return result;
end

function Hr4(index::Float64,radius::Float64,neff::Float64,r::Union{Float64,Vector{Float64}},k0::Float64,nu::Int64)
    facteur=1.0/(k0^2*index^2-k0^2*neff^2)*k0*neff;
    U=sqrt(abs((neff^2)-(index^2)))*k0*radius;
    if (neff<index)
        result=facteur*bessely_prime.(nu,U*r/radius)/bessely(nu,U)*U/radius;
    elseif (neff==index)
        if (nu==0)
            result=facteur./r;
        else
            result=-facteur*(r/radius).^(-nu-1)*nu/radius;
        end
    else
        result=facteur*besselk_prime.(nu,U*r/radius)/besselk(nu,U)*U/radius;
    end
    return result;
end

function scalarField_test(r::Float64,nu::Int64,neff::Float64,k0::Float64,radius::Union{Vector{Float64},Float64},index::Vector{Float64},A::Vector{Float64},B::Vector{Float64})
    nb_couches=length(index);
    pos=1;
    @inbounds for i=1:length(radius)
        pos=pos+(r>radius[i]);
    end
    if (pos==1)
        E=A[1]*phi1(index[1],radius[1],neff,r,k0,nu);
    elseif (pos==nb_couches)
        E=B[nb_couches]*phi2(index[nb_couches],radius[nb_couches-1],neff,r,k0,nu);
    else
        E=A[pos].*phi1(index[pos],radius[pos],neff,r,k0,nu)+B[pos].*phi2(index[pos],radius[pos],neff,r,k0,nu);
    end
    return E;
end

function scalarField(r::Union{Array{Float64,2},Vector{Float64},Float64},nu::Int64,neff::Float64,k0::Float64,radius::Union{Vector{Float64},Float64},index::Vector{Float64},A::Vector{Float64},B::Vector{Float64})
    nb_couches=length(index);
    pos=fill(1,size(r));
    @inbounds for i=1:length(radius)
        @. pos=pos+(r>radius[i]);
    end
    E=Array{Float64}(undef,size(r));
    E[pos.==1].=A[1].*phi1(index[1],radius[1],neff,r[pos.==1],k0,nu);
    E[pos.==nb_couches].=B[nb_couches].*phi2(index[nb_couches],radius[nb_couches-1],neff,r[pos.==nb_couches],k0,nu);
    @inbounds for i=2:nb_couches-1
        #E[pos.==i].=A[i].*phi1(index[i],radius[i],neff,r[pos.==i],k0,nu).+B[i].*phi2(index[i],radius[i],neff,r[pos.==i],k0,nu);
        Ei=@view E[pos.==i];
        ri=@view r[pos.==i];
        @. Ei=A[i]*phi1(index[i],radius[i],neff,ri,k0,nu)+B[i]*phi2(index[i],radius[i],neff,ri,k0,nu);
    end
    return E;
end

function cossin(phi::Float64,nu::Int64,sc::Bool)
    if (sc)
        return sin(nu*phi);
    else
        return cos(nu*phi);
    end
end

function cossin2(phi::Float64,nu::Int64,sc::Bool)
    if (sc)
        return cos(nu*phi);
    else
        return -sin(nu*phi);
    end
end

function Er(r::Union{Array{Float64,2},Vector{Float64},Float64},phi::Union{Array{Float64,2},Vector{Float64},Float64},nu::Int64,neff::Float64,k0::Float64,radius::Union{Vector{Float64},Float64},index::Vector{Float64},A1::Vector{Float64},A2::Vector{Float64},A3::Vector{Float64},A4::Vector{Float64},sc::Bool)
    nb_couches=length(index);
    pos=fill(1,size(r));
    @inbounds for i=1:length(radius)
        pos=pos.+(r.>radius[i]);
    end
    Er=Array{Float64}(undef,size(r));
    Er1_task=Threads.@spawn Er1(index[1],radius[1],neff,r[pos.==1],k0,nu);
    Er2_task=Threads.@spawn Er2(index[end],radius[end],neff,r[pos.==nb_couches],k0,nu);
    Er3_task=Threads.@spawn Er3(index[1],radius[1],neff,r[pos.==1],k0,nu);
    Er4_task=Threads.@spawn Er4(index[end],radius[end],neff,r[pos.==nb_couches],k0,nu);
    #Er[pos.==1].=A1[1]*Er1(index[1],radius[1],neff,r[pos.==1],k0,nu)+A3[1]*Er3(index[1],radius[1],neff,r[pos.==1],k0,nu);
    Er[pos.==1].=A1[1]*fetch(Er1_task)+A3[1]*fetch(Er3_task);
    #Er[pos.==nb_couches].=A2[end]*Er2(index[end],radius[end],neff,r[pos.==nb_couches],k0,nu)+A4[end]*Er4(index[end],radius[end],neff,r[pos.==nb_couches],k0,nu);
    Er[pos.==nb_couches].=A2[end]*fetch(Er2_task)+A4[end]*fetch(Er4_task);
    @inbounds for i=2:nb_couches-1
        Er1_task=Threads.@spawn Er1(index[i],radius[i],neff,r[pos.==i],k0,nu);
        Er2_task=Threads.@spawn Er2(index[i],radius[i],neff,r[pos.==i],k0,nu);
        Er3_task=Threads.@spawn Er3(index[i],radius[i],neff,r[pos.==i],k0,nu);
        Er4_task=Threads.@spawn Er4(index[i],radius[i],neff,r[pos.==i],k0,nu);
        Er[pos.==i].=A1[i]*fetch(Er1_task)+A2[i]*fetch(Er2_task)+A3[i]*fetch(Er3_task)+A4[i]*fetch(Er4_task);
        #Er[pos.==i].=A1[i]*Er1(index[i],radius[i],neff,r[pos.==i],k0,nu)+A2[i]*Er2(index[i],radius[i],neff,r[pos.==i],k0,nu)+A3[i]*Er3(index[i],radius[i],neff,r[pos.==i],k0,nu)+A4[i]*Er4(index[i],radius[i],neff,r[pos.==i],k0,nu);
    end
    return cossin.(phi,nu,sc).*Er;
end

function Hr(r::Union{Array{Float64,2},Vector{Float64},Float64},phi::Union{Array{Float64,2},Vector{Float64},Float64},nu::Int64,neff::Float64,k0::Float64,radius::Union{Vector{Float64},Float64},index::Vector{Float64},A1::Vector{Float64},A2::Vector{Float64},A3::Vector{Float64},A4::Vector{Float64},sc::Bool)
    nb_couches=length(index);
    pos=fill(1,size(r));
    @inbounds for i=1:length(radius)
        pos=pos.+(r.>radius[i]);
    end
    Hr=Array{Float64}(undef,size(r));
    Hr1_task=Threads.@spawn Hr1(index[1],radius[1],neff,r[pos.==1],k0,nu);
    Hr2_task=Threads.@spawn Hr2(index[end],radius[end],neff,r[pos.==nb_couches],k0,nu);
    Hr3_task=Threads.@spawn Hr3(index[1],radius[1],neff,r[pos.==1],k0,nu);
    Hr4_task=Threads.@spawn Hr4(index[end],radius[end],neff,r[pos.==nb_couches],k0,nu);
    #Er[pos.==1].=A1[1]*Er1(index[1],radius[1],neff,r[pos.==1],k0,nu)+A3[1]*Er3(index[1],radius[1],neff,r[pos.==1],k0,nu);
    Hr[pos.==1].=A1[1]*fetch(Hr1_task)+A3[1]*fetch(Hr3_task);
    #Er[pos.==nb_couches].=A2[end]*Er2(index[end],radius[end],neff,r[pos.==nb_couches],k0,nu)+A4[end]*Er4(index[end],radius[end],neff,r[pos.==nb_couches],k0,nu);
    Hr[pos.==nb_couches].=A2[end]*fetch(Hr2_task)+A4[end]*fetch(Hr4_task);
    @inbounds for i=2:nb_couches-1
        Hr1_task=Threads.@spawn Hr1(index[i],radius[i],neff,r[pos.==i],k0,nu);
        Hr2_task=Threads.@spawn Hr2(index[i],radius[i],neff,r[pos.==i],k0,nu);
        Hr3_task=Threads.@spawn Hr3(index[i],radius[i],neff,r[pos.==i],k0,nu);
        Hr4_task=Threads.@spawn Hr4(index[i],radius[i],neff,r[pos.==i],k0,nu);
        Hr[pos.==i].=A1[i]*fetch(Hr1_task)+A2[i]*fetch(Hr2_task)+A3[i]*fetch(Hr3_task)+A4[i]*fetch(Hr4_task);
        #Er[pos.==i].=A1[i]*Er1(index[i],radius[i],neff,r[pos.==i],k0,nu)+A2[i]*Er2(index[i],radius[i],neff,r[pos.==i],k0,nu)+A3[i]*Er3(index[i],radius[i],neff,r[pos.==i],k0,nu)+A4[i]*Er4(index[i],radius[i],neff,r[pos.==i],k0,nu);
    end
    return cossin2.(phi,nu,sc).*Hr;
end

function Ephi(r::Union{Array{Float64,2},Vector{Float64},Float64},phi::Union{Array{Float64,2},Vector{Float64},Float64},nu::Int64,neff::Float64,k0::Float64,radius::Union{Vector{Float64},Float64},index::Vector{Float64},A1::Vector{Float64},A2::Vector{Float64},A3::Vector{Float64},A4::Vector{Float64},sc::Bool)
    nb_couches=length(index);
    pos=fill(1,size(r));
    @inbounds for i=1:length(radius)
        pos=pos.+(r.>radius[i]);
    end
    Ephi=Array{Float64}(undef,size(r));
    Ephi1_task=Threads.@spawn Ephi1(index[1],radius[1],neff,r[pos.==1],k0,nu);
    Ephi2_task=Threads.@spawn Ephi2(index[end],radius[end],neff,r[pos.==nb_couches],k0,nu);
    Ephi3_task=Threads.@spawn Ephi3(index[1],radius[1],neff,r[pos.==1],k0,nu);
    Ephi4_task=Threads.@spawn Ephi4(index[end],radius[end],neff,r[pos.==nb_couches],k0,nu);
    #Ephi[pos.==1].=A1[1]*Ephi1(index[1],radius[1],neff,r[pos.==1],k0,nu)+A3[1]*Ephi3(index[1],radius[1],neff,r[pos.==1],k0,nu);
    Ephi[pos.==1].=A1[1]*fetch(Ephi1_task)+A3[1]*fetch(Ephi3_task);
    #Ephi[pos.==nb_couches].=A2[end]*Ephi2(index[end],radius[end],neff,r[pos.==nb_couches],k0,nu)+A4[end]*Ephi4(index[end],radius[end],neff,r[pos.==nb_couches],k0,nu);
    Ephi[pos.==nb_couches].=A2[end]*fetch(Ephi2_task)+A4[end]*fetch(Ephi4_task);
    @inbounds for i=2:nb_couches-1
        Ephi1_task=Threads.@spawn Ephi1(index[i],radius[i],neff,r[pos.==i],k0,nu);
        Ephi2_task=Threads.@spawn Ephi2(index[i],radius[i],neff,r[pos.==i],k0,nu);
        Ephi3_task=Threads.@spawn Ephi3(index[i],radius[i],neff,r[pos.==i],k0,nu);
        Ephi4_task=Threads.@spawn Ephi4(index[i],radius[i],neff,r[pos.==i],k0,nu);
        Ephi[pos.==i].=A1[i]*fetch(Ephi1_task)+A2[i]*fetch(Ephi2_task)+A3[i]*fetch(Ephi3_task)+A4[i]*fetch(Ephi4_task);
        #Ephi[pos.==i].=A1[i]*Ephi1(index[i],radius[i],neff,r[pos.==i],k0,nu)+A2[i]*Ephi2(index[i],radius[i],neff,r[pos.==i],k0,nu)+A3[i]*Ephi3(index[i],radius[i],neff,r[pos.==i],k0,nu)+A4[i]*Ephi4(index[i],radius[i],neff,r[pos.==i],k0,nu);
    end
    return cossin2.(phi,nu,sc).*Ephi;
end

function Hphi(r::Union{Array{Float64,2},Vector{Float64},Float64},phi::Union{Array{Float64,2},Vector{Float64},Float64},nu::Int64,neff::Float64,k0::Float64,radius::Union{Vector{Float64},Float64},index::Vector{Float64},A1::Vector{Float64},A2::Vector{Float64},A3::Vector{Float64},A4::Vector{Float64},sc::Bool)
    nb_couches=length(index);
    pos=fill(1,size(r));
    @inbounds for i=1:length(radius)
        pos=pos.+(r.>radius[i]);
    end
    Hphi=Array{Float64}(undef,size(r));
    Hphi1_task=Threads.@spawn Hphi1(index[1],radius[1],neff,r[pos.==1],k0,nu);
    Hphi2_task=Threads.@spawn Hphi2(index[end],radius[end],neff,r[pos.==nb_couches],k0,nu);
    Hphi3_task=Threads.@spawn Hphi3(index[1],radius[1],neff,r[pos.==1],k0,nu);
    Hphi4_task=Threads.@spawn Hphi4(index[end],radius[end],neff,r[pos.==nb_couches],k0,nu);
    #Ephi[pos.==1].=A1[1]*Ephi1(index[1],radius[1],neff,r[pos.==1],k0,nu)+A3[1]*Ephi3(index[1],radius[1],neff,r[pos.==1],k0,nu);
    Hphi[pos.==1].=A1[1]*fetch(Hphi1_task)+A3[1]*fetch(Hphi3_task);
    #Ephi[pos.==nb_couches].=A2[end]*Ephi2(index[end],radius[end],neff,r[pos.==nb_couches],k0,nu)+A4[end]*Ephi4(index[end],radius[end],neff,r[pos.==nb_couches],k0,nu);
    Hphi[pos.==nb_couches].=A2[end]*fetch(Hphi2_task)+A4[end]*fetch(Hphi4_task);
    @inbounds for i=2:nb_couches-1
        Hphi1_task=Threads.@spawn Hphi1(index[i],radius[i],neff,r[pos.==i],k0,nu);
        Hphi2_task=Threads.@spawn Hphi2(index[i],radius[i],neff,r[pos.==i],k0,nu);
        Hphi3_task=Threads.@spawn Hphi3(index[i],radius[i],neff,r[pos.==i],k0,nu);
        Hphi4_task=Threads.@spawn Hphi4(index[i],radius[i],neff,r[pos.==i],k0,nu);
        Hphi[pos.==i].=A1[i]*fetch(Hphi1_task)+A2[i]*fetch(Hphi2_task)+A3[i]*fetch(Hphi3_task)+A4[i]*fetch(Hphi4_task);
        #Ephi[pos.==i].=A1[i]*Ephi1(index[i],radius[i],neff,r[pos.==i],k0,nu)+A2[i]*Ephi2(index[i],radius[i],neff,r[pos.==i],k0,nu)+A3[i]*Ephi3(index[i],radius[i],neff,r[pos.==i],k0,nu)+A4[i]*Ephi4(index[i],radius[i],neff,r[pos.==i],k0,nu);
    end
    return cossin.(phi,nu,sc).*Hphi;
end

function EzHz(r::Union{Array{Float64,2},Vector{Float64},Float64},phi::Union{Array{Float64,2},Vector{Float64},Float64},nu::Int64,neff::Float64,k0::Float64,radius::Union{Vector{Float64},Float64},index::Vector{Float64},A1::Vector{Float64},A2::Vector{Float64},A3::Vector{Float64},A4::Vector{Float64},sc::Bool)
    nb_couches=length(index);
    pos=fill(1,size(r));
    @inbounds for i=1:length(radius)
        pos=pos.+(r.>radius[i]);
    end
    Ez=Array{Float64}(undef,size(r));
    Hz=Array{Float64}(undef,size(r));
    EHz1_task=Threads.@spawn EHz1(index[1],radius[1],neff,r[pos.==1],k0,nu);
    EHz2_task=Threads.@spawn EHz2(index[end],radius[end],neff,r[pos.==nb_couches],k0,nu);
    Ez[pos.==1].=A1[1]*fetch(EHz1_task);
    Ez[pos.==nb_couches].=A2[end]*fetch(EHz2_task);
    Hz[pos.==1].=A3[1]*fetch(EHz1_task);
    Hz[pos.==nb_couches].=A4[end]*fetch(EHz2_task);
    @inbounds for i=2:nb_couches-1
        EHz1_task=Threads.@spawn EHz1(index[i],radius[i],neff,r[pos.==i],k0,nu);
        EHz2_task=Threads.@spawn EHz2(index[i],radius[i],neff,r[pos.==i],k0,nu);
        Ez[pos.==i].=A1[i]*fetch(EHz1_task)+A2[i]*fetch(EHz2_task);
        Hz[pos.==i].=A3[i]*fetch(EHz1_task)+A4[i]*fetch(EHz2_task);
    end
    return cossin.(phi,nu,sc).*Ez,cossin2.(phi,nu,sc).*Hz;
end

function ExEy(r::Union{Array{Float64,2},Vector{Float64},Float64},phi::Union{Array{Float64,2},Vector{Float64},Float64},nu::Int64,neff::Float64,k0::Float64,radius::Union{Vector{Float64},Float64},index::Vector{Float64},A1::Vector{Float64},A2::Vector{Float64},A3::Vector{Float64},A4::Vector{Float64},sc::Bool)
    Er_task=Threads.@spawn Er(r,phi,nu,neff,k0,radius,index,A1,A2,A3,A4,sc);
    Ephi_task=Threads.@spawn Ephi(r,phi,nu,neff,k0,radius,index,A1,A2,A3,A4,sc)
    return cos.(phi).*fetch(Er_task)-sin.(phi).*fetch(Ephi_task),sin.(phi).*fetch(Er_task)+cos.(phi).*fetch(Ephi_task);
end

function HxHy(r::Union{Array{Float64,2},Vector{Float64},Float64},phi::Union{Array{Float64,2},Vector{Float64},Float64},nu::Int64,neff::Float64,k0::Float64,radius::Union{Vector{Float64},Float64},index::Vector{Float64},A1::Vector{Float64},A2::Vector{Float64},A3::Vector{Float64},A4::Vector{Float64},sc::Bool)
    Hr_task=Threads.@spawn Hr(r,phi,nu,neff,k0,radius,index,A1,A2,A3,A4,sc);
    Hphi_task=Threads.@spawn Hphi(r,phi,nu,neff,k0,radius,index,A1,A2,A3,A4,sc)
    return cos.(phi).*fetch(Hr_task)-sin.(phi).*fetch(Hphi_task),sin.(phi).*fetch(Hr_task)+cos.(phi).*fetch(Hphi_task);
end

function findAB_scalar(k0::Float64,nu::Int64,neff::Float64,radius::Union{Vector{Float64},Float64},index::Vector{Float64})
    nb_couches=length(index);
    Asol=Vector{Float64}(undef,nb_couches);
    Bsol=Vector{Float64}(undef,nb_couches);
    Asol[1]=1.0;
    Bsol[1]=0.0;
    phi_prec=phi1(index[1],radius[1],neff,radius[1],k0,nu);
    phi_prime_prec=phi1_prime(index[1],radius[1],neff,radius[1],k0,nu);
    @inbounds for couche=2:nb_couches-1
        p1=phi1(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        p1_prime=phi1_prime(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        p2=phi2(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        p2_prime=phi2_prime(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        det=1.0/(p1*p2_prime-p2*p1_prime);
        A=det*(p2_prime*phi_prec-p2*phi_prime_prec);
        B=det*(-p1_prime*phi_prec+p1*phi_prime_prec);
        phi_prec=A+B;
        phi_prime_prec=A*phi1_prime(index[couche],radius[couche],neff,radius[couche],k0,nu)+B*phi2_prime(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Asol[couche]=A;
        Bsol[couche]=B;
    end
    Asol[nb_couches]=0.0;
    Bsol[nb_couches]=phi_prec;
    return Asol,Bsol;
end

function findA_TM(k0::Float64,nu::Int64,neff::Float64,radius::Union{Vector{Float64},Float64},index::Vector{Float64})
    nb_couches=length(index);
    A1=zeros(nb_couches);
    A2=zeros(nb_couches);
    A3=zeros(nb_couches);
    A4=zeros(nb_couches);
    A1[1]=1.0;
    Ez_prec=1.0;
    Hphi_prec=Hphi1(index[1],radius[1],neff,radius[1],k0,nu);
    @inbounds for couche=2:nb_couches-1
        a=EHz1(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        b=EHz2(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        c=Hphi1(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        d=Hphi2(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        det=1.0/(a*d-b*c);
        A=det*(d*Ez_prec-b*Hphi_prec);
        B=det*(a*Hphi_prec-c*Ez_prec);
        Ez_prec=A+B;
        Hphi_prec1=A*Hphi1(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Hphi_prec2=B*Hphi2(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Hphi_prec=Hphi_prec1+Hphi_prec2;
        A1[couche]=A;
        A2[couche]=B;
    end
    A2[nb_couches]=Ez_prec;
    return A1,A2,A3,A4;
end

function findA_TE(k0::Float64,nu::Int64,neff::Float64,radius::Union{Vector{Float64},Float64},index::Vector{Float64})
    nb_couches=length(index);
    A1=zeros(nb_couches);
    A2=zeros(nb_couches);
    A3=zeros(nb_couches);
    A4=zeros(nb_couches);
    A3[1]=1.0;
    Hz_prec=1.0;
    Ephi_prec=Ephi3(index[1],radius[1],neff,radius[1],k0,nu);
    @inbounds for couche=2:nb_couches-1
        a=EHz1(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        b=EHz2(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        c=Ephi3(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        d=Ephi4(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        det=1.0/(a*d-b*c);
        A=det*(d*Hz_prec-b*Ephi_prec);
        B=det*(a*Ephi_prec-c*Hz_prec);
        Hz_prec=A+B;
        Ephi_prec1=A*Ephi3(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Ephi_prec2=B*Ephi4(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Ephi_prec=Ephi_prec1+Ephi_prec2;
        A3[couche]=A;
        A4[couche]=B;
    end
    A4[nb_couches]=Hz_prec;
    return A1,A2,A3,A4;
end

function findA_HE_EH(k0::Float64,nu::Int64,neff::Float64,radius::Union{Vector{Float64},Float64},index::Vector{Float64})
    nb_couches=length(index);
    Ez_prec_1=1.0;
    Hz_prec_1=0.0;
    Ephi_prec_1=Ephi1(index[1],radius[1],neff,radius[1],k0,nu);
    Hphi_prec_1=Hphi1(index[1],radius[1],neff,radius[1],k0,nu);
    Ez_prec_2=1.0;
    Hz_prec_2=1.0;
    Ephi_prec_2=Ephi_prec_1+Ephi3(index[1],radius[1],neff,radius[1],k0,nu);
    Hphi_prec_2=Hphi_prec_1+Hphi3(index[1],radius[1],neff,radius[1],k0,nu);
    @inbounds for couche=2:nb_couches-1
        m11=EHz1(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        m12=EHz2(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        m31=Ephi1(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        m32=Ephi2(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        m33=Ephi3(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        m34=Ephi4(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        m41=Hphi1(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        m42=Hphi2(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        m43=Hphi3(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        m44=Hphi4(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        m23=m11;
        m24=m12;
        det=1.0/(m11*m23*(m32*m44-m34*m42)-m11*m24*(m32*m43-m33*m42)-m12*m23*(m31*m44-m34*m41)+m12*m24*(m31*m43-m33*m41));
        sol1_1=det*((m23*m32*m44-m23*m34*m42-m24*m32*m43+m24*m33*m42)*Ez_prec_1+m12*(m33*m44-m34*m43)*Hz_prec_1-m12*(m23*m44-m24*m43)*Ephi_prec_1+m12*(m23*m34-m24*m33)*Hphi_prec_1);
		sol2_1=det*(-(m23*m31*m44-m23*m34*m41-m24*m31*m43+m24*m33*m41)*Ez_prec_1-m11*(m33*m44-m34*m43)*Hz_prec_1+m11*(m23*m44-m24*m43)*Ephi_prec_1-m11*(m23*m34-m24*m33)*Hphi_prec_1);
		sol3_1=det*(-m24*(m31*m42-m32*m41)*Ez_prec_1+(m11*m32*m44-m11*m34*m42-m12*m31*m44+m12*m34*m41)*Hz_prec_1+m24*(m11*m42-m12*m41)*Ephi_prec_1-m24*(m11*m32-m12*m31)*Hphi_prec_1);
		sol4_1=det*(m23*(m31*m42-m32*m41)*Ez_prec_1-(m11*m32*m43-m11*m33*m42-m12*m31*m43+m12*m33*m41)*Hz_prec_1-m23*(m11*m42-m12*m41)*Ephi_prec_1+m23*(m11*m32-m12*m31)*Hphi_prec_1);
		sol1_2=det*((m23*m32*m44-m23*m34*m42-m24*m32*m43+m24*m33*m42)*Ez_prec_2+m12*(m33*m44-m34*m43)*Hz_prec_2-m12*(m23*m44-m24*m43)*Ephi_prec_2+m12*(m23*m34-m24*m33)*Hphi_prec_2);
		sol2_2=det*(-(m23*m31*m44-m23*m34*m41-m24*m31*m43+m24*m33*m41)*Ez_prec_2-m11*(m33*m44-m34*m43)*Hz_prec_2+m11*(m23*m44-m24*m43)*Ephi_prec_2-m11*(m23*m34-m24*m33)*Hphi_prec_2);
		sol3_2=det*(-m24*(m31*m42-m32*m41)*Ez_prec_2+(m11*m32*m44-m11*m34*m42-m12*m31*m44+m12*m34*m41)*Hz_prec_2+m24*(m11*m42-m12*m41)*Ephi_prec_2-m24*(m11*m32-m12*m31)*Hphi_prec_2);
		sol4_2=det*(m23*(m31*m42-m32*m41)*Ez_prec_2-(m11*m32*m43-m11*m33*m42-m12*m31*m43+m12*m33*m41)*Hz_prec_2-m23*(m11*m42-m12*m41)*Ephi_prec_2+m23*(m11*m32-m12*m31)*Hphi_prec_2);
        Ez_prec_1=sol1_1+sol2_1;
        Hz_prec_1=sol3_1+sol4_1;
        Ephi1v=Ephi1(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Ephi2v=Ephi2(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Ephi3v=Ephi3(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Ephi4v=Ephi4(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Hphi1v=Hphi1(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Hphi2v=Hphi2(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Hphi3v=Hphi3(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Hphi4v=Hphi4(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Ephi_prec_1=sol1_1*Ephi1v+sol2_1*Ephi2v+sol3_1*Ephi3v+sol4_1*Ephi4v;
        Hphi_prec_1=sol1_1*Hphi1v+sol2_1*Hphi2v+sol3_1*Hphi3v+sol4_1*Hphi4v;
        Ez_prec_2=sol1_2+sol2_2;
        Hz_prec_2=sol3_2+sol4_2;
        Ephi_prec_2=sol1_2*Ephi1v+sol2_2*Ephi2v+sol3_2*Ephi3v+sol4_2*Ephi4v;
        Hphi_prec_2=sol1_2*Hphi1v+sol2_2*Hphi2v+sol3_2*Hphi3v+sol4_2*Hphi4v;
    end
    Ephi2v=Ephi2(index[nb_couches],radius[nb_couches-1],neff,radius[nb_couches-1],k0,nu);
    Ephi4v=Ephi4(index[nb_couches],radius[nb_couches-1],neff,radius[nb_couches-1],k0,nu);
    Hphi2v=Hphi2(index[nb_couches],radius[nb_couches-1],neff,radius[nb_couches-1],k0,nu);
    Hphi4v=Hphi4(index[nb_couches],radius[nb_couches-1],neff,radius[nb_couches-1],k0,nu);
    E0=Ephi_prec_1-Ez_prec_1*Ephi2v-Hz_prec_1*Ephi4v;
    H0=Hphi_prec_1-Ez_prec_1*Hphi2v-Hz_prec_1*Hphi4v;
    E1=Ephi_prec_2-Ez_prec_2*Ephi2v-Hz_prec_2*Ephi4v;
    H1=Hphi_prec_2-Ez_prec_2*Hphi2v-Hz_prec_2*Hphi4v;
    e2=E1-E0;
    h2=H1-H0;
    e1=E0;
    h1=H0;
    alpha1=-e1/e2;
    alpha2=-h1/h2;
    alpha=(alpha1+alpha2)/2.0;

    A1=zeros(nb_couches);
    A2=zeros(nb_couches);
    A3=zeros(nb_couches);
    A4=zeros(nb_couches);
    Ez_prec=1.0;
    Hz_prec=alpha;
    Ephi_prec=Ephi1(index[1],radius[1],neff,radius[1],k0,nu)+alpha*Ephi3(index[1],radius[1],neff,radius[1],k0,nu);
    Hphi_prec=Hphi1(index[1],radius[1],neff,radius[1],k0,nu)+alpha*Hphi3(index[1],radius[1],neff,radius[1],k0,nu);
    A1[1]=1.0;
    A3[1]=alpha;
    @inbounds for couche=2:nb_couches-1
        m11=EHz1(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        m12=EHz2(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        m31=Ephi1(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        m32=Ephi2(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        m33=Ephi3(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        m34=Ephi4(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        m41=Hphi1(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        m42=Hphi2(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        m43=Hphi3(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        m44=Hphi4(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        m23=m11;
        m24=m12;
        det=1.0/(m11*m23*(m32*m44-m34*m42)-m11*m24*(m32*m43-m33*m42)-m12*m23*(m31*m44-m34*m41)+m12*m24*(m31*m43-m33*m41));
        sol1=det*((m23*m32*m44-m23*m34*m42-m24*m32*m43+m24*m33*m42)*Ez_prec+m12*(m33*m44-m34*m43)*Hz_prec-m12*(m23*m44-m24*m43)*Ephi_prec+m12*(m23*m34-m24*m33)*Hphi_prec);
		sol2=det*(-(m23*m31*m44-m23*m34*m41-m24*m31*m43+m24*m33*m41)*Ez_prec-m11*(m33*m44-m34*m43)*Hz_prec+m11*(m23*m44-m24*m43)*Ephi_prec-m11*(m23*m34-m24*m33)*Hphi_prec);
		sol3=det*(-m24*(m31*m42-m32*m41)*Ez_prec+(m11*m32*m44-m11*m34*m42-m12*m31*m44+m12*m34*m41)*Hz_prec+m24*(m11*m42-m12*m41)*Ephi_prec-m24*(m11*m32-m12*m31)*Hphi_prec);
		sol4=det*(m23*(m31*m42-m32*m41)*Ez_prec-(m11*m32*m43-m11*m33*m42-m12*m31*m43+m12*m33*m41)*Hz_prec-m23*(m11*m42-m12*m41)*Ephi_prec+m23*(m11*m32-m12*m31)*Hphi_prec);

        Ez_prec=sol1+sol2;
        Hz_prec=sol3+sol4;
        Ephi1v=Ephi1(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Ephi2v=Ephi2(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Ephi3v=Ephi3(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Ephi4v=Ephi4(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Hphi1v=Hphi1(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Hphi2v=Hphi2(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Hphi3v=Hphi3(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Hphi4v=Hphi4(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Ephi_prec=sol1*Ephi1v+sol2*Ephi2v+sol3*Ephi3v+sol4*Ephi4v;
        Hphi_prec=sol1*Hphi1v+sol2*Hphi2v+sol3*Hphi3v+sol4*Hphi4v;
        A1[couche]=sol1;
        A2[couche]=sol2;
        A3[couche]=sol3;
        A4[couche]=sol4;
    end
    A2[nb_couches]=Ez_prec;
    A4[nb_couches]=Hz_prec;

    return A1,A2,A3,A4;
end

function findroot_scalar(k0::Float64,nu::Int64,neff::Vector{Float64},radius::Union{Vector{Float64},Float64},index::Vector{Float64})
    nb_couches=length(index);
    phi_prec_task=Threads.@spawn phi1(index[1],radius[1],neff,radius[1],k0,nu);
    phi_prime_prec_task=Threads.@spawn phi1_prime(index[1],radius[1],neff,radius[1],k0,nu);
    phi_prec=fetch(phi_prec_task);
    phi_prime_prec=fetch(phi_prime_prec_task);
    @inbounds for couche=2:nb_couches-1
        p1_task=Threads.@spawn phi1(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        p1_prime_task=Threads.@spawn phi1_prime(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        p2_task=Threads.@spawn phi2(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        p2_prime_task=Threads.@spawn phi2_prime(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        p1=fetch(p1_task);
        p1_prime=fetch(p1_prime_task);
        p2=fetch(p2_task);
        p2_prime=fetch(p2_prime_task);
        det=@. 1.0/(p1*p2_prime-p2*p1_prime);
        A=@. det*(p2_prime*phi_prec-p2*phi_prime_prec);
        B=@. det*(-p1_prime*phi_prec+p1*phi_prime_prec);
        phi_prec=A+B;
        phi_prime_prec1_task=Threads.@spawn A.*phi1_prime(index[couche],radius[couche],neff,radius[couche],k0,nu);
        phi_prime_prec2_task=Threads.@spawn B.*phi2_prime(index[couche],radius[couche],neff,radius[couche],k0,nu);
        phi_prime_prec=fetch(phi_prime_prec1_task)+fetch(phi_prime_prec2_task);
    end
    phi_ext_task=Threads.@spawn phi2(index[nb_couches],radius[nb_couches-1],neff,radius[nb_couches-1],k0,nu);
    phi_prime_ext_task=Threads.@spawn phi2_prime(index[nb_couches],radius[nb_couches-1],neff,radius[nb_couches-1],k0,nu);
    phi_ext=fetch(phi_ext_task);
    phi_prime_ext=fetch(phi_prime_ext_task);
    delta=phi_ext.*phi_prime_prec-phi_prime_ext.*phi_prec;
    t=((delta[2:end-1].*delta[1:end-2]).<=0) .& (((delta[3:end]-delta[2:end-1]).*delta[2:end-1]).>=0);
    pos=findall(x->x>0,t).+1;
    return pos;
end

function findroot_TM(k0::Float64,nu::Int64,neff::Vector{Float64},radius::Union{Vector{Float64},Float64},index::Vector{Float64})
    nb_couches=length(index);
    Ez_prec=ones(size(neff));
    Hphi_prec=Hphi1(index[1],radius[1],neff,radius[1],k0,nu);
    @inbounds for couche=2:nb_couches-1
        a_task=Threads.@spawn EHz1(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        b_task=Threads.@spawn EHz2(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        c_task=Threads.@spawn Hphi1(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        d_task=Threads.@spawn Hphi2(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        a=fetch(a_task);
        b=fetch(b_task);
        c=fetch(c_task);
        d=fetch(d_task);
        det=@. 1.0/(a*d-b*c);
        A=@. det*(d*Ez_prec-b*Hphi_prec);
        B=@. det*(a*Hphi_prec-c*Ez_prec);
        Ez_prec=A+B;
        Hphi_prec1_task=Threads.@spawn A.*Hphi1(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Hphi_prec2_task=Threads.@spawn B.*Hphi2(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Hphi_prec=fetch(Hphi_prec1_task)+fetch(Hphi_prec2_task);
    end
    delta=Hphi_prec-Hphi2(index[nb_couches],radius[nb_couches-1],neff,radius[nb_couches-1],k0,nu).*Ez_prec;
    t=((delta[2:end-1].*delta[1:end-2]).<=0) .& (((delta[3:end]-delta[2:end-1]).*delta[2:end-1]).>=0);
    pos=findall(x->x>0,t).+1;
    return pos;
end

function findroot_TE(k0::Float64,nu::Int64,neff::Vector{Float64},radius::Union{Vector{Float64},Float64},index::Vector{Float64})
    nb_couches=length(index);
    Hz_prec=ones(size(neff));
    Ephi_prec=Ephi3(index[1],radius[1],neff,radius[1],k0,nu);
    @inbounds for couche=2:nb_couches-1
        a_task=Threads.@spawn EHz1(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        b_task=Threads.@spawn EHz2(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        c_task=Threads.@spawn Ephi3(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        d_task=Threads.@spawn Ephi4(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        a=fetch(a_task);
        b=fetch(b_task);
        c=fetch(c_task);
        d=fetch(d_task);
        det=@. 1.0/(a*d-b*c);
        A=@. det*(d*Hz_prec-b*Ephi_prec);
        B=@. det*(a*Ephi_prec-c*Hz_prec);
        Hz_prec=A+B;
        Ephi_prec1_task=Threads.@spawn A.*Ephi3(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Ephi_prec2_task=Threads.@spawn B.*Ephi4(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Ephi_prec=fetch(Ephi_prec1_task)+fetch(Ephi_prec2_task);
    end
    delta=Ephi_prec-Ephi4(index[nb_couches],radius[nb_couches-1],neff,radius[nb_couches-1],k0,nu).*Hz_prec;
    t=((delta[2:end-1].*delta[1:end-2]).<=0) .& (((delta[3:end]-delta[2:end-1]).*delta[2:end-1]).>=0);
    pos=findall(x->x>0,t).+1;
    return pos;
end

function findroot_HE_EH(k0::Float64,nu::Int64,neff::Vector{Float64},radius::Union{Vector{Float64},Float64},index::Vector{Float64})
    nb_couches=length(index);
    Ez_prec_1=ones(size(neff));
    Hz_prec_1=zeros(size(neff));
    Ephi_prec_1=Ephi1(index[1],radius[1],neff,radius[1],k0,nu);
    Hphi_prec_1=Hphi1(index[1],radius[1],neff,radius[1],k0,nu);
    Ez_prec_2=ones(size(neff));
    Hz_prec_2=ones(size(neff));
    Ephi_prec_2=Ephi_prec_1+Ephi3(index[1],radius[1],neff,radius[1],k0,nu);
    Hphi_prec_2=Hphi_prec_1+Hphi3(index[1],radius[1],neff,radius[1],k0,nu);
    @inbounds for couche=2:nb_couches-1
        m11_task=Threads.@spawn EHz1(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        m12_task=Threads.@spawn EHz2(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        m31_task=Threads.@spawn Ephi1(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        m32_task=Threads.@spawn Ephi2(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        m33_task=Threads.@spawn Ephi3(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        m34_task=Threads.@spawn Ephi4(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        m41_task=Threads.@spawn Hphi1(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        m42_task=Threads.@spawn Hphi2(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        m43_task=Threads.@spawn Hphi3(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        m44_task=Threads.@spawn Hphi4(index[couche],radius[couche],neff,radius[couche-1],k0,nu);
        m11=fetch(m11_task);
        m12=fetch(m12_task);
        m23=m11;
        m24=m12;
        m31=fetch(m31_task);
        m32=fetch(m32_task);
        m33=fetch(m33_task);
        m34=fetch(m34_task);
        m41=fetch(m41_task);
        m42=fetch(m42_task);
        m43=fetch(m43_task);
        m44=fetch(m44_task);
        det=@. 1.0/(m11*m23*(m32*m44-m34*m42)-m11*m24*(m32*m43-m33*m42)-m12*m23*(m31*m44-m34*m41)+m12*m24*(m31*m43-m33*m41));
        sol1_1=@. det*((m23*m32*m44-m23*m34*m42-m24*m32*m43+m24*m33*m42)*Ez_prec_1+m12*(m33*m44-m34*m43)*Hz_prec_1-m12*(m23*m44-m24*m43)*Ephi_prec_1+m12*(m23*m34-m24*m33)*Hphi_prec_1);
		sol2_1=@. det*(-(m23*m31*m44-m23*m34*m41-m24*m31*m43+m24*m33*m41)*Ez_prec_1-m11*(m33*m44-m34*m43)*Hz_prec_1+m11*(m23*m44-m24*m43)*Ephi_prec_1-m11*(m23*m34-m24*m33)*Hphi_prec_1);
		sol3_1=@. det*(-m24*(m31*m42-m32*m41)*Ez_prec_1+(m11*m32*m44-m11*m34*m42-m12*m31*m44+m12*m34*m41)*Hz_prec_1+m24*(m11*m42-m12*m41)*Ephi_prec_1-m24*(m11*m32-m12*m31)*Hphi_prec_1);
		sol4_1=@. det*(m23*(m31*m42-m32*m41)*Ez_prec_1-(m11*m32*m43-m11*m33*m42-m12*m31*m43+m12*m33*m41)*Hz_prec_1-m23*(m11*m42-m12*m41)*Ephi_prec_1+m23*(m11*m32-m12*m31)*Hphi_prec_1);
		sol1_2=@. det*((m23*m32*m44-m23*m34*m42-m24*m32*m43+m24*m33*m42)*Ez_prec_2+m12*(m33*m44-m34*m43)*Hz_prec_2-m12*(m23*m44-m24*m43)*Ephi_prec_2+m12*(m23*m34-m24*m33)*Hphi_prec_2);
		sol2_2=@. det*(-(m23*m31*m44-m23*m34*m41-m24*m31*m43+m24*m33*m41)*Ez_prec_2-m11*(m33*m44-m34*m43)*Hz_prec_2+m11*(m23*m44-m24*m43)*Ephi_prec_2-m11*(m23*m34-m24*m33)*Hphi_prec_2);
		sol3_2=@. det*(-m24*(m31*m42-m32*m41)*Ez_prec_2+(m11*m32*m44-m11*m34*m42-m12*m31*m44+m12*m34*m41)*Hz_prec_2+m24*(m11*m42-m12*m41)*Ephi_prec_2-m24*(m11*m32-m12*m31)*Hphi_prec_2);
		sol4_2=@. det*(m23*(m31*m42-m32*m41)*Ez_prec_2-(m11*m32*m43-m11*m33*m42-m12*m31*m43+m12*m33*m41)*Hz_prec_2-m23*(m11*m42-m12*m41)*Ephi_prec_2+m23*(m11*m32-m12*m31)*Hphi_prec_2);
        Ez_prec_1=sol1_1+sol2_1;
        Hz_prec_1=sol3_1+sol4_1;
        Ephi1_task=Threads.@spawn Ephi1(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Ephi2_task=Threads.@spawn Ephi2(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Ephi3_task=Threads.@spawn Ephi3(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Ephi4_task=Threads.@spawn Ephi4(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Hphi1_task=Threads.@spawn Hphi1(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Hphi2_task=Threads.@spawn Hphi2(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Hphi3_task=Threads.@spawn Hphi3(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Hphi4_task=Threads.@spawn Hphi4(index[couche],radius[couche],neff,radius[couche],k0,nu);
        Ephi_prec_1=sol1_1.*fetch(Ephi1_task)+sol2_1.*fetch(Ephi2_task)+sol3_1.*fetch(Ephi3_task)+sol4_1.*fetch(Ephi4_task);
        Hphi_prec_1=sol1_1.*fetch(Hphi1_task)+sol2_1.*fetch(Hphi2_task)+sol3_1.*fetch(Hphi3_task)+sol4_1.*fetch(Hphi4_task);
        Ez_prec_2=sol1_2+sol2_2;
        Hz_prec_2=sol3_2+sol4_2;
        Ephi_prec_2=sol1_2.*fetch(Ephi1_task)+sol2_2.*fetch(Ephi2_task)+sol3_2.*fetch(Ephi3_task)+sol4_2.*fetch(Ephi4_task);
        Hphi_prec_2=sol1_2.*fetch(Hphi1_task)+sol2_2.*fetch(Hphi2_task)+sol3_2.*fetch(Hphi3_task)+sol4_2.*fetch(Hphi4_task);
    end
    Ephi2_task=Threads.@spawn Ephi2(index[nb_couches],radius[nb_couches-1],neff,radius[nb_couches-1],k0,nu);
    Ephi4_task=Threads.@spawn Ephi4(index[nb_couches],radius[nb_couches-1],neff,radius[nb_couches-1],k0,nu);
    Hphi2_task=Threads.@spawn Hphi2(index[nb_couches],radius[nb_couches-1],neff,radius[nb_couches-1],k0,nu);
    Hphi4_task=Threads.@spawn Hphi4(index[nb_couches],radius[nb_couches-1],neff,radius[nb_couches-1],k0,nu);
    E0=Ephi_prec_1-Ez_prec_1.*fetch(Ephi2_task)-Hz_prec_1.*fetch(Ephi4_task);
    H0=Hphi_prec_1-Ez_prec_1.*fetch(Hphi2_task)-Hz_prec_1.*fetch(Hphi4_task);
    E1=Ephi_prec_2-Ez_prec_2.*fetch(Ephi2_task)-Hz_prec_2.*fetch(Ephi4_task);
    H1=Hphi_prec_2-Ez_prec_2.*fetch(Hphi2_task)-Hz_prec_2.*fetch(Hphi4_task);
    delta=E0.*H1-E1.*H0;
    t=((delta[2:end-1].*delta[1:end-2]).<=0) .& (((delta[3:end]-delta[2:end-1]).*delta[2:end-1]).>=0);
    pos=findall(x->x>0,t).+1;
    return pos;
end

function findroot(k0::Float64,nu::Int64,neff::Vector{Float64},radius::Union{Vector{Float64},Float64},index::Vector{Float64},type::Int64)
    if (type==1)
        return findroot_scalar(k0,nu,neff,radius,index);
    elseif (type==2)
        return findroot_TM(k0,nu,neff,radius,index);
    elseif (type==3)
        return findroot_TE(k0,nu,neff,radius,index);
    elseif (type==4)
        return findroot_HE_EH(k0,nu,neff,radius,index);
    else
        throw(ErrorException("Type must be between 1 and 4 in findroot"));
        return;
    end
end

#lambda : wavelength ; nu : azimuthal number ; radius : raddii vector ; index : refractive indices vector : maxPosition : maximum position vector for field computation ; numberofPoints : number of points between 0 and maxposition ; precision : precision on effective index ; type : Scalar/Vector ; firstDivision : number of divsion between min(n) and max(n) for the first step
"""
    multi_step_fiber_modes(lambda::Real,nu::Integer,radius::Union{Vector{<:Real},Real},index::Vector{<:Real};maxPosition::Real=0,numberofPoints::Integer=100,precision::Float64=1E-12,type::Symbol=:Scalar,firstDivision::Integer=0)

Returns a vector of `ScalarMode1D` if type=:Scalar or a vector of `VectorMode` if type=:Vector. If maxPosition=0, the modes do not include the electromagnetic field.

- lambda: wavelength
- nu: azimuthal number
- radius: outer radius of each layer (the cladding is inifinite and has no radius)
- index: refractive index of each layer (the cladding is included so that length(index) must be equal to length(radius)+1)
- maxPosition: maximum value of the radial coordinate r for scalar modes or cartesian coordinates x/y for vector modes. If maxPosition=0, the field is not computed.
- numberofPoints: number of points between 0 and maxPosition
- precision: precision required on the effective index 
- type: must be :Scalar or :Vector
- firstDivision: in the case of a very multimode fiber, you can increase this number if some modes are missing. If firstDivision=0, the value of firstDivision is approximated by the solver.
"""
function multi_step_fiber_modes(lambda::Real,nu::Integer,radius::Union{Vector{<:Real},Real},index::Vector{<:Real};maxPosition::Real=0,numberofPoints::Integer=100,precision::Float64=1E-12,type::Symbol=:Scalar,firstDivision::Integer=0)
    if (length(index) != (length(radius)+1))
        throw(DimensionMismatch("dim(radii) must be equal to dim(n)-1"));
    end
    if (lambda<=0)
        throw(DomainError(lambda, "Wavelength must be strictly positive"));
    end
    if (nu<0)
        throw(DomainError(nu, "The azimuthal number must be positive"));
    end
    if (minimum(radius)<=0)
        throw(DomainError(radius, "Radii must be positive"));
    end
    if (minimum(index)<=0)
        throw(DomainError(index, "Refractive index must be positive"));
    end
    if (precision<=0)
        throw(DomainError(precision, "Precision must be positive"));
    end
    if (firstDivision<0)
        throw(DomainError(firstDivision, "First precision must be positive"));
    end
    if (type!=:Scalar) & (type!=:Vector)
        throw(DomainError(firstDivision, "type must be :Scalar or :Vector"));
    end
    if (maxPosition<0)
        throw(DomainError(maxPosition, "The maximum position must be positive"));
    end
    if (numberofPoints<2)
        throw(DomainError(numberofPoints, "The number of points must be greater than 2"));
    end
    if (maxPosition>0)
        if (type==:Scalar)
            position=collect(LinRange(0,maxPosition,numberofPoints));
        else
            position=collect(LinRange(-maxPosition,maxPosition,2*numberofPoints-1));
        end
    else
        position=Vector{Float64}();
    end
    radius=convert.(Float64,radius);
    index=convert.(Float64,index);
    k0=2*pi/lambda;
    n_max=maximum(index);
    n_min=index[end];
    V=k0*radius[end]*sqrt(maximum(index)^2-index[end]^2);

    if (type==:Scalar)
        modes=Vector{ScalarMode1D}(undef,0);
    else
        modes=Vector{VectorMode}(undef,0);
    end

    if (type==:Scalar)
        t=1;
    elseif ((type==:Vector) & (nu==0))
        t=[2,3];
    elseif ((type==:Vector) & (nu>0))
        t=4;
    end

    pos0=0;
    for it=1:length(t)
        if (firstDivision==0)
            first=trunc(Int,maximum([V^2*10,100]));
            n1=length(findroot(k0,nu,collect(LinRange(n_min,n_max,first)),radius,index,t[it]));
            first=first*4;
            n2=length(findroot(k0,nu,collect(LinRange(n_min,n_max,first)),radius,index,t[it]));
            first=first*4;
            neff=collect(LinRange(n_min,n_max,first));
            pos=findroot(k0,nu,collect(LinRange(n_min,n_max,first)),radius,index,t[it]);
            n3=length(pos);
            while ((n1!=n2) & (n1!=n3))
                n1=n2;
                n2=n3;
                first=first*4;
                neff=collect(LinRange(n_min,n_max,first));
                pos=findroot(k0,nu,collect(LinRange(n_min,n_max,first)),radius,index,t[it]);
                n3=length(pos);
            end
        else
            first=maximum([firstDivision,100]);
            neff=collect(LinRange(n_min,n_max,first));
            pos=findroot(k0,nu,collect(LinRange(n_min,n_max,first)),radius,index,t[it]);
        end
        p=(n_max-n_min)/(first-1);
        neff_min=neff[pos.-1];
        neff_max=neff[pos.+1];
        if (!isempty(pos))
            while (p>precision)
                nb1=minimum([2001,trunc(Int,2*p/precision)+2]);
                p=(neff_max[1]-neff_min[1])/(nb1-1);
                for n_neff=1:length(neff_min)
                    neff=collect(LinRange(neff_min[n_neff],neff_max[n_neff],nb1));
                    pos=findroot(k0,nu,neff,radius,index,t[it]);
                    if (length(pos)!=1)
                        throw(ErrorException("Cannot find one root, try to increase firstDivision"));
                        return;
                    end
                    neff_min[n_neff]=neff[pos[1].-1];
                    neff_max[n_neff]=neff[pos[1].+1];
                end
            end
        end

        if (length(position)==0)
            Nm=length(neff_min);
            resize!(modes,pos0+Nm);
            for n=1:length(neff_min)
                neff_result=(neff_min[n]+neff_max[n])/2.0;
                if (t[it]==1)
                    name=string("LP ",string(nu),",",string(length(neff_min)-n+1));
                elseif (t[it]==2)
                    name=string("TM 0",",",string(length(neff_min)-n+1));
                elseif (t[it]==3)
                    name=string("TE 0",",",string(length(neff_min)-n+1));
                else
                    nb=length(neff_min)-n+1;
                    if (nb%2==1)
                        nb2=trunc(Int,(nb+1)/2);
                        name=string("HE ",string(nu),",",string(nb2));
                    else
                        nb2=trunc(Int,nb/2);
                        name=string("EH ",string(nu),",",string(nb2));
                    end
                end
                if (type==:Scalar)
                    modes[pos0+n]=ScalarMode1D(name,neff_result,lambda,nu,[],[]);
                else
                    empty_array=Array{Float64}(undef,0,2);
                    modes[pos0+n]=VectorMode(name,neff_result,lambda,[],[],empty_array,empty_array,empty_array,empty_array,empty_array,empty_array);
                end
            end
        else
            if (type==:Scalar)
                Nm=length(neff_min);
                resize!(modes,Nm);
                Threads.@threads for n=1:length(neff_min)
                    neff_result=(neff_min[n]+neff_max[n])/2.0;
                    name=string("LP ",string(nu),",",string(length(neff_min)-n+1));
                    A,B=findAB_scalar(k0,nu,neff_result,radius,index);
                    E=scalarField(position,nu,neff_result,k0,radius,index,A,B);
                    modes[n]=ScalarMode1D(name,neff_result,lambda,nu,position,E);
                end
            else
                #Vector case
                phi=atan.(position',position);
                r=hypot.(position,position');
                if (t[it]==2)
                    #TM
                    Nm=length(neff_min);
                    resize!(modes,Nm);
                    Threads.@threads for n=1:length(neff_min)
                        neff_result=(neff_min[n]+neff_max[n])/2.0;
                        A1,A2,A3,A4=findA_TM(k0,nu,neff_result,radius,index);
                        #=TM_Ex,TM_Ey=ExEy(r,phi,nu,neff_result,k0,radius,index,A1,A2,A3,A4,false);
                        TM_Ez,TM_Hz=EzHz(r,phi,nu,neff_result,k0,radius,index,A1,A2,A3,A4,false);
                        TM_Hx,TM_Hy=HxHy(r,phi,nu,neff_result,k0,radius,index,A1,A2,A3,A4,false);=#
                        TM_Ex_Ey_task=Threads.@spawn ExEy(r,phi,nu,neff_result,k0,radius,index,A1,A2,A3,A4,false);
                        TM_Ez_Hz_task=Threads.@spawn EzHz(r,phi,nu,neff_result,k0,radius,index,A1,A2,A3,A4,false);
                        TM_Hx_Hy_task=Threads.@spawn HxHy(r,phi,nu,neff_result,k0,radius,index,A1,A2,A3,A4,false);
                        TM_Ex,TM_Ey=fetch(TM_Ex_Ey_task);
                        TM_Ez,TM_Hz=fetch(TM_Ez_Hz_task);
                        TM_Hx,TM_Hy=fetch(TM_Hx_Hy_task);
                        name=string("TM 0,",string(length(neff_min)-n+1));
                        modes[n]=VectorMode(name,neff_result,lambda,position,position,TM_Ex,TM_Ey,-TM_Ez,TM_Hx,TM_Hy,-TM_Hz);
                    end
                elseif t[it]==3
                    #TE
                    Nm=length(neff_min);
                    resize!(modes,pos0+Nm);
                    Threads.@threads for n=1:length(neff_min)
                        neff_result=(neff_min[n]+neff_max[n])/2.0;
                        A1,A2,A3,A4=findA_TE(k0,nu,neff_result,radius,index);
                        #=TE_Ex,TE_Ey=ExEy(r,phi,nu,neff_result,k0,radius,index,A1,A2,A3,A4,true);
                        TE_Ez,TE_Hz=EzHz(r,phi,nu,neff_result,k0,radius,index,A1,A2,A3,A4,true);
                        TE_Hx,TE_Hy=HxHy(r,phi,nu,neff_result,k0,radius,index,A1,A2,A3,A4,true);=#
                        TE_Ex_Ey_task=Threads.@spawn ExEy(r,phi,nu,neff_result,k0,radius,index,A1,A2,A3,A4,true);
                        TE_Ez_Hz_task=Threads.@spawn EzHz(r,phi,nu,neff_result,k0,radius,index,A1,A2,A3,A4,true);
                        TE_Hx_Hy_task=Threads.@spawn HxHy(r,phi,nu,neff_result,k0,radius,index,A1,A2,A3,A4,true);
                        TE_Ex,TE_Ey=fetch(TE_Ex_Ey_task);
                        TE_Ez,TE_Hz=fetch(TE_Ez_Hz_task);
                        TE_Hx,TE_Hy=fetch(TE_Hx_Hy_task);
                        name=string("TE 0,",string(length(neff_min)-n+1));
                        modes[pos0+n]=VectorMode(name,neff_result,lambda,position,position,TE_Ex,TE_Ey,-TE_Ez,TE_Hx,TE_Hy,-TE_Hz);
                    end
                else
                    #HE-EH
                    Nm=2*length(neff_min);
                    resize!(modes,Nm);
                    Threads.@threads for n=1:length(neff_min)
                        nb=length(neff_min)-n+1;
                        neff_result=(neff_min[n]+neff_max[n])/2.0;
                        A1,A2,A3,A4=findA_HE_EH(k0,nu,neff_result,radius,index);
                        if (nb%2==1)
                            nb2=trunc(Int,(nb+1)/2);
                            name1=string("HE ",string(nu),",",string(nb2),'a');
                            name2=string("HE ",string(nu),",",string(nb2),'b');
                        else
                            nb2=trunc(Int,nb/2);
                            name1=string("EH ",string(nu),",",string(nb2),'a');
                            name2=string("EH ",string(nu),",",string(nb2),'b');
                        end
                        HEEH_Ex_Ey_task=Threads.@spawn ExEy(r,phi,nu,neff_result,k0,radius,index,A1,A2,A3,A4,true);
                        HEEH_Ez_Hz_task=Threads.@spawn EzHz(r,phi,nu,neff_result,k0,radius,index,A1,A2,A3,A4,true);
                        HEEH_Hx_Hy_task=Threads.@spawn HxHy(r,phi,nu,neff_result,k0,radius,index,A1,A2,A3,A4,true);
                        HEEH_Ex2_Ey2_task=Threads.@spawn ExEy(r,phi,nu,neff_result,k0,radius,index,A1,A2,A3,A4,false);
                        HEEH_Ez2_Hz2_task=Threads.@spawn EzHz(r,phi,nu,neff_result,k0,radius,index,A1,A2,A3,A4,false);
                        HEEH_Hx2_Hy2_task=Threads.@spawn HxHy(r,phi,nu,neff_result,k0,radius,index,A1,A2,A3,A4,false);
                        HEEH_Ex,HEEH_Ey=fetch(HEEH_Ex_Ey_task);
                        HEEH_Ez,HEEH_Hz=fetch(HEEH_Ez_Hz_task);
                        HEEH_Hx,HEEH_Hy=fetch(HEEH_Hx_Hy_task);
                        HEEH_Ex2,HEEH_Ey2=fetch(HEEH_Ex2_Ey2_task);
                        HEEH_Ez2,HEEH_Hz2=fetch(HEEH_Ez2_Hz2_task);
                        HEEH_Hx2,HEEH_Hy2=fetch(HEEH_Hx2_Hy2_task);
                        modes[2*n-1]=VectorMode(name1,neff_result,lambda,position,position,HEEH_Ex,HEEH_Ey,-HEEH_Ez,HEEH_Hx,HEEH_Hy,-HEEH_Hz);
                        modes[2*n]=VectorMode(name2,neff_result,lambda,position,position,HEEH_Ex2,HEEH_Ey2,-HEEH_Ez2,HEEH_Hx2,HEEH_Hy2,-HEEH_Hz2);
                    end
                end
            end
        end
        pos0=Nm;
    end
    return reverse(modes);
end
