module OpticalFibersMumpsExt

    using OpticalFibers
    using LinearMaps
    using LinearAlgebra
    using SparseArrays
    using ArnoldiMethod
    using MPI
    using MUMPS

    struct ShiftAndInvert_MUMPS{TA,TB,TT}
        m::TA
        B::TB
        temp::TT
    end

    function (M::ShiftAndInvert_MUMPS)(y,x)
        MUMPS.set_job!(M.m,4)
        mul!(M.temp,M.B,x)
        associate_rhs!(M.m,M.temp)
        MUMPS.solve!(M.m)
        MUMPS.get_rhs!(y,M.m)
    end

    function OpticalFibers.eigs_MUMPS(A,B;sigma=0,nev::Int64=1,tol::Float64=0.0,restarts::Int64=200,verbose::Bool=false)
        N=size(A,1)
        if (!MPI.Initialized())
            MPI.Init();
        end
        icntl = get_icntl(verbose=false)
        m = Mumps{eltype(A)}(mumps_unsymmetric, icntl, default_cntl64)
        associate_matrix!(m,SparseMatrixCSC(A-sigma*B));
        factorize!(m);
        a = ShiftAndInvert_MUMPS(m,B,Vector{eltype(B)}(undef, N))
        map_mumps=LinearMap{eltype(A)}(a, N, ismutating=true)
        if (tol!=0)
            decomp, history  = partialschur(map_mumps, nev=nev, tol=tol, restarts=restarts, which=:LM)
        else
            decomp, history  = partialschur(map_mumps, nev=nev, restarts=restarts, which=:LM)
        end
        if (verbose)
            @show history
        end
        λs_inv, X = partialeigen(decomp);
        MUMPS.finalize!(m);
        λs=(1 ./λs_inv).+sigma
        return λs,X;
    end

    function OpticalFibers.eigs_MUMPS(A;sigma=0,nev::Int64=1,tol::Float64=0.0,restarts::Int64=200,verbose::Bool=false)
        OpticalFibers.eigs_MUMPS(A,eltype(A).(SparseMatrixCSC(I,size(A,1),size(A,2)));sigma=sigma,nev=nev,tol=tol,restarts=restarts,verbose=verbose)
    end

    struct ShiftAndInvert_quadratic_MUMPS{solver,matrix,vector}
        m::solver
        A::matrix
        C::matrix
        sigma::ComplexF64
        temp::vector
        temp2::vector
        N::Int64
    end

    function (M::ShiftAndInvert_quadratic_MUMPS)(y,x)
        MUMPS.set_job!(M.m,4)
        x1=@view x[1:M.N]
        x2=@view x[M.N+1:end]
        y1=@view y[1:M.N]
        y2=@view y[M.N+1:end]
        mul!(M.temp2, M.C, x2)
        mul!(M.temp, M.A, x1)
        axpy!(-M.sigma,M.temp2,M.temp)
        associate_rhs!(M.m,M.temp)
        MUMPS.solve!(M.m)
        MUMPS.get_rhs!(y2,M.m)
        @. y1=(y2-x1)/M.sigma
    end

    function OpticalFibers.eigs_quadratic_MUMPS(A,B,C;sigma=0,nev::Int64=1,tol::Float64=0.0,restarts::Int64=200,verbose::Bool=false)
        if (!MPI.Initialized())
            MPI.Init();
        end
        N=size(A,1)
        icntl = get_icntl(verbose=false)
        m = Mumps{eltype(A)}(mumps_unsymmetric, icntl, default_cntl64)
        associate_matrix!(m,SparseMatrixCSC(A+sigma*B+sigma^2*C));
        factorize!(m);
        a = ShiftAndInvert_quadratic_MUMPS(m,A,C,ComplexF64(sigma),Vector{eltype(B)}(undef, N),Vector{eltype(B)}(undef, N),N)
        map_mumps=LinearMap{eltype(A)}(a, 2*N, ismutating=true)
        if (tol!=0)
            decomp, history  = partialschur(map_mumps, nev=nev, tol=tol, restarts=restarts, which=:LM)
        else
            decomp, history  = partialschur(map_mumps, nev=nev, restarts=restarts, which=:LM)
        end
        if (verbose)
            @show history
        end
        λs_inv, X = partialeigen(decomp);
        MUMPS.finalize!(m);
        λs=(1 ./λs_inv).+sigma
        return λs,X;
    end

end