module OpticalFibersMumpsExt

    using OpticalFibers
    using LinearMaps
    using LinearAlgebra
    using SparseArrays
    using ArnoldiMethod
    using MPI
    using MUMPS

    struct ShiftAndInvert_MUMPS
        m::Mumps
        B
        temp
    end

    function (M::ShiftAndInvert_MUMPS)(y,x)
        MUMPS.set_job!(M.m,4)
        mul!(M.temp,M.B,x)
        associate_rhs!(M.m,M.temp)
        MUMPS.solve!(M.m)
        MUMPS.get_rhs!(y,M.m)
    end

    function OpticalFibers.eigs_MUMPS(A,B;sigma=0,nev::Int64=1,tol::Float64=0.0,restarts::Int64=200,verbose::Bool=false)
        if (!MPI.Initialized())
            MPI.Init();
        end
        icntl = get_icntl(verbose=false)
        m = Mumps{eltype(A)}(mumps_unsymmetric, icntl, default_cntl64)
        associate_matrix!(m,SparseMatrixCSC(A-sigma*B));
        factorize!(m);
        a = ShiftAndInvert_MUMPS(m,B,Vector{eltype(B)}(undef, size(B,1)))
        map_mumps=LinearMap{eltype(A)}(a, size(A,1), ismutating=true)
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

end