module OpticalFibersCudssExt

    using OpticalFibers
    using LinearMaps
    using LinearAlgebra
    using SparseArrays
    using ArnoldiMethod
    using CUDA,CUDA.CUSPARSE
    using CUDSS

    struct ShiftAndInvert_CUDA{TA,TB,TT,TT2}
        A_lu::TA
        B::TB
        temp::TT
        y2::TT
        temp_cpu::TT2
    end

    function (M::ShiftAndInvert_CUDA)(y,x)
        mul!(M.temp_cpu, M.B, x)
        copyto!(M.temp,M.temp_cpu)
        cudss("solve", M.A_lu, M.y2, M.temp)
        y.=Vector(M.y2)
    end

    function OpticalFibers.eigs_CUDA(A,B;sigma=0,nev::Int64=1,tol::Float64=0.0,restarts::Int64=200,verbose::Bool=false,ir_n_steps::Int64=10)
        solver=CudssSolver(CuSparseMatrixCSR(A-sigma*B),"G",'F')
        cudss_set(solver, "use_matching", 1)
        cudss_set(solver,"ir_n_steps",ir_n_steps)
        x = CUDA.zeros(eltype(A),10)#sans importance pour les 2 phases suivantes
        cudss("analysis", solver, x, x)
        cudss("factorization", solver, x, x)
        a = ShiftAndInvert_CUDA(solver,B,CuVector{eltype(A)}(undef, size(A,1)),CuVector{eltype(A)}(undef, size(A,1)),Vector{eltype(A)}(undef, size(A,1)))
        map_CUDA=LinearMap{eltype(A)}(a, size(A,1), ismutating=true);
        if (tol!=0)
            decomp,history  = partialschur(map_CUDA, nev=nev, tol=tol, restarts=restarts, which=:LM)
        else
            decomp,history  = partialschur(map_CUDA, nev=nev, restarts=restarts, which=:LM)
        end
        if (verbose)
            @show history
        end
        λs_inv, X = partialeigen(decomp);
        λs=(1 ./λs_inv).+sigma
        return λs,X;
    end

    function OpticalFibers.eigs_CUDA(A;sigma=0,nev::Int64=1,tol::Float64=0.0,restarts::Int64=200,verbose::Bool=false,ir_n_steps::Int64=10)
        OpticalFibers.eigs_CUDA(A,eltype(A).(SparseMatrixCSC(I,size(A,1),size(A,2)));sigma=sigma,nev=nev,tol=tol,restarts=restarts,verbose=verbose,ir_n_steps=ir_n_steps)
    end

end