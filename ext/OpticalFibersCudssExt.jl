module OpticalFibersCudssExt

    using OpticalFibers
    using LinearMaps
    using LinearAlgebra
    using SparseArrays
    using ArnoldiMethod
    using CUDA,CUDA.CUSPARSE
    using CUDSS

    struct ShiftAndInvert_CUDA{gpu_solver,gpu_matrix,gpu_vector}
        M::gpu_matrix #stockée pour éviter les pb de mémoire (M passe par le GC alors qu'elle ne devrait pas)
        A_lu::gpu_solver
        B::gpu_matrix
        x_gpu::gpu_vector
        temp::gpu_vector
        N::Int64
    end

    function (M::ShiftAndInvert_CUDA)(y,x)
        copyto!(M.x_gpu,unsafe_wrap(Vector{eltype(x)}, pointer(x), M.N))
        mul!(M.temp, M.B, M.x_gpu) 
        cudss("solve", M.A_lu, M.x_gpu, M.temp)
        copyto!(unsafe_wrap(Vector{eltype(y)}, pointer(y), M.N),M.x_gpu)
    end

    function OpticalFibers.eigs_CUDA(A,B;sigma=0,nev::Int64=1,tol::Float64=0.0,restarts::Int64=200,verbose::Bool=false,ir_n_steps::Int64=10,matching_alg::String="algo6",ir_tol::Float64=1E-9)
        M=CuSparseMatrixCSR(A-sigma*B)
        N=size(A,1)
        solver=CudssSolver(M,"G",'F')
        cudss_set(solver, "matching_alg", matching_alg)
        cudss_set(solver,"ir_tol",ir_tol)
        cudss_set(solver,"ir_n_steps",ir_n_steps)
        x = CUDA.zeros(eltype(A),10)#sans importance pour les 2 phases suivantes
        cudss("analysis", solver, x, x)
        cudss("factorization", solver, x, x)
        a = ShiftAndInvert_CUDA(M,solver,CuSparseMatrixCSR(B),CuVector{eltype(A)}(undef, N),CuVector{eltype(A)}(undef, N),N)
        map_CUDA=LinearMap{eltype(A)}(a, N, ismutating=true);
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

    struct ShiftAndInvert_quadratic_CUDA{gpu_solver,gpu_matrix,gpu_vector,gpu_vector_view}
        M::gpu_matrix
        A_lu::gpu_solver
        A::gpu_matrix
        C::gpu_matrix
        sigma::ComplexF64
        x_gpu::gpu_vector
        temp_gpu::gpu_vector
        temp_gpu2::gpu_vector
        x1_gpu::gpu_vector_view
        x2_gpu::gpu_vector_view
        N::Int64
    end

    function ShiftAndInvert_quadratic_CUDA(M,A_lu,A,C,sigma,x_gpu,temp_gpu,temp_gpu2,N)
        x1_gpu=@view x_gpu[1:N]
        x2_gpu=@view x_gpu[N+1:end]
        return ShiftAndInvert_quadratic_CUDA(M,A_lu,A,C,sigma,x_gpu,temp_gpu,temp_gpu2,x1_gpu,x2_gpu,N)
    end

    function (M::ShiftAndInvert_quadratic_CUDA)(y,x)
        sigma=M.sigma
        copyto!(M.x_gpu,unsafe_wrap(Vector{eltype(x)}, pointer(x), 2*M.N))
        mul!(M.temp_gpu2, M.C, M.x2_gpu)
        mul!(M.temp_gpu, M.A, M.x1_gpu)
        axpy!(-sigma,M.temp_gpu2,M.temp_gpu)#A voir si ce n'est pas mieux de faire du @.
        #@.  M.temp_gpu .-= M.temp_gpu2 .*sigma
        cudss("solve", M.A_lu, M.x2_gpu, M.temp_gpu)
        @. M.x1_gpu=(M.x2_gpu-M.x1_gpu)/sigma
        copyto!(unsafe_wrap(Vector{eltype(y)}, pointer(y), 2*M.N),M.x_gpu)
    end

    function OpticalFibers.eigs_quadratic_CUDA(A,B,C;sigma=0,nev::Int64=1,tol::Float64=0.0,restarts::Int64=200,verbose::Bool=false,ir_n_steps::Int64=10,matching_alg::String="algo6",ir_tol::Float64=1E-9)
        N=size(A,1)
        #Mettre M dans la structure permet de faire en sorte qu'elle reste bien en mémoire durant tout le processus
        M=CuSparseMatrixCSR(A+sigma*B+sigma^2*C)
        solver=CudssSolver(M,"G",'F')
        cudss_set(solver, "matching_alg", matching_alg)
        cudss_set(solver,"ir_tol",ir_tol)
        cudss_set(solver,"ir_n_steps",ir_n_steps)
        x = CUDA.zeros(eltype(A),10)#sans importance pour les 2 phases suivantes
        cudss("analysis", solver, x, x)
        cudss("factorization", solver, x, x)
        A_gpu=CuSparseMatrixCSR(A)
        C_gpu=CuSparseMatrixCSR(C)
        a = ShiftAndInvert_quadratic_CUDA(M,solver,A_gpu,C_gpu,ComplexF64(sigma),CuVector{eltype(A)}(undef, 2*N),CuVector{eltype(A)}(undef, N),CuVector{eltype(A)}(undef, N),N)
        map_CUDA=LinearMap{eltype(A)}(a, 2*N, ismutating=true);
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

end