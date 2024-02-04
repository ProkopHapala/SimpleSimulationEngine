# Load the essential modules
#using LinearAlgebra
#using SparseArrays
#using Plots

# =========== Functions

function Jacobi_step_Dens( A::Array{Float64,2}, b::Array{Float64,2}, x::Array{Float64,2}, xnew )
    n = size(b,1)
    #println("n : ", n)
    s   = zeros(3)
    err = 0.0
    for i=1:n
        s[:] .= 0.0
        for j=1:n
            if( j != i )
                s += x[j,:] .* A[i,j]
            end
        end
        y   = ( b[i,:] - s[:] ) ./ A[i,i]
        e = x[i,:]-y[:]
        err += sum( e.*e )
        xnew[i,:] = y
    end
    return err
end

function GaussSeidel_step_Sparse( A::Array{Float64,2}, b::Array{Float64,2}, x::Array{Float64,2} )
    n = length(b)
    s = zeros(3)
    err = 0.0
    for i=1:n
        s[:] = 0.0
        for j=1:n
            if( j != i )
                s += x[j,:] .* A[i,j]
            end
        end
        y   = ( b[i,:] - s[:] ) ./ A[i,i]
        e = x[i,:]-y[:]
        err += sum( e.*e )
        x[i,:] = y
    end
    return err
end

function SolveIterative( A::Array{Float64,2}, b::Array{Float64,2}, x::Array{Float64,2}; niter::Int=10, tol::Float64=1.e-3, method::Int=1 )
    xnew = zeros( size(x) )
    for iter=1:niter
        if     method == 1
            err = Jacobi_step_Dens( A, b, x, xnew )
            x[:,:] = xnew[:,:]
        elseif method == 2
            err = GaussSeidel_step_Sparse( A, b, x )
        end
        if err < tol^2
            println("Converged after ", iter, " iterations |err|=", sqrt(err) )
            return xnew
        else
            println("iter, err : ", iter, " ", sqrt(err) )
        end
        x[:] = xnew
    end 
end

function CholeskyDecomp( A::Array{Float64,2} )
    # seems to be: The Cholesky–Banachiewicz algorithm  (https://en.wikipedia.org/wiki/Cholesky_decomposition)
    n = size(A,1)
    L = zeros(n,n)
    for i in 1:n
        for j in 1:i
            sum1 = 0;
            if (j == i)  # diagonal
                for k in 1:j
                    sum1 += (L[j,k])^2
                end
                L[j,j] = sqrt(  A[j,j] - sum1 )
            else
                # Evaluating L(i, j)  using L(j, j)
                for k in 1:j
                    sum1 += (L[i,k] * L[j,k]);
                end
                if( L[j,j] > 0 )
                    L[i,j] = ( A[i,j] - sum1) / L[j,j]
                end
            end
        end
    end
    return L
end

function CholeskyDecomp_Crout( A::Array{Float64,2} )
    # Cholesky–Crout algorithm
    n = size(A,1)
    L = zeros(n,n)
    nop = 0 
    for j in 1:n
        sum1 = 0.0

        ngs = []
        for k in 1:j
            # put j,j inside string
            #if(L[j,k]>0) println( "[$j,$k] $(L[j,k]) " ) end
            if(L[j,k]>0) push!(ngs, k ) end
            sum1 += (L[j,k])^2
            nop += 1
        end
        println( "ngs[$j] n=", length(ngs), ngs )
        L[j,j] = sqrt(  A[j,j] - sum1 )

        for i in j+1:n
            sum1 = 0.0
            for k in 1:j
                sum1 += (L[i,k] * L[j,k]);
                nop += 1
            end
            L[i,j] = ( A[i,j] - sum1) / L[j,j]
        end
    end
    #println("CholeskyDecomp_Crout() nops=", nop );
    return L
end

function CholeskyDecomp_sparse( A::Array{Float64,2}, neighs::Array{Vector{Int},1}, neighs2::Array{Vector{Int},1}=nothing )
    n = size(A,1)
    L = zeros(n,n)
    nop = 0
    for j in 1:n
        sum1 = 0.0
        ngs = neighs[j]
        ngs_ = []
        for k in ngs
            if k <= j
                #println( "[$j,$k] $(L[j,k]) " )
                push!(ngs, k )
                sum1 += (L[j,k])^2
                nop += 1
            end
        end
        println( "sparse ngs[$j] n=", length(ngs), ngs )
        L[j,j] = sqrt(  A[j,j] - sum1 )

        for i in j+1:n
            sum1 = 0.0
            for k in 1:j
                sum1 += (L[i,k] * L[j,k]);
                nop += 1
            end
            L[i,j] = ( A[i,j] - sum1) / L[j,j]
        end

        # if neighs2 === nothing
        #     for i in j+1:n
        #         sum1 = 0.0
        #         for k in ngs
        #             if k <= j
        #                 sum1 += (L[i,k] * L[j,k]);
        #                 nop += 1
        #             end
        #         end
        #         L[i,j] = ( A[i,j] - sum1) / L[j,j]
        #     end
        # else
        #     ngs2 = neighs2[j]
        #     for i in ngs2
        #         if i > j
        #             sum1 = 0.0
        #             for k in ngs
        #                 if k <= j
        #                     sum1 += (L[i,k] * L[j,k]);
        #                     nop += 1
        #                 end
        #             end
        #             L[i,j] = ( A[i,j] - sum1) / L[j,j]
        #         end
        #     end
        # end
    end

    #println("CholeskyDecomp_sparse() nops=", nop );
    return L
end

# https://fncbook.github.io/fnc/linsys/linear-systems.html

function forwardsub(L::Array{Float64,2},b::Array{Float64,2})
    n = size(L,1)
    x = zeros( size(b)   )
    s = zeros( size(b,2) )
    x[1,:] = b[1,:]/L[1,1]
    for i = 2:n
        s[:] .= 0.0
        for j = 1:i-1
            s[:] .+= L[i,j]*x[j,:]
        end
        x[i,:] = ( b[i,:] .- s[:] ) ./ L[i,i]
    end    
    return x
end

function backsub(U::Array{Float64,2},b::Array{Float64,2})
    n = size(U,1)
    x = zeros( size(b)   )
    s = zeros( size(b,2) )
    x[n] = b[n]/U[n,n]
    for i = n-1:-1:1
        s[:] .= 0.0
        for j=i+1:n
            s[:] .+= U[i,j]*x[j,:]
        end
        x[i,:] = ( b[i] .- s[:] ) ./ U[i,i]
    end    
    return x
end