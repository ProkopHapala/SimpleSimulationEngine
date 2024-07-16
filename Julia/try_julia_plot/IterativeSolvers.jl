# Load the essential modules
#using LinearAlgebra
#using SparseArrays
#using Plots

using DataStructures

# =========== Functions

function Jacobi_step_Dens( A::Array{Float64,2}, b::Array{Float64,2}, x::Array{Float64,2}, xnew::Array{Float64,2} )
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
    #println( "CholeskyDecomp_Crout()" )
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
        #println( "ngs[$j] ", ngs )
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

function CholeskyDecomp_sparse( A::Array{Float64,2}, neighs::Array{Vector{Int},1}, neighs2::Array{Vector{Int},1}=nothing, tol::Float64=1.e-16, bSortedSets::Bool=true )
    #neighsets = [ Set(ng) for ng in neighs ] 
    # How can I initialize set with reserved size in Julia?
    #s = Set{Int}(100)  # set with reserved size 100
    # if bSortedSets
    #     neighsets = [ SortedSet(ng) for ng in neighs ]
    # else
    #     neighsets = [ Set(ng) for ng in neighs ] 
    # end
    #neighsets = [ SortedSet(ng) for ng in neighs ] 
    neighsets = [ Set(ng) for ng in neighs ]          # this seems faster
    n = size(A,1)
    L = zeros(n,n)
    nop = 0
    #println( "CholeskyDecomp_sparse()" )

    nngmax = 0
    for j in 1:n
        sum1 = 0.0
        ngs  = neighsets[j]
        ngs_ = []
        for k in ngs
            if k <= j
                #println( "[$j,$k] $(L[j,k]) " )
                #push!(ngs_, k )
                sum1 += (L[j,k])^2
                nop += 1
            end
        end
        #println( "ngs[$j] ", ngs_, "       ", ngs )
        Ljj    = sqrt(  A[j,j] - sum1 )
        L[j,j] = Ljj
        invLjj = 1.0/Ljj

        for i in j+1:n
            sum1 = 0.0
            for k in ngs
                if k <= j
                    sum1 += L[i,k] * L[j,k];
                    nop += 1
                end
            end
            Lij = ( A[i,j] - sum1) * invLjj
            if abs(Lij) > tol
                L[i,j] = Lij
                push!(ngs,          i )
                push!(neighsets[i], j )
            end
        end

        # for i in j+1:n
        #     sum1 = 0.0
        #     for k in ngs
        #         if k <= j
        #             sum1 += (L[i,k] * L[j,k]);
        #             nop += 1
        #         end
        #     end
        #     Lij = ( A[i,j] - sum1) / L[j,j]
        #     if abs(Lij) > tol
        #         L[i,j] = Lij
        #         push!(ngs, i )
        #         push!(neighsets[i], j )
        #     end
        # end

        # ngs2 = neighs2[j]
        # for i in ngs2
        #     if i > j
        #         sum1 = 0.0
        #         for k in ngs
        #             if k <= j
        #                 sum1 += (L[i,k] * L[j,k]);
        #                 nop += 1
        #             end
        #         end
        #         Lij = ( A[i,j] - sum1) / L[j,j]
        #         if abs(Lij) > tol
        #             L[i,j] = Lij
        #             push!(ngs, i )
        #             push!(neighsets[i], j )
        #         end
        #     end
        # end
        
        # for i in j+1:n
        #     sum1 = 0.0
        #     for k in 1:j
        #         sum1 += (L[i,k] * L[j,k]);
        #         nop += 1
        #     end
        #     Lij = ( A[i,j] - sum1) / L[j,j]
        #     if abs(Lij) > tol
        #         L[i,j] = Lij
        #         push!(ngs, i )
        #         push!(neighsets[i], j )
        #     end
        # end

        nngmax = max( nngmax, length(ngs) )
    end
    neighsCh = [ sort!(collect(ngs)) for ngs in neighsets ]
    #println("CholeskyDecomp_sparse() nops=", nop, " nngmax=", nngmax );
    return L, neighsCh
end

function IncompleteCholeskyDecomp(A::Array{Float64,2} )
    # https://en.wikipedia.org/wiki/Incomplete_Cholesky_factorization
	n = size(a,1);
	for k = 1:n
		A[k,k] = sqrt(A[k,k]);
		for i = (k+1):n
		    if (A[i,k] != 0)
		        A[i,k] = A[i,k]/A[k,k];            
		    end
		end
		for j = (k+1):n
		    for i = j:n
		        if (A[i,j] != 0)
		            A[i,j] = A[i,j] - A[i,k]*A[j,k];  
		        end
		    end
		end
	end
    # for i = 1:n
    #     for j = i+1:n
    #         a[i,j] = 0;
    #     end
    # end       
end

function IncompleteCholeskyDecompSparse(a::Array{Float64,2}, neighs::Array{Vector{Int},1}, tol = 1.e-10 )
    # https://en.wikipedia.org/wiki/Incomplete_Cholesky_factorization
	n = size(a,1);
	for k = 1:n
        Akk = sqrt(A[k,k]);
		A[k,k] = Akk
        invAkk = 1.0/Akk;
        ngk = neighs[k];
        for i in ngk
            if i > k
                aik = A[i,k];
                if ( abs(aik)>tol )
                    A[i,k] = aik*invAkk          
                end
            end
        end
        for j in ngk
            if j > k
                ajk = A[j,k]
                ngj = neighs[j];
                for i in ngj
                    if i >= j
                        aij = A[i,j];
                        if ( abs(aij)>tol )
                            A[i,j] = aij - A[i,k]*ajk;  
                        end
                    end
                end
            end
        end
	end
    # for i = 1:n
    #     for j = i+1:n
    #         a[i,j] = 0;
    #     end
    # end       
end

function mapMatrixNeighs( A::Array{Float64,2}, tol::Float64 = 1.e-16 )
    # find which rows i are coupled with witch columns j by matrix elements Aij 
    n = size(A,1)
    neighs  = [ Vector{Int}() for _ in 1:n ]   
    #neighsT = [ Vector{Int}() for _ in 1:n ]
    nngmax = 0
    for i in 1:n
        ngi = neighs[i]
        for j in 1:i-1
            if abs(A[i,j]) > tol
                push!(ngi,        j )
                #push!(neighsT[j], i )
            end
        end
        nngi   = length(ngi)
        nngmax = max( nngmax, nngi )
        #println("[$i] nngi : ", nngi )
    end
    #println("mapMatrixNeighs() nngmax : ", nngmax )
    return neighs #, neighsT
end

function invNeighs( neighs::Array{Vector{Int},1} )
    n = length(neighs)
    neighsT = [ Vector{Int}() for _ in 1:n ]
    for i in 1:n
        #ngT = neighsT[i]
        for j in neighs[i]
            push!(neighsT[j], i )
        end
    end
    return neighsT
end

# https://fncbook.github.io/fnc/linsys/linear-systems.html

function forwardsub(L::Array{Float64,2},b::Array{Float64,2})
    n = size(L,1)
    x = zeros( size(b)   )
    s = zeros( size(b,2) )
    x[1,:] = b[1,:]/L[1,1]
    nop = 0
    for i = 2:n
        s[:] .= 0.0
        for j = 1:i-1
            s[:] .+= L[i,j]*x[j,:]
            nop += 1
        end
        x[i,:] = ( b[i,:] .- s[:] ) ./ L[i,i]
    end    
    #println("forwardsub() nops=", nop );   
    return x
end

function forwardsub_sparse(L::Array{Float64,2}, b::Array{Float64,2}, neighs::Array{Vector{Int},1} )
    n = size(L,1)
    x = zeros( size(b)   )
    s = zeros( size(b,2) )
    x[1,:] = b[1,:]/L[1,1]
    nop = 0
    for i = 2:n
        s[:] .= 0.0
        for j in neighs[i]
            if j < i
                s[:] .+= L[i,j]*x[j,:]
                nop += 1
            end
        end
        x[i,:] = ( b[i,:] .- s[:] ) ./ L[i,i]
    end  
    #println("forwardsub_sparse() nops=", nop );     
    return x
end

function backsub_bad(U::Array{Float64,2},b::Array{Float64,2})
    n = size(U,1)
    x = zeros( size(b)   )
    s = zeros( size(b,2) )
    x[n] = b[n]/U[n,n]
    nop = 0
    for i = n-1:-1:1
        s[:] .= 0.0
        for j=i+1:n
            s[:] .+= U[i,j]*x[j,:]
            nop += 1
        end
        x[i,:] = ( b[i] .- s[:] ) ./ U[i,i]
    end    
    #println("backsub_sparse() nops=", nop );   
    return x
end

function backsub(U::Array{Float64,2}, b::Array{Float64,2})
    n = size(U, 1)
    x = zeros(size(b))
    s = zeros(size(b, 2))
    x[n, :] = b[n, :] / U[n, n]  # Corrected to element-wise operation
    nop = 0
    for i = n-1:-1:1
        s[:] .= 0.0
        for j = i+1:n
            s[:] .+= U[i, j] * x[j, :]
            nop += 1
        end
        x[i, :] = (b[i, :] .- s[:]) ./ U[i, i]  # Corrected to element-wise operation
    end
    return x
end

function backsub_sparse(U::Array{Float64,2},b::Array{Float64,2}, neighs::Array{Vector{Int},1} )
    n = size(U,1)
    x = zeros( size(b)   )
    s = zeros( size(b,2) )
    x[n] = b[n]/U[n,n]
    nop = 0
    for i = n-1:-1:1
        s[:] .= 0.0
        for j in neighs[i]
            if j > i
                s[:] .+= U[i,j]*x[j,:]
                nop += 1
            end
        end
        x[i,:] = ( b[i] .- s[:] ) ./ U[i,i]
    end 
    #println("backsub_sparse() nops=", nop );   
    return x
end

