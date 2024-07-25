# Load the essential modules
#using LinearAlgebra
#using SparseArrays
#using Plots

using DataStructures

# function write_matrix_to_file(A::Matrix{T}, filename::String; delimiter::Char='\t', format::String="%.6f") where T
#     open(filename, "w") do io
#         rows, cols = size(A)
#         for i in 1:rows
#             for j in 1:cols
#                 @printf(io, format, A[i,j])
#                 if j < cols
#                 print(io, delimiter)
#             end
#         end
#         println(io)
#         end
#     end
#     println("Matrix written to $filename")
# end

# =========== Functions

function Jacobi_step_Dens( A::Matrix{Float64}, b::Matrix{Float64}, x::Matrix{Float64}, xnew::Matrix{Float64} )
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

function GaussSeidel_step_Sparse( A::Matrix{Float64}, b::Matrix{Float64}, x::Matrix{Float64} )
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

function SolveIterative( A::Matrix{Float64}, b::Matrix{Float64}, x::Matrix{Float64}; niter::Int=10, tol::Float64=1.e-3, method::Int=1 )
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

function CholeskyDecomp( A::Matrix{Float64} )
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

function CholeskyDecomp_Crout( A::Matrix{Float64} )
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

function CholeskyDecomp_sparse( A::Matrix{Float64}, neighs::Array{Vector{Int},1}, neighs2::Array{Vector{Int},1}=nothing, tol::Float64=1.e-16, bSortedSets::Bool=true )
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

function IncompleteCholeskyDecomp(A::Matrix{Float64} )
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

function IncompleteCholeskyDecompSparse(a::Matrix{Float64}, neighs::Array{Vector{Int},1}, tol = 1.e-10 )
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

function mapMatrixNeighs( A::Matrix{Float64}, tol::Float64 = 1.e-16 )
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

function forwardsub(L::Matrix{Float64},b::Matrix{Float64})
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

function backsub(U::Matrix{Float64}, b::Matrix{Float64})
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

function forwardsub_sparse(L::Matrix{Float64}, b::Matrix{Float64}, neighs::Array{Vector{Int},1} )
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

function backsub_sparse(U::Matrix{Float64},b::Matrix{Float64}, neighs::Array{Vector{Int},1} )
    n = size(U,1)
    x = zeros( size(b)   )
    s = zeros( size(b,2) )
    #x[n]   = b[n]    / U[n, n] # this was probably error 
    x[n, :] = b[n, :] / U[n, n] # Correced 
    nop = 0
    for i = n-1:-1:1
        s[:] .= 0.0
        for j in neighs[i]
            if j > i
                s[:] .+= U[i,j]*x[j,:]
                nop += 1
            end
        end
        x[i,:] = ( b[i,:] .- s[:] ) ./ U[i,i]
    end 
    #println("backsub_sparse() nops=", nop );   
    return x
end




function forwardsub_sparse_f(L::Matrix{Float32}, b::Matrix{Float32}, neighs::Array{Vector{Int},1} )
    n = size(L,1)
    x = zeros( Float32, size(b)   )
    s = zeros( Float32, size(b,2) )
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

function backsub_sparse_f(U::Matrix{Float32},b::Matrix{Float32}, neighs::Array{Vector{Int},1} )
    n = size(U,1)
    x = zeros( Float32, size(b)   )
    s = zeros( Float32, size(b,2) )
    #x[n]   = b[n]    / U[n, n] # this was probably error 
    x[n, :] = b[n, :] / U[n, n] # Correced 
    nop = 0
    for i = n-1:-1:1
        s[:] .= 0.0
        for j in neighs[i]
            if j > i
                s[:] .+= U[i,j]*x[j,:]
                nop += 1
            end
        end
        x[i,:] = ( b[i,:] .- s[:] ) ./ U[i,i]
    end 
    #println("backsub_sparse() nops=", nop );   
    return x
end



function CholeskyDecomp_LDLT(A::Matrix{T}) where T<:AbstractFloat
    n = size(A, 1)
    L = Matrix{T}(I, n, n)  # Initialize L as identity matrix
    D = zeros(T, n)
    for j in 1:n
        D[j] = A[j,j] - dot(L[j,1:j-1].^2, D[1:j-1])
        for i in j+1:n
            L[i,j] = (A[i,j] - dot(L[i,1:j-1] .* L[j,1:j-1], D[1:j-1])) / D[j]
        end
    end
    return L, D
end

function forward_substitution(L::Matrix{T}, b::Vector{T}) where T<:AbstractFloat
    n = size(L, 1)
    y = zeros(T, n)
    for i in 1:n
        y[i] = b[i] - dot(L[i,1:i-1], y[1:i-1])
    end
    return y
end

function forward_substitution_transposed(L::Matrix{T}, b::Vector{T}) where T<:AbstractFloat
    n = size(L, 1)
    x = zeros(T, n)
    for i in n:-1:1
        x[i] = b[i] - dot(L[i+1:n,i], x[i+1:n])
    end
    return x
end

function solve_LDLT( L::Matrix{T}, D::Vector{T}, b::Vector{T}) where T<:AbstractFloat
    #L, D = CholeskyLDLDecomp(A)
    z = forward_substitution(L, b)
    #y = diagonal_solve(D, z)
    y = z ./ D
    x = forward_substitution_transposed(L, y)
    return x
end


function CholeskyDecomp_LDLT_sparse(A::Matrix{T}, neighs::Array{Vector{Int},1}, tol::T=1e-16) where T<:AbstractFloat
    neighsets = [Set(ng) for ng in neighs]
    n = size(A, 1)
    L = Matrix{T}(I, n, n)
    D = zeros(T, n)
    nop = 0
    nngmax = 0

    for j in 1:n
        sum1 = 0.0
        ngs = neighsets[j]
        for k in ngs
            if k < j
                sum1 += (L[j,k]^2) * D[k]
                nop += 1
            end
        end
        D[j] = A[j,j] - sum1

        for i in j+1:n
            sum1 = 0.0
            for k in ngs
                if k < j
                    sum1 += L[i,k] * L[j,k] * D[k]
                    nop += 1
                end
            end
            Lij = (A[i,j] - sum1) / D[j]
            if abs(Lij) > tol
                L[i,j] = Lij
                push!(ngs, i)
                push!(neighsets[i], j)
            end
        end
        nngmax = max(nngmax, length(ngs))
    end

    neighsLDLT = [sort!(collect(ngs)) for ngs in neighsets]
    return L, D, neighsLDLT
end

function forward_substitution_sparse(L::Matrix{T}, b::Vector{T}, neighs::Array{Vector{Int},1}) where T<:AbstractFloat
    n = size(L, 1)
    y = zeros(T, n)
    for i in 1:n
        sum1 = 0.0
        for j in neighs[i]
            if j < i
                sum1 += L[i,j] * y[j]
            end
        end
        y[i] = b[i] - sum1
    end
    return y
end

function forward_substitution_transposed_sparse(L::Matrix{T}, b::Vector{T}, neighs::Array{Vector{Int},1}) where T<:AbstractFloat
    n = size(L, 1)
    x = zeros(T, n)
    for i in n:-1:1
        sum1 = 0.0
        for j in neighs[i]
            if j > i
                sum1 += L[j,i] * x[j]
            end
        end
        x[i] = b[i] - sum1
    end
    return x
end

function solve_LDLT_sparse(L::Matrix{T}, D::Vector{T},  neighs::Array{Vector{Int},1}, b::Vector{T} ) where T<:AbstractFloat
    z = forward_substitution_sparse(L, b, neighs)
    y = z ./ D
    x = forward_substitution_transposed_sparse(L, y, neighs)
    return x
end


const N_MAX_NEIGH = 32

function find_or_add_to_neighlist!(neighs::Matrix{Int}, i::Int, n::Int)
    for j in 1:N_MAX_NEIGH
        if neighs[j,i] == n
            return j  # n is already in the list, return its position
        elseif neighs[j,i] == 0
            neighs[j,i] = n  # n is not in the list, add it here
            return j  # return the position where n was added
        end
    end
    return 0  # list is full, couldn't add n
end

function CholeskyDecomp_LDLT_sparse_m(A::Matrix{Float64}, neighs::Matrix{Int}, tol::Float64=1e-16)
    n = size(A, 1)
    L = Matrix{Float64}(I, n, n)
    D = zeros(Float64, n)
    nop = 0
    nngmax = 0

    for j in 1:n
        sum1 = 0.0
        for k in 1:N_MAX_NEIGH
            if neighs[k,j] == 0
                break
            end
            if neighs[k,j] < j
                sum1 += (L[j,neighs[k,j]]^2) * D[neighs[k,j]]
                nop += 1
            end
        end
        D[j] = A[j,j] - sum1

        for i in j+1:n
            sum1 = 0.0
            for k in 1:N_MAX_NEIGH
                if neighs[k,j] == 0
                    break
                end
                if neighs[k,j] < j
                    sum1 += L[i,neighs[k,j]] * L[j,neighs[k,j]] * D[neighs[k,j]]
                    nop += 1
                end
            end
            Lij = (A[i,j] - sum1) / D[j]
            if abs(Lij) > tol
                L[i,j] = Lij
                find_or_add_to_neighlist!(neighs, j, i)
                find_or_add_to_neighlist!(neighs, i, j)
            end
        end
        nngmax = max(nngmax, count(!iszero, view(neighs,:,j)))
    end

    return L, D, neighs
end

function forward_substitution_sparse_m(L::Matrix{Float64}, b::Vector{Float64}, neighs::Matrix{Int})
    n = size(L, 1)
    y = zeros(Float64, n)
    for i in 1:n
        sum1 = 0.0
        for k in 1:N_MAX_NEIGH
            if neighs[k,i] == 0
                break
            end
            j = neighs[k,i]
            if j < i
                sum1 += L[i,j] * y[j]
            end
        end
        y[i] = b[i] - sum1
    end
    return y
end

function forward_substitution_transposed_sparse_m(L::Matrix{Float64}, b::Vector{Float64}, neighs::Matrix{Int})
    n = size(L, 1)
    x = zeros(Float64, n)
    for i in n:-1:1
        sum1 = 0.0
        for k in 1:N_MAX_NEIGH
            if neighs[k,i] == 0
                break
            end
            j = neighs[k,i]
            if j > i
                sum1 += L[j,i] * x[j]
            end
        end
        x[i] = b[i] - sum1
    end
    return x
end

function solve_LDLT_sparse_m(L::Matrix{Float64}, D::Vector{Float64}, b::Vector{Float64}, neighs::Matrix{Int})
    z = forward_substitution_sparse(L, b, neighs)
    y = z ./ D
    x = forward_substitution_transposed_sparse(L, y, neighs)
    return x
end