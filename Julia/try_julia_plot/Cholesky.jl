# Load the essential modules
#using LinearAlgebra
#using SparseArrays
#using Plots

using DataStructures

include("MatrixUtils.jl")

# =========== Functions

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

function solve_LDLT( L::Matrix{T}, D::Vector{T}, b::Vector{T}) where T<:AbstractFloat
    #L, D = CholeskyLDLDecomp(A)
    z = forward_substitution(L, b)
    #y = diagonal_solve(D, z)
    y = z ./ D
    x = forward_substitution_transposed(L, y)
    return x
end

# function inset_to_set(i::Int64,ng::Int64, st::Set{Int64})
#     if !( ng in st )
#         push!(st,ng)
#         println("insert neigh[",i,"].insert(",ng,")" )
#     end
# end

function CholeskyDecomp_LDLT_sparse(A::Matrix{T}, neighs::Array{Vector{Int},1}, tol::T=1e-16) where T<:AbstractFloat
    neighsets = [Set(ng) for ng in neighs]
    n = size(A, 1)
    L = Matrix{T}(I, n, n)
    D = zeros(T, n)
    nop = 0
    nngmax = 0

    for j in 1:n
        #println(" --- CholeskyDecomp_LDLT_sparse[",j,"]" );
        sum1 = 0.0
        ngs = neighsets[j]
        for k in ngs
            if k < j
                val = (L[j,k]^2) * D[k]
                sum1 += val
                #println("sum[$j,$k] $(lpad(val, 20, ' '))   $(lpad(sum1, 20, ' '))  ")
                nop += 1
            end
        end
        D[j] = A[j,j] - sum1
        #println(  "Ch[%i] D,A,sum  %20.10f   %20.10f   %20.10f ",  D[j],  A[j*n+j],  sum  );
        #println("Ch[$j] D,A,sum  $(lpad(D[j], 20, ' '))   $(lpad(A[j,j], 20, ' '))   $(lpad(sum1, 20, ' '))")
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

                #inset_to_set(i,j,neighsets[i])
                #inset_to_set(i,j,ngs)
                #push!(ngs,          i)

                push!(ngs,          i)
                push!(neighsets[i], j)
            end
        end
        nngmax = max(nngmax, length(ngs))
    end

    neighsLDLT = [sort!(collect(ngs)) for ngs in neighsets]
    return L, D, neighsLDLT
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