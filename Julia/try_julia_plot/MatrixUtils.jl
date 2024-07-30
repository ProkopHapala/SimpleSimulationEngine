# Load the essential modules
#using LinearAlgebra
#using SparseArrays
#using Plots

using DataStructures

# =========== Functions

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

function forward_substitution(L::Matrix{T}, b::Vector{T}; bPrint::Bool=:false) where T<:AbstractFloat
    n = size(L, 1)
    x = zeros(T, n)
    for i in 1:n
        sum  = dot(L[i,1:i-1], x[1:i-1])
        x[i] = b[i] - sum
        if(bPrint) 
            #println("fwsub()[",i,"](x=",x[i],",b=",b[i],",sum=",sum,")")
            println("fwsub()sum[",i,"]=   ",sum )
        end
    end
    return x
end

function forward_substitution_transposed(L::Matrix{T}, b::Vector{T}) where T<:AbstractFloat
    n = size(L, 1)
    x = zeros(T, n)
    for i in n:-1:1
        x[i] = b[i] - dot(L[i+1:n,i], x[i+1:n])
    end
    return x
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

# function forward_substitution_sparse_m(L::Matrix{Float64}, b::Vector{Float64}, neighs::Matrix{Int})
#     n = size(L, 1)
#     y = zeros(Float64, n)
#     for i in 1:n
#         sum1 = 0.0
#         for k in 1:N_MAX_NEIGH
#             if neighs[k,i] == 0
#                 break
#             end
#             j = neighs[k,i]
#             if j < i
#                 sum1 += L[i,j] * y[j]
#             end
#         end
#         y[i] = b[i] - sum1
#     end
#     return y
# end

# function forward_substitution_transposed_sparse_m(L::Matrix{Float64}, b::Vector{Float64}, neighs::Matrix{Int})
#     n = size(L, 1)
#     x = zeros(Float64, n)
#     for i in n:-1:1
#         sum1 = 0.0
#         for k in 1:N_MAX_NEIGH
#             if neighs[k,i] == 0
#                 break
#             end
#             j = neighs[k,i]
#             if j > i
#                 sum1 += L[j,i] * x[j]
#             end
#         end
#         x[i] = b[i] - sum1
#     end
#     return x
# end

