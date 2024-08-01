# Load the essential modules
#using LinearAlgebra
#using SparseArrays
#using Plots

using LinearAlgebra
using Random

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

function SOR!(A::AbstractMatrix, res::AbstractVector, z::AbstractVector, omega::Float64)
    n = length(res)
    for i in 1:n
        sum = 0.0
        for j in 1:n
            if i != j
                sum += A[i, j] * z[j]
            end
        end
        z[i] = (1 - omega) * z[i] + omega * (res[i] - sum) / A[i, i]
    end
end

function conjugate_gradient!(  A::AbstractMatrix,  b::AbstractVector,   x::AbstractVector;  tol=1e-6, niter::Int=10, bPrint::Bool=:false )

    # Initialize residual vector
    res = b - A * x
    dir = copy(res)     # Initialize search direction vector
    err = norm(res)     # Compute initial squared residual norm

    for iter in 1:niter
        A_dir = A * dir
        alpha = err^2 / dot(dir, A_dir )
        
        x   .+= dir   .* alpha
        res .-= A_dir .* alpha

        err_new = norm(res)
        if bPrint 
            println("CG[",iter,"] err=", err_new )
        end
        if err_new < tol
            return x, iter
        end
        beta = (err_new / err)^2
        dir .= res .+ ( dir .* beta  ) 
        err = err_new
    end
    return x, -1
end

function conjugate_gradient_jac!(  A::AbstractMatrix,  b::AbstractVector,   x::AbstractVector;  tol=1e-6, niter::Int=10, bPrint::Bool=:false )

    jacobi_precond = ( 1 ./ diag(A) ) #.+ 1e-16
    #M =  Diagonal(1 ./ diag(A))
    # Initialize residual vector
    res = b - A * x
    #z .*= jacobi_precond
    z = res .* jacobi_precond
    #z = M * res
    dir = copy(z)     # Initialize search direction vector
    err = dot(z,res)  # Compute initial squared residual norm

    for iter in 1:niter
        A_dir = A * dir
        alpha = err / dot(dir, A_dir )
        
        x   .+= dir   .* alpha
        res .-= A_dir .* alpha

        z = res .* jacobi_precond
        #z = M * res

        err_new = dot(res,z)     # reduction
        
        if bPrint 
            println("CG[",iter,"] err=", err_new )
        end
        if sqrt(err_new) < tol
            return x, iter
        end
        beta = err_new / err
        dir .= z .+ ( dir .* beta  ) 
        err = err_new
    end
    return x, niter
end


function conjugate_gradient_sor!(A::AbstractMatrix, b::AbstractVector, x::AbstractVector; tol=1e-6, niter::Int=10, omega::Float64=1.5, bPrint::Bool=false)
    n = length(b)
    res = b - A * x
    z = zeros(eltype(x), n)
    dir = zeros(eltype(x), n)
    err = 0.0

    jacobi_precond = ( 1 ./ diag(A) ) #.+ 1e-16

    # Initial SOR preconditioning
    SOR!(A, res, z, omega)
    #SOR!(A, res .* jacobi_precond , z, omega)

    dir .= z
    err = dot(z, res)

    for iter in 1:niter
        A_dir = A * dir
        alpha = err / dot(dir, A_dir)
        
        x   .+= dir   .* alpha
        res .-= A_dir .* alpha

        # SOR preconditioning step
        SOR!(A, res  , z, omega)
        #SOR!(A, res .* jacobi_precond  , z, omega)

        err_new = dot(res, z)
        
        if bPrint
            println("CG[", iter, "] err=", err_new)
        end
        if sqrt(err_new) < tol
            return x, iter
        end
        beta = err_new / err
        dir .= z .+ (dir .* beta)
        err = err_new
    end
    return x, niter
end


function conjugate_gradient_precond!(A::AbstractMatrix, b::AbstractVector, x::AbstractVector, M::AbstractMatrix; tol=1e-6, niter::Int=10, bPrint::Bool=:false )
    r = b - A * x
    z = M * r        # Preconditioning step as matrix multiplication
    p   = copy(z)
    err = dot(r, z)

    for iter in 1:niter
        Ap    = A * p
        alpha = err / dot(p, Ap)

        x .+= alpha .* p
        r .-= alpha .* Ap

        z = M * r        # Preconditioning step as matrix multiplication
        err_new = dot(r, z)
        if bPrint 
            println("CG[",iter,"] err=", err_new )
        end
        if sqrt(err_new) < tol
            return x, iter
        end

        beta = err_new / err
        p .= z .+ beta .* p

        err = err_new
    end

    return x, niter
end

# Jacobi preconditioner as a matrix
function jacobi_preconditioner_matrix(A::AbstractMatrix)
    return Diagonal(1 ./ diag(A))
end

# Symmetric Gauss-Seidel preconditioner approximation
function sgs_preconditioner_matrix(A::AbstractMatrix)
    D = Diagonal(A)
    L = LowerTriangular(A) - D
    return (D + L) * inv(D) * (D + L')
end


function generate_test_matrix( n::Int64, off::Float64; diag::Float64=1.0)
    A = randn(n, n)
    return (A + A')*off + I*diag
end

# Helper function to check if a matrix is SPD
function is_spd(A::AbstractMatrix)
    isapprox(A, A') || return false  # Check symmetry
    try
        cholesky(A)  # Check positive definiteness
        return true
    catch
        return false
    end
end

function test_CG(  n::Int64=3; bPre::Bool=:false, off::Float64=0.1, nIterMax::Int64=100, tol::Float64=1e-6 )
    A = generate_test_matrix(n,off)
    b = randn(n)

    println("A: "); display(A)
    
    # Solve using Julia's direct solver
    x_direct = A \ b
    x_cg = zeros(n)
    if bPre
        M_pre = jacobi_preconditioner_matrix(A)
        #M_pre = sgs_preconditioner_matrix(A)
        #conjugate_gradient_precond!(A, b, x_cg, M_pre; tol=tol, niter=nIterMax, bPrint=true)
        conjugate_gradient_sor!(A, b, x_cg, tol=tol, niter=nIterMax, bPrint=true)
    else
        conjugate_gradient_jac!(A, b, x_cg, tol=tol, niter=nIterMax, bPrint=true)
        #conjugate_gradient!(A, b, x_cg, tol=tol, niter=nIterMax, bPrint=true)
    end
    # Compare results
    rel_error = norm(x_cg - x_direct) / norm(x_direct)
    
    println("Matrix size:    n=$n  off=$off")
    println("Relative error:            $rel_error")
    println("Direct solver residual: ", norm(A * x_direct - b) ) 
    println("CG     solver residual: ", norm(A * x_cg - b    ) )
end