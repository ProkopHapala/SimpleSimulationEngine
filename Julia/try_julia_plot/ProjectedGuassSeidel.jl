# Load the essential modules
#using LinearAlgebra
#using SparseArrays
#using Plots

# ======= Bulding system as dense matrices ( A, J, M )

function build_Jacobian(bonds::Array{Tuple{Int,Int},1}, points::Array{Float64,2} )
    hs, l0s = process_bonds( bonds, points )
    nb      = length(bonds)
    np      = size(points, 1)
    J = zeros(nb, 3*np)
    #print("Jacobian size ", size(J) )
    for (ib, (i, j)) in enumerate(bonds)
        h = hs[ib,:]
        J[ib, (3*i-2):(3*i)] = -h
        J[ib, (3*j-2):(3*j)] =  h
    end
    #print("J:"); display(J)
    return J, l0s, hs
end

function build_PGSMatrix(bonds, points, masses)
    J, l0s, hs = build_Jacobian(bonds, points)
    M          = build_MassMatrix(masses)
    A          = J * M * J'
    return A, J, M
end

# ======= Bulding system as sparse solver

function build_PGSMatrix_direct( bonds::Array{Tuple{Int,Int},1}, hs::Array{Float64,2}, ms::Array{Float64,1}, neighBs::Array{Vector{Int},1} )
    """
    Calculate the dot product of A = J*M*J^T with a vector x
    k,l index of bonds
    i,j index of points 
    Akk =  (mi + mj) * dot( hi, hi ) = mi + mj
    Akl =  -mi       * dot( hk, hl )   ... only for k,l which share a point i
    """
    nb = size(hs, 1)
    A  = zeros(nb, nb)
    for k in 1:nst
        (i,j) = bonds[k]
        mi = ms[i]
        mj = ms[j]
        A[k,k] = mi + mj
        hk     = hs[k,:]
        for l in neighBs[i]
            if k < l
                akl = -mi * dot( hk, hs[l,:] )
                A[k,l] = akl
                A[l,k] = akl
            end
        end
        for l in neighBs[j]
            if k < l
                akl = -mj * dot( hk, hs[l,:] )
                A[k, l] = akl
                A[l, k] = akl
            end
        end
    end
    return A
end