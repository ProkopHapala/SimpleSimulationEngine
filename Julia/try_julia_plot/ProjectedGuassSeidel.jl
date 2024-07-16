# Load the essential modules
using LinearAlgebra
#using SparseArrays
#using Plots

# =========== Functions

function buildRope( n::Int, m::Float64=1.0, m_end::Float64=1000.0, l::Float64=1.0 )
    bonds = [ (i,i+1) for i=1:(n-1) ]
    masses = [m_end; ones(n-2)*m; m_end]
    points = zeros(n, 3)
    x0     = -l*(n)/2
    for i=1:n
        points[i, 1] = (i-0.5)*l + x0
    end
    return bonds, points, masses
end

function buildGrid2D( nx::Int,ny::Int; m::Float64=1.0, m_end::Float64=1000.0, l::Float64=1.0, k::Float64=1.0, kDiag::Float64=-1.0 )
    np = (nx+1)*(ny+1)
    #print("np :"); display(np)
    masses = ones(np)*m
    masses[ 1    ] = m_end
    masses[ nx+1 ] = m_end
    points = zeros(np, 3)
    #print("masses :"); display(masses)
    # points loop over iy,ix
    for iy=1:ny+1
        for ix=1:nx+1
            i = (iy-1)*(nx+1) + ix
            masses[i] = m
            points[i, 1] = (ix-1.)*l
            points[i, 2] = (iy-1.)*-l
        end
    end
    #print("points :"); display(points)
    # bonds loop over iy,ix
    bonds = Array{Tuple{Int,Int},1}(undef, 0)
    ks    = Array{Float64,1}(undef, 0)
    #fixed = Set{Int}()
    fixed = [1, nx+1]
    for iy=1:ny+1
        for ix=1:nx+1
            i = (iy-1)*(nx+1) + ix
            # horizontal
            if ix <= nx
                j = i+1
                #println("i,j : ", i," ", j, "  iy,ix ", iy," ",ix  )
                push!( bonds, (i,j) )
                push!( ks   , k )
            end
            # vertical
            if iy <= ny
                j = i + nx + 1
                push!( bonds, (i,j) )
                push!( ks   , k )
            end
            if( kDiag > 0 ) # diagonal
                if (ix <= nx) && (iy <= ny)
                    j = i + nx + 2
                    push!( bonds, (i,j) )
                    push!( ks   , kDiag )
                end
                if (ix > 1) && (iy <= ny)
                    j = i + nx
                    push!( bonds, (i,j) )
                    push!( ks   , kDiag )
                end
            end 
        end
    end
    #print("bonds :"); display(bonds)

    println( "bonds:",  typeof(bonds) )
    println( "points:", typeof(points) )
    println( "masses:", typeof(masses) )
    println( "ks:", typeof(ks) )

    println("buildGrid2D() DONE !!!!!!!!!!!!!!!!!!!!! ")

    return bonds, points, masses, ks, fixed
end

function process_bonds( bonds::Array{Tuple{Int,Int},1}, points::Array{Float64,2} )
    #println("process_bonds");
    n  = length(bonds)
    ls = zeros(n)
    hs = zeros(n, 3)
    #print("points :"); display(points)
    for ib=1:n
        (i,j) = bonds[ib]
        d = points[i,:] - points[j,:]
        l = norm(d)
        ls[ib]  = l
        hs[ib,:] = d / l
    end
    #println("hs:"); display(hs)
    return hs, ls
end

function point_neighBs( bonds::Array{Tuple{Int,Int},1}, np::Int )
    #np = size(points, 1)
    #nb = size(bonds,  1)
    #neighBs = Vector{Vector{Int}}(undef, nps)
    neighBs = [ Vector{Int}() for _ in 1:np ]    
    for (ib, (i,j)) in enumerate(bonds)
        #println("ib, (i,j) : ", ib, " ", (i,j) )
        push!( neighBs[i], ib )
        push!( neighBs[j], ib )
    end
    return neighBs
end

function point_neighs(bonds::Array{Tuple{Int,Int},1}, np::Int, bSort::Bool=true)
    #np = size(points, 1)
    #nb = size(bonds,  1)
    neighs = [ Vector{Int}() for _ in 1:np ]    
    for (ib, (i,j)) in enumerate(bonds)
        #println("ib, (i,j) : ", ib, " ", (i,j) )
        push!( neighs[i], j )
        push!( neighs[j], i )
    end
    # sort vectors inside neighs
    if bSort
        for i in 1:np
            sort!(neighs[i])
        end
    end
    return neighs
end  

function point_neighs2( neighs::Array{Vector{Int},1} )
    # this function finds the neighbors of the neighbors ( i.e. second order neighbors )
    np = length(neighs)
    #ngmax = maximum( [length(ng) for ng in neighs] )
    #mask = zeros(Bool, ngmax)
    found = Set{Int}()
    neighs2 = [ Vector{Int}() for _ in 1:np ]
    for i in 1:np
        ngi=neighs[i]
        empty!(found)
        for j in ngi
            ngj=neighs[j]
            for k in ngj
                if( k != i )
                    #found.add(k)
                    push!(found, k)
                end
            end
        end
        neighs2[i] = collect(found) # convert Set to Vector
        sort!(neighs2[i])
    end
    return neighs2
end
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

function build_MassMatrix(ms::Array{Float64,1} )
    #nps    = length(ms)
    M_diag = repeat(ms, inner=3 )
    M      = Diagonal(M_diag)
    return M
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