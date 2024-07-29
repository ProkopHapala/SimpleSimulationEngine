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
    #println( "bonds:",  typeof(bonds) )
    #println( "points:", typeof(points) )
    #println( "masses:", typeof(masses) )
    #println( "ks:", typeof(ks) )
    #println("buildGrid2D() DONE !!!!!!!!!!!!!!!!!!!!! ")
    return bonds, points, masses, ks, fixed
end


struct TrussEdge
    i::Int64
    j::Int64
    kind::Int64
end

# mutable struct Truss
#    points::Vector{Vector{Float64}}
#    edges::Vector{TrussEdge}
#    blocks::Vector{Tuple{Int64, Int64}}
# end

function make_ortho_u!(ax::Vector{Float64}, dir::Vector{Float64})
    ax -= (ax ⋅ dir) * dir
    return ax / norm(ax) 
end

# function from_angle!(angle::Float64)
#     return [cos(angle), sin(angle)]
# end

# function mul_cmplx!(v::Vector{Float64}, drot::Vector{Float64})
#     a = v[1] * drot[1] - v[2] * drot[2]
#     b = v[1] * drot[2] + v[2] * drot[1]
#     return [a, b]
# end


function ngonTruss(p0::Vector{Float64}, p1::Vector{Float64}, ax::Vector{Float64}; n::Int=8, k::Float64=1.00 )
    #println("Truss::ngonTruss() n=$n p0=($p0[1],$p0[2],$p0[3]) p1=($p1[1],$p1[2],$p1[3]) ax=($ax[1],$ax[2],$ax[3])")

    dir  = p1 - p0
    r    = norm(dir)
    dir  = dir / norm(dir)
    ax   = make_ortho_u!(ax, dir)
    side = cross(dir, ax)

    rot  = 1.0 + 0.0im
    drot = cis(2*π / n)  # Using complex number representation for rotation

    points = Vector{Vector{Float64}}()
    edges  = Vector{TrussEdge}()

    io = n
    for i in 1:n
        R = dir * real(rot) + side * imag(rot)
        push!(points, p0 + R )
        push!(edges, TrussEdge(i, io, 1))
        rot *= drot  # Complex multiplication
        io=i
    end

    points = Matrix(hcat(points...)')

    bonds = Array{Tuple{Int, Int}, 1}(undef, 0)
    ks = Array{Float64, 1}(undef, 0)
    for e in edges
        push!(bonds, (e.i, e.j))
        push!(ks, k )
    end
    return points, bonds, ks
end



function wheel(p0::Vector{Float64}, p1::Vector{Float64}, ax::Vector{Float64}, width::Float64; n::Int=8,  kscale::Float64=1.0, kdict::Vector{Float64} = [1.0, 1.0, 1.0, 1.0] )
    #println("Truss::wheel() n=$n p0=($p0[1],$p0[2],$p0[3]) p1=($p1[1],$p1[2],$p1[3]) ax=($ax[1],$ax[2],$ax[3])")
    kind_long   = 1
    kind_perp   = 2
    kind_zigIn  = 3
    kind_zigOut = 4

    dir  = p1 - p0
    r    = norm(dir)
    dir  = dir / norm(dir)
    ax   = make_ortho_u!(ax, dir)
    side = cross(dir, ax)

    dnp  = 4
    i00  = 1
    i000 = i00

    rot  = 1.0 + 0.0im
    drot = cis(π / n)  # Using complex number representation for rotation

    points = Vector{Vector{Float64}}()
    edges  = Vector{TrussEdge}()

    for i in 1:n
        i01 = i00 + 1
        i10 = i00 + 2
        i11 = i00 + 3

        R = dir * real(rot) + side * imag(rot)
        push!(points, p0 + R * (r + width))
        push!(points, p0 + R * (r - width))
        rot *= drot  # Complex multiplication
        R = dir * real(rot) + side * imag(rot)
        push!(points, p0 + ax *  width + R * r)
        push!(points, p0 + ax * -width + R * r)
        rot *= drot  # Complex multiplication

        push!(edges, TrussEdge(i00, i01, kind_perp))
        push!(edges, TrussEdge(i10, i11, kind_perp))
        push!(edges, TrussEdge(i00, i10, kind_zigIn))
        push!(edges, TrussEdge(i00, i11, kind_zigIn))
        push!(edges, TrussEdge(i01, i10, kind_zigIn))
        push!(edges, TrussEdge(i01, i11, kind_zigIn))
        if i < n
            push!(edges, TrussEdge(i10, i00 + dnp, kind_zigOut))
            push!(edges, TrussEdge(i10, i01 + dnp, kind_zigOut))
            push!(edges, TrussEdge(i11, i00 + dnp, kind_zigOut))
            push!(edges, TrussEdge(i11, i01 + dnp, kind_zigOut))
            push!(edges, TrussEdge(i00, i00 + dnp, kind_long))
            push!(edges, TrussEdge(i01, i01 + dnp, kind_long))
            push!(edges, TrussEdge(i10, i10 + dnp, kind_long))
            push!(edges, TrussEdge(i11, i11 + dnp, kind_long))
        else
            push!(edges, TrussEdge(i10, i000 + 0, kind_zigOut))
            push!(edges, TrussEdge(i10, i000 + 1, kind_zigOut))
            push!(edges, TrussEdge(i11, i000 + 0, kind_zigOut))
            push!(edges, TrussEdge(i11, i000 + 1, kind_zigOut))
            push!(edges, TrussEdge(i00, i000 + 0, kind_long))
            push!(edges, TrussEdge(i01, i000 + 1, kind_long))
            push!(edges, TrussEdge(i10, i000 + 2, kind_long))
            push!(edges, TrussEdge(i11, i000 + 3, kind_long))
        end
        i00 += dnp
    end

    points = Matrix(hcat(points...)')

    bonds = Array{Tuple{Int, Int}, 1}(undef, 0)
    ks = Array{Float64, 1}(undef, 0)
    for e in edges
        push!(bonds, (e.i, e.j))
        push!(ks, kdict[e.kind]*kscale )
    end
    return points, bonds, ks
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

function build_MassMatrix(ms::Array{Float64,1} )
    #nps    = length(ms)
    M_diag = repeat(ms, inner=3 )
    M      = Diagonal(M_diag)
    return M
end