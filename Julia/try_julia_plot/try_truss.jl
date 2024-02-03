include("ProjectedGuassSeidel.jl")
include("plot_utils.jl")

# =========== Data

"""
# ==== Sytem 1 : triangle (3 points, 3 sticks)

masses = [10.0, 1.0, 1.5]

points = [
    0.0  0.0  0.0;
    1.0 -1.0  0.0;
   -1.0 -1.0  0.0
]

bonds = [
    (1, 2),
    (2, 3),
    (3, 1)
]
"""

# ==== Sytem 2 : rope with 10 segments (11 points, 10 sticks)

#bonds, points, masses = buildRope( 11 )
nx = 3; ny = 3
bonds, points0, masses, ks = buildGrid2D( nx,ny, kDiag=1.0 )
_, l0s = process_bonds( bonds, points0 );

#print("bonds :");  display(bonds)
#print("points :"); display(points)
#print("masses :"); display(masses)

plt = plot(legend=false, aspect_ratio = :equal )

#plot_truss( plt, bonds, points0, c=:blue )

# add random displacement to the points
points = points0 + 0.1*randn(size(points0))
points[1   ,:] = points0[1   ,:]
points[nx+1,:] = points0[nx+1,:]


# =========== Main Body

neighBs = point_neighBs( bonds, points );
hs,ls   = process_bonds( bonds, points );

dl     = ls - l0s
strain = dl ./ l0s; #print("strain : "); display(strain)

#A_      = build_PGSMatrix_direct( bonds, hs, l0s, masses, neighBs );

# A, J, M = build_PGSMatrix(bonds, points, masses)

# Output the final matrix J and A to the console
#print("J:"); display(J)
#print("A :"); display(A )
#print("A_:"); display(A_)

plot_truss( plt, bonds, points, strain=strain*10.0, lw=3.0 )
display( plt )
