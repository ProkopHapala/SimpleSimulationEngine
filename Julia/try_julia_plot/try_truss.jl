using Random

include("ProjectedGuassSeidel.jl")
include("plot_utils.jl")
include("IterativeSolvers.jl")
include("ProjectiveDynamics.jl")

# =========== Data

function main()

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


# ==== Sytem 2 : rope with 10 segments (11 points, 10 sticks)

#bonds, points, masses = buildRope( 11 )
#nx = 3; ny = 3
#nx = 3; ny = 3

#nx = 10; ny = 10
nx = 30; ny = 30
#nx = 50; ny = 50

print("buildGrid2D()"); @time   bonds, points0, masses, ks = buildGrid2D( nx,ny, kDiag=1.0 )
print("process_bonds()"); @time   _, l0s = process_bonds( bonds, points0 );


#print("bonds :");  display(bonds)
#print("points :"); display(points)
#print("masses :"); display(masses)

plt = plot(legend=false, aspect_ratio = :equal )

#plot_truss( plt, bonds, points0, c=:blue )

# add random displacement to the points
Random.seed!(1234)  # Sets the seed to 1234
points = points0 + 0.1*randn(size(points0))
points[1   ,:] = points0[1   ,:]
points[nx+1,:] = points0[nx+1,:]

# =========== Main Body
np = size(points,1)
print("point_neighs()"); @time neighs  = point_neighs( bonds,  np )     #;print("neighs : "); display(neighs)
print("point_neighBs()"); @time neighBs = point_neighBs( bonds, np )      #;
print("process_bonds()"); @time hs,ls   = process_bonds( bonds, points )  #;

#neighsets::Array{Set{Int}} = [ Set(neighs[i]) for i in 1:np ] 

print("point_neighs2()");@time neighs2 = point_neighs2( neighs )    #; print("neighs2 : "); display(neighs2)

#return

dl     = ls - l0s
strain = dl ./ l0s; #print("strain : "); display(strain)

#Ad      = build_PGSMatrix_direct( bonds, hs, l0s, masses, neighBs );

# linear solve for the displacement x = A_ \cdot b
#x = Ad \ dl
dt = 0.1
#print("x : "); display(x)
gravity = [0.0, 0.0, -9.81] 
pnew    = pnew = points .+ gravity' .* (dt^2)

dt = 0.1
print("make_PDMatrix()");@time  A = make_PDMatrix( neighBs, bonds, masses, dt, ks )                     # ;print( "A : "); display(A)
print("make_PD_rhs()");@time    b = make_PD_rhs(   neighBs, bonds, masses, dt, ks, points, l0s, pnew )  # ;print( "b : "); display(b)

#Ch = cholesky(A);             #;print("Cholesky.L : "); display(Ch.L)
#L  = CholeskyDecomp( A )       ;print("L : "); display(L)
#L_ = CholeskyDecomp_Crout( A ) ;print("L_ : "); display(L_)
#Ch = cholesky(A);              ;print("Cholesky.L : "); display(Ch.L)
print("CholeskyDecomp_Crout()");@time   L_C = CholeskyDecomp_Crout( A )          #;print("L_C: "); display(L_C)


# L_   = CholeskyDecomp_sparse( A, neighs ) #;print("L_: "); display(L_)


print("CholeskyDecomp_sparse()");@time L,neighsCh  = CholeskyDecomp_sparse( A, neighs, neighs2 ) #;print("L: "); display(L)


#print("neighsCh"); display(neighsCh)

#plot_matrix_log(plt, L-L_C )
plot_matrix_log(plt, abs.(L) .> 1e-16 )
#plot_matrix_log(plt, L-L_ )



U  = copy(L') 
println("CholeskyDecomp_sparse() max_error=max(abs(L*L^T-A)) : ", maximum(abs.(L*U      - A)) );
println("CholeskyDecomp_Crout()  max_error=max(abs(L*L^T-A)) : ", maximum(abs.(L_C*L_C' - A)) );

#plot_matrix_log(plt, Ch.L )
#plot_matrix_log(plt, L )
#plot_matrix_bwr( plt, err_Ch )
#plot_matrix_log(plt, Ch.U )

#y_ch = forwardsub_bak(L,b[:,1])  #;print("y_chb : "); display(y_ch)
#x_ch = backsub_bak(  U,y_ch)     ;print("x_chb : "); display(x_ch)


print("forwardsub()");@time y_ch = forwardsub(L,b)  #;print("y_ch : "); display(y_ch)
print("backsub()");@time    x_ch = backsub(U,y_ch)  #;print("x_ch : "); display(x_ch)

neighsL = mapMatrixNeighs(L)    #;print("neighsL"); display(neighsL)   
neighsU = invNeighs( neighsL )

print("forwardsub()");@time y_ch_ = forwardsub_sparse(L,b   ,neighsL)  #;print("y_ch : "); display(y_ch)
print("backsub()");@time    x_ch_ = backsub_sparse(   U,y_ch,neighsU)  #;print("x_ch : "); display(x_ch)

println("max_error=max(abs(y_ch-y_ch_)) : ", maximum(abs.(y_ch-y_ch_)) );
println("max_error=max(abs(x_ch-x_ch_)) : ", maximum(abs.(x_ch-x_ch_)) );



display( plt ); return

x_ref = A \ b           #;print("x_ref : ");display(x_ref)

x_err = x_ch - x_ref    #;print("x_err : ");display(x_err)




#==


x0 = zeros(size(b)) 
x  = SolveIterative( A, b, x0, niter=10, tol=1.e-3, method=1 )   #; print("x : "); display(x)

# solve A*x=b using direct solver
x_ref = A \ b     #;print("x_ref : "); display(x_ref)


# invert matrix A_
invA = inv(A)      #;print("invA : "); display(invA)
x_inv = invA * b   #;print("x_inv : "); display(x_inv)

# if element of invA<trahold, set it to zero
threshold = 1.e-6
invA[ abs.(invA) .< threshold ] .= 0.0

# cholensky factorization of A

#plot_matrix_bwr(plt, A )
#plot_matrix_log(plt, A )
#plot_matrix_log(plt, invA )
plot_matrix_log(plt, L )

# A, J, M = build_PGSMatrix(bonds, points, masses)

# Output the final matrix J and A to the console
#print("J:"); display(J)
#print("A :"); display(A )
#print("A_:"); display(A_)

plot_truss( plt, bonds, points, strain=strain*10.0, lw=1.0 )
display( plt )

==#
display( plt )
end # function mai


main()