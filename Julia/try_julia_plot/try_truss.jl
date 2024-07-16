using Random

include("ProjectedGuassSeidel.jl")
include("plot_utils.jl")
include("IterativeSolvers.jl")
include("ProjectiveDynamics.jl")


"""
This is example running dynamical simulation of 2D truss (cloth) using projective dynamics.
"""



"""

ToDo:
 * Algorithm seems to work fine with  Julia's native Cholesky factorization (  AFact = cholesky(A) ) 
 * Test what is the problem with using our factorization algorithm
    * Is the problem with factorization, or with forward/back subsititution ?
      * seems problem is with substitution, because error is small (<1e-13), and using "\" operator insted of forwardsub() and  backsub() works even with L,U (from our Cholesky factorization ) as well as with  AFact (from Julia's cholesky(A) )

"""

# =========== Data

function test_Factorization( A::Array{Float64,2} )
plt = plot(legend=false, aspect_ratio = :equal )
#Ch = cholesky(A);             #;print("Cholesky.L : "); display(Ch.L)
#L  = CholeskyDecomp( A )       ;print("L : "); display(L)
#L_ = CholeskyDecomp_Crout( A ) ;print("L_ : "); display(L_)
#Ch = cholesky(A);              ;print("Cholesky.L : "); display(Ch.L)
print("CholeskyDecomp_Crout()");@time   L_C = CholeskyDecomp_Crout( A )          #;print("L_C: "); display(L_C)
# L_   = CholeskyDecomp_sparse( A, neighs ) #;print("L_: "); display(L_)
print("CholeskyDecomp_sparse()");@time L,neighsCh  = CholeskyDecomp_sparse( A, neighs, neighs2 ) #;print("L: "); display(L)
#print("neighsCh"); display(neighsCh)
#plot_matrix_log(plt, L-L_C )
#plot_matrix_log(plt, abs.(L) .> 1e-16 )
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
print("forwardsub()");        @time y_ch    = forwardsub(L,b)  #;print("y_ch : "); display(y_ch)
print("backsub()");           @time x_ch    = backsub(U,y_ch)  #;print("x_ch : "); display(x_ch)
print("mapMatrixNeighs()");   @time neighsL = mapMatrixNeighs(L)    #;print("neighsL"); display(neighsL)   
print("invNeighs()"      );   @time neighsU = invNeighs( neighsL )
print("forwardsub_sparse()"); @time y_ch_   = forwardsub_sparse(L,b   ,neighsL)  #;print("y_ch : "); display(y_ch)
print("backsub_sparse()");    @time x_ch_   = backsub_sparse(   U,y_ch,neighsU)  #;print("x_ch : "); display(x_ch)
println("max_error=max(abs(y_ch-y_ch_)) : ", maximum(abs.(y_ch-y_ch_)) );
println("max_error=max(abs(x_ch-x_ch_)) : ", maximum(abs.(x_ch-x_ch_)) );
display( plt ); return
x_ref = A \ b           #;print("x_ref : ");display(x_ref)
x_err = x_ch - x_ref    #;print("x_err : ");display(x_err)
end


function buildTruss_1()
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
return bonds, points, masses
end 


###################
###    main()   ###
###################

function main( nx=5, ny=5, dt=5.0 )
#function main( nx=5, ny=5, dt=5.0 )    
# ==== Sytem 1 : triangle (3 points, 3 sticks)
# ==== Sytem 2 : rope with 10 segments (11 points, 10 sticks)
#bonds, points, masses = buildRope( 11 )
#nx = 3; ny = 3
#nx = 3; ny = 3
#nx = 5; ny = 5
#nx = 10; ny = 10
#nx = 30; ny = 30
#nx = 50; ny = 50
#dt = 5.0
gravity = [0.0, -9.81, 0.0 ]  # gravity

# ===========  Generate System Gometry
# ---- Generate 2D truss, points, bonds, masses, stiffness, relaxed lenghts
print("buildGrid2D()");   @time   bonds, points0, masses, ks,fixed = buildGrid2D( nx,ny, k=1000.0, kDiag=10.0 )
print("process_bonds()"); @time   _, l0s = process_bonds( bonds, points0 );

#println( size(points), size(gravity) )
#f_grav = [ gravity .* m for m in masses ]   
#println( "typeof(masses) ", typeof(masses) )

f_grav = masses .* gravity'    #;print("f_grav : "); display(f_grav)

plt = plot(legend=false, aspect_ratio = :equal )
#plot_truss( plt, bonds, points0, c=:blue )

# add random displacement to the points
Random.seed!(1234)  # Sets the seed to 1234
points = points0 .+ 0.1*randn(size(points0))
points[fixed,:] .= points0[fixed,:]
masses[fixed  ] .= 1000000.0
#points[nx+1,:] = points0[nx+1,:]

# ===========  Generate system matrix ( Projective Dynamics Matrix ) using neighborhood information

np = size(points,1)
print("point_neighs()");   @time neighs  = point_neighs(  bonds, np )     #;print("neighs : "); display(neighs)      # finds neighbors
print("point_neighBs()");  @time neighBs = point_neighBs( bonds, np )     #;                                         # finds bonds to neighbors
print("point_neighs2()");  @time neighs2 = point_neighs2( neighs    )     #; print("neighs2 : "); display(neighs2)   # finds neighbors of neighbors ( i.e. second order neighbors )

#print("bonds :");  display(bonds)
#print("points :"); display(points)
#print("masses :"); display(masses)
print("make_PDMatrix()");  @time  A = make_PD_Matrix( neighBs, bonds, masses, dt, ks )     # ;print( "A : "); display(A)         # creates projective dynamics matrix
#make_PD_Matrix( neighBs::Array{Vector{Float64},1}, bonds::Array{Tuple{Int,Int},1}, masses::Array{Float64,1}, dt::Float64, ks::Array{Float64,1} )
#print("make_PDMatrix()"); @time  A = make_PD_Matrix( neighBs, bonds, masses, dt, ks )                     # ;print( "A : "); display(A)
#print("make_PD_rhs()");   @time  b = make_PD_rhs(   neighBs, bonds, masses, dt, ks, points, l0s, pnew )  # ;print( "b : "); display(b)
#print(" ==== projective dynamics matrix A: "); display(A)

# =========== Matrix Factorization (or inversion) for efficient solution of linear system later

print("CholeskyDecomp_sparse()"); @time L,neighsCh  = CholeskyDecomp_sparse( A, neighs, neighs2 ) #;print("L: "); display(L)    # evaluate Cholensky decomposition for sparse matrix (efficient)
print("CholeskyDecomp_Crout() "); @time L_C         = CholeskyDecomp_Crout( A ) #;print("L: "); display(L)                      # evaluate Cholensky decomposition using Cholesky–Crout algorithm
U  = copy(L') 

"""
# ----- Check our Cholensky decomposition versus Julia's native cholesky() algorithm
println("CholeskyDecomp_sparse() max_error=max(abs(L*L^T-A)) : ", maximum(abs.(L*U      - A)) );                     #  evaluate error of CholeskyDecomp_sparse()
println("CholeskyDecomp_Crout()  max_error=max(abs(L*L^T-A)) : ", maximum(abs.(L_C*L_C' - A)) );                     #  evaluate error of CholeskyDecomp_Crout() 
#AFact = factorize(A) ;print("AFact : "); display(AFact)
AFact = cholesky(A) #;print("AFact : "); display(AFact)                                                              #  evaluate Cholensky decomposition using Julia's native cholesky() algorithm
dL = L - AFact.L     # evaluate error of CholeskyDecomp_sparse()  with respect to Julia's native cholesky() algorithm
dU = U - AFact.U     # evaluate error of CholeskyDecomp_sparse()  with respect to Julia's native cholesky() algorithm
plot_matrix_log(plt, dL )  #  plot error of  CholeskyDecomp_sparse()  with respect to Julia's native cholesky() algorithm
plot_matrix_log(plt, dU )  #  plot error of  CholeskyDecomp_sparse()  with respect to Julia's native cholesky() algorithm                                    
"""

print("mapMatrixNeighs()"); @time neighsL = mapMatrixNeighs(L)    #;print("neighsL"); display(neighsL)   # find which elements of L are non-zero ( for each row    i map all non-zero Aij to columns j )
print("invNeighs()");       @time neighsU = invNeighs( neighsL )      #;print("neighsU"); display(neighsU)     # find non-zero back-neighbors          ( for each column j map all non-zero Aij to rows    i )

#plot_matrix_log(plt, L )
#plot_matrix_log(plt, abs.(L) .> 1e-16 )

print("inv(A)");          @time invA = inv(A)      # direct marix inversion (for speed comparison)
#;print("invA : ");    display(invA)
#plot_matrix_log(plt, invA )
#plot_matrix_log(plt, abs.(invA) .> 1e-6 )

# =========== Matrix Factorization (or inversion) for efficient solution of linear system later

print("process_bonds()"); @time hs,ls = process_bonds( bonds, points )  #;

dl     = ls - l0s
strain = dl ./ l0s; #print("strain : "); display(strain)

#Ad      = build_PGSMatrix_direct( bonds, hs, l0s, masses, neighBs );

# linear solve for the displacement x = A_ \cdot b
#x = Ad \ dl

#print("x : "); display(x)
println( size(points), size(gravity) )
f_ext  = f_grav .* 1.0
println( " size(f_ext) 1 ",  size(f_ext) )
#f_ext   = f_ext .+ gravity #.* masses   
#println( " size(f_ext) 2 ",  size(f_ext) )

v       = points .*0.0 
pnew    = points .+ v*dt .+ f_ext .* (dt^2)*0.0

# ========== Solver 

plot_truss( plt, bonds, points0, lw=1.0, c=:black )
#plot_truss( plt, bonds, points,  lw=1.0, strain=strain*10.0, )
plot_truss( plt, bonds, points,  lw=1.0, c=:red )

#bSparseCholesky = false
bSparseCholesky = true
niter = 100
for i=1:niter
    pnew   = points .+ v*dt .+ f_ext*(dt^2)
    pnew[fixed,:] .= points0[fixed,:]
    
    b      = make_PD_rhs( neighBs, bonds, masses, dt, ks, points, l0s, pnew )  # ;print( "b : "); display(b)
            
        # ---- Method 1)
        #y    = forwardsub(L,b   )  #;print("y_ch : "); display(y_ch)
        #pnew = backsub(   U,y_ch)  #;print("x_ch : "); display(x_ch)

        # ---- Method 2)
        y    = forwardsub_sparse(L,b,neighsL)  #;print("y_ch : "); display(y_ch)
        pnew = backsub_sparse(   U,y,neighsU)  #;print("x_ch : "); display(x_ch)

        # ---- Method 3)   using  forward-and-backsubstitution with Our Cholensky factorization
        #y    = forwardsub(L,b)
        #pnew = backsub(   U,y)

        # ---- Method 4)   using  forward-and-backsubstitution with Our Cholensky factorization from Julia's native algorithm
        #L = Matrix(AFact.L)
        #U = Matrix(AFact.U)
        #y    = forwardsub(L,b)
        #pnew = backsub(   U,y)

        # ---- Method 5b)   using factorization by Julia's native algorithm
        #y    = L \ b
        #y    = forwardsub(L,b)
        #pnew = U \ y
        #pnew    = forwardsub(U,y)
        #pnew = backsub(   U,y)

        # ---- Method 5)   using factorization by Julia's native algorithm
        #y    = AFact.L \ b
        #pnew = AFact.U \ y

        # ---- Method 6)   multiply by inverse matrix
        #pnew = AFact \ b

        # ---- Method 8)   multiply by inverse matrix
        #pnew = invA * b 

        # ---- Method 9)   direct sover
        #pnew = A \ b

    dpos   =  pnew - points
    #print("dpos[$i] " ); display(dpos)
    res    = maximum( abs.(dpos) )   ;println("residual[$i] : ", res );
    v      = dpos/dt
    
    points[:,:] .= pnew[:,:]
    #points = points .+ x
    #plot!( plt, [points0[i,1],points[i,1]], [points0[i,2],points[i,2]], color=:blue, lw=0.5 )

    #plot_truss( plt, bonds, points, lw=1.0, c=:blue )
    #plot_truss( plt, bonds, points0, lw=1.0, c=:blue )
end

plot_truss( plt, bonds, points, lw=1.0, c=:blue )

display( plt ); return

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
==#
display( plt )
end # function mai


main()