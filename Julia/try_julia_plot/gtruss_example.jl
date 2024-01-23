# Load the essential modules
using LinearAlgebra
using SparseArrays
using Plots


# How to run:
# from bash: 
#  > julia gtruss_example.jl
# from vscode REPL:
# cd("./Julia/try_julia_plot")
# > include("gtruss_example.jl")

#
# Given connectivities and coordinates, evaluate both the
# length, the  angle between the local and
# the global reference systems for a given element ele and the global dofs
#
function Ele_Data( ele::Int64, coord::Array{Float64,2}, connect::Array{Int64,2} )

    # recover the nodes
    node1, node2 = connect[ele,:]

    # recover the coordinates
    x1,y1 = coord[node1,:]
    x2,y2 = coord[node2,:]

    # Evaluates L for this element
    L = sqrt((x2-x1)^2 + (y2-y1)^2)

    # Evaluates theta for this element. This usage of atan
    # is equivalent to the usual atan2 command used in matlab.
    theta = atan(y2-y1,x2-x1)

    # Find the global DOFs for this element
    dofs = [2*(node1-1)+1 ; 2*(node1-1)+2 ; 2*(node2-1)+1 ; 2*(node2-1)+2  ]

    # Return L, theta and dofs for this element
    return L, theta, dofs

end

#
# Returns the local stiffness matrix of an element
#
function Local_Stiffness( E::Float64, A::Float64, L::Float64 )
         (E*A/L)*[1.0 -1.0 ; -1.0 1.0 ]
end

#
# Returns the rotation matrix of an element
#
function Rotation_Matrix(theta::Float64)
        c = cos(theta)
        s = sin(theta)
        [c s 0.0 0.0 ; 0.0 0.0 c s]
end

#
# Assembly the global stiffness matrix
#
function Global_Stiffness( ne::Int64, connect::Array{Int64,2}, coord::Array{Float64,2},  E::Array{Float64,1},  A::Array{Float64,1})

    # Define vectors I, J and V for the sparse assembly
    # The number of entries will be ne * (4*4)
    nent = ne*4*4
    I = zeros(Int64,nent)
    J = zeros(Int64,nent)
    V = zeros(nent)

    # Global counter for I,J and V
    counter = 0

    # For each element
    for ele=1:ne

       # Evaluates L, theta and dofs
       Le, theta, dofs = Ele_Data(ele,coord,connect)

       # Extract both E and A for this element
       Ee = E[ele]
       Ae = A[ele]

       # Assemble the local stiffness
       ke = Local_Stiffness(Ee,Ae,Le)

       # Assemble the rotation matrix
       Re = Rotation_Matrix(theta)

       # Rotate the local element matrix to the global system
       Ke = Re'*ke*Re

       # Loop to assamble I,J and V with the information of this element
       for i=1:4
           dof_i = dofs[i]
           for j=1:4
               dof_j = dofs[j]
               counter += 1
               I[counter] = dof_i
               J[counter] = dof_j
               V[counter] = Ke[i,j]
           end #j
      end #i

  end # ele

  # Assembly the global sparse matrix and return
  return sparse(I,J,V)

end


#
# Assembly the load vector
#
function Load_Vector( nn::Int64, nnbc::Int64, nbc::Array{Float64,2} )

     # Declare the force vector
     F = zeros(2*nn)

     # Loop over nbc and stores the information
     for i=1:nnbc

         # Find the node, the dof and the value
         node = nbc[i,1]
         dof  = nbc[i,2]
         val  = nbc[i,3]

         # Since the array ebc is double precision, all operations with node and
         # dof will also be double precision. Thus, in order to acess an
         # array position, we have to convert it to integer.
         pos = Int(2*(node-1)+dof)

         # Add valur to this position
         F[pos] += val

     end

     # Return the force vector
     return F

end


#
# Apply the essential boundary conditions (just homogeneous by now)
# by zeroing the row and the collumn, with 1.0 in the diagonal. This is NOT
# THE PROPER WAY TO DEAL WITH THE PROBLEM but this is a first example.
#
function Apply_EBC!(nebc::Int64,ebc::Array{Int64,2},K,F)

    # Loop for the natural boundary conditions
    for i=1:nebc

        # Find the node and the dof
        node = ebc[i,1]
        dof = ebc[i,2]

        # evaluates the global position
        pos = Int(2*(node-1)+dof)

        # Apply 0.0 to row and column (equation)
        K[pos,:] .= 0.0
        K[:,pos] .= 0.0
        F[pos] = 0.0

        # Set the diagonal to 1.0
        K[pos,pos] = 1.0

    end

end


#
# Plot the structure using lines
#
function Plot_Structure(ne::Int64, connect::Array{Int64,2}, coord::Array{Float64,2}, U=[],dscale=0.0)

    # Starts the plot
    display(plot());

    # Loop over the elements, extracting the nodes
    # and their coordinates
    for ele=1:ne

        # recover the nodes
        node1, node2 = connect[ele,:]
        #println(node1," ",node2)

        # recover the coordinates
        x1,y1 = coord[node1,:]
        x2,y2 = coord[node2,:]

        # Adds a line to the plot
        display( plot!([x1,x2],[y1,y2],  color=:blue,legend=false) );

    end #ele

    # If the displacement vector is given
    if length(U)>0


        # Loop over the elements. Now, the coordinates of each
        # node are given by the initial coordinates
        # plus the displacement*dscale
        for ele=1:ne

            # recover the nodes
            node1, node2 = connect[ele,:]
            #println(node1," ",node2)

            # recover the coordinates
            x1,y1 = coord[node1,:]
            x2,y2 = coord[node2,:]

            # Add the displacements
            pos = 2*(node1-1)+1; x1 = x1 + U[pos]*dscale
            pos = 2*(node1-1)+2; y1 = y1 + U[pos]*dscale
            pos = 2*(node2-1)+1; x2 = x2 + U[pos]*dscale
            pos = 2*(node2-1)+2; y2 = y2 + U[pos]*dscale

            # Adds a line to the plot
            display(plot!([x1,x2],[y1,y2], color=:red,legend=false));

        end #ele

    end #if length

end


#
# Truss
#
function Analysis(ne,nn,nebc,nnbc,connect,coord,ebc,nbc,E,A)

    print("1.1: Assembly the global force vector\n")
    F = Load_Vector(nn,nnbc,nbc)

    print("1.2: Assembly the global Stiffnes Matrix\n")
    K = Global_Stiffness(ne,connect,coord,E,A)

    print("1.3: Apply the essential boundary conditions\n")
    Apply_EBC!(nebc,ebc,K,F)

    print("1.4: Use Cholesky to solve the linear system\n")
    C = cholesky(K)

    print("--- solve for displacement U = C/F \n")
    U = C\F

end



function Truss()

    # Input Data
    ne = 6
    nn = 4
    nebc = 3
    nnbc = 2
    coord = [0.0 0.0 ;
             1.0 0.0 ;
             0.0 0.5 ;
             1.0 0.5]

    connect = [1  2 ;
               2  3 ;
               3  4 ;
               4  1 ;
               1  3 ;
               4  2 ]

    ebc = [1 1 ;
           1 2 ;
           3 1 ]
    nbc = [2 1 300.0 ;
           3 2 -200.0]

    E = 210E9*ones(ne)
    A = 1E-4*ones(ne)

    print("1. Performs the analysis\n")
    U = Analysis(ne,nn,nebc,nnbc,connect,coord,ebc,nbc,E,A)

    print("2. Plots the undeformed and the deformed strucure\n")
    Plot_Structure(ne,connect,coord,U,1000.0)

end



Truss()
#gui()
readline()