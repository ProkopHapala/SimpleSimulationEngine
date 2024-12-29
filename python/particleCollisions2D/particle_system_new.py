import numpy as np
import matplotlib.pyplot as plt

CONST_COULOMB = 14.3996448915

def getMorseQ( d, R0, E0, Qij, alpha=1.7 ):
    r    = np.linalg.norm(d)
    ir   = 1.0/r

    E_Coulomb = CONST_COULOMB * Qij * ir
    f_Coulomb = d * ( ir*irE_Coulomb ) 

    e         = np.exp(-alpha * (r - R0))
    E_Morse   =            E0 *           (e*e - 2*e)
    f_Morse   = d * ( ir * E0 * 2*alpha * (e*e -   e) ) 
    
    return f_Coulomb + f_Morse , E_Coulomb + E_Morse

def rhs_ij( pi, pj, k, l0 ):
    '''
    b_i = \sum_j (K_{ij} p_{ij})
    where p_{ij} is ideal position of particle i which satisfy the constraint between i and j
    p_{ij} = p_j + (p_i-p_j)/|p_i-p_j| * l0
    '''
    d = pi - pj
    l = np.linalg.norm(d)
    return d * (k * l0/l)

def make_PD_RHS( ps, neighs, kngs, l0s, masses=None, dt=1.0, b=None ):
    '''
    b_i = \sum_j (K_{ij} p_{ij})  + M/dt^2 p'_i
    '''
    n = len(ps)
    if b is None:
        b = np.zeros((n,2))
    inv_dt2 = 1.0/(dt*dt);
    for i in range(n):
        pi = ps[i]
        bi = np.zeros(2)
        ngsi  = neighs[i]
        ksi   = kngs  [i]
        l0i   = l0s   [i]
        ni    = len(ngsi) 
        for jj in range(ni):
            j  = ngsi[jj]
            if j<0: break
            k  = ksi [jj]
            l0 = l0i [jj] 
            pj = ps  [j ]
            bi += rhs_ij( pi, pj, k, l0 )
        if masses is not None:
            bi += masses[i] * inv_dt2    # for the moment we neglect the inertial term
        b[i] = bi
    return b

def update_PD_matrix( apos, neighs, kngs, dt=1.0 ):
    """
    Update the diagonal matrix for Projective Dynamics
    """
    n = len(apos)
    Aii = np.zeros(n)
    for ip in range(n):
        # Inertial term: m_i/dt^2
        Aii[ip] = 1.0 / (dt**2)  # Using unit mass
        # Sum of stiffness for all neighbors
        for j, k in zip(neighs[ip], kngs[ip]):
            Aii[ip] += k
    return Aii

def jacobi_iteration_sparse(ps, b, neighs, kngs, Aii, ps_out=None  ):
    """One iteration of Jacobi method using sparse operations"""
    n      = len(ps)
    if ps_out is None:
        ps_out = np.zeros_like(ps)
    #r    = np.zeros_like(x)
    err2  = 0.0 
    for i in range(n):
        #print("CPU i: %i Aii[i]: %f b[i]: %f " %(i, Aii[i], b[i]) );
        sum_j  = 0  # RHS term
        ngsi = neighs[i]
        ksi  = kngs[i] 
        ni = len(ngsi)
        for jj in range(ni):
            j      = ngsi[jj]
            if j<0: break
            k      = ksi[jj]
            sum_j -= k * ps[j]                 # Off-diagonal contribution
        ps_out[i] =  (b[i] - sum_j) / Aii[i]   # solution x_new = (b - sum_(j!=i){ Aij * x[j] } ) / Aii
        y_i       = Aii[i]*ps_out[i] + sum_j   # y = Ax
        ri = b[i] - y_i                        # Residual r = b - Ax
        err2 += ri*ri
        #print("i: %i r_i: %f pi_: %f x_i: %f " %(i, r[i], Api_/Aii[i], x[i] ) )
    return ps_out, err2

def jacobi_iteration_sparse_fly( ps, ps_out, neighs, kngs, l0s, masses=None, dt=1.0 ):
    """
    This version of Jacobi solver updates b and Aii on the fly
    """
    n = len(ps)
    #ps_out = np.zeros_like(ps)
    #r      = np.zeros_like(ps)
    err2  = 0.0 
    inv_dt2 = 1.0/(dt*dt);
    for i in range(n):
        #print("CPU i: %i Aii[i]: %f b[i]: %f " %(i, Aii[i], b[i]) );
        sum_j = 0     # RHS term
        mi    = masses[i] * inv_dt2   # intertial term M_i/dt^2
        bi    = mi*ps[i]              # bi  = \sum_j (K_{ij} d_{ij})  + M_i/dt^2 p'_i
        Aii   = mi                    # Aii = \sum_j (K_{ij} +          M_i/dt^2
        
        pi    = ps[i] 
        ngsi  = neighs[i]
        ksi   = kngs[i] 
        l0i   = l0s[i]
        ni    = len(ngsi)
        for jj in range(ni):
            j      = ngsi[jj]
            if j<0: break
            k      = ksi[jj]
            sum_j -= k * ps_in[j]   # Off-diagonal contribution

            l0   = l0i[jj]
            pj   = ps_in[j]
            bi  += rhs_ij( pi, pj, k, l0 ) # NOTE: updating bi on the fly here goes beyond linear solver
            Aii += k

        ps_out[i] =  (bi - sum_j) / Aii[i]   # solution x_new = (b - sum_(j!=i){ Aij * x[j] } ) / Aii
        y_i  = Aii[i]*pi + sum_j             # This is Ax_i
        ri   = bi - y_i                      # Residual r = b - Ax
        err2 += ri*ri
        #print("i: %i r_i: %f pi_: %f x_i: %f " %(i, r[i], Api_/Aii[i], x[i] ) )
    return ps_out, err2


def jacobi_iteration_sparse_debug(ps, b, neighs, kngs, Aii, bPlot=True, l0s=None):
    """Debug version of Jacobi iteration that plots p'_ij points"""

    print("jacobi_iteration_sparse_debug")
    print("b: ",      b      )
    print("neighs: ", neighs )
    print("kngs: ",   kngs   )
    print("Aii: ",    Aii    )
    
    n = len(ps)
    ps_out = np.zeros_like(ps)
    r = np.zeros_like(ps)
    
    if bPlot: 
        plt.scatter(ps[:,0], ps[:,1], c='blue', label='Current positions')
    
    # For each particle
    for i in range(n):
        pi    = ps[i]
        bi    = b [i]    #   bi = \sum_j (K_{ij} d_{ij})
        ngsi  = neighs[i]
        ksi   = kngs[i]
        ni    = len(ngsi)
        sum_j = np.zeros(2)  

        sum_pij = np.zeros(2)  
        for jj in range(ni):
            j  = ngsi[jj]
            if j<0: break
            k  = ksi [jj]
            pj = ps  [j]
            sum_j += k * pj
            if bPlot:
                l0   = l0s[i][jj]
                d_ij = rhs_ij( pi, pj, 1.0, l0 )
                p_ij = pj + d_ij

                sum_pij += p_ij * k
                plt.plot(p_ij[0], p_ij[1], '.r', alpha=0.5)                     # Plot p'_ij point
                plt.plot([pj[0], p_ij[0]], [pj[1], p_ij[1]], 'k--', lw=0.5, alpha=0.3 ) # Draw line from current position to p'_ij
                
        # Note 
        # d_{ij}  = (p_i - p_j)/|p_i-p_j| * l0
        # bi      = \sum_j K_{ij} d_{ij}
        # sum_j   = \sum_j K_{ij} p_j
        # p'_{ij} = p_j + d_{ij}
        # A_{ii}  = \sum_j K_{ij}
        # p^{new}_i = (bi + sum_j) / A_{ii} 
        #           = \sum_j K_{ij} (d_{ij} + p_j) / A_{ii} 
        #           = K_{ij} p'_{ij} / A_{ii} 
        #           = {sum_j K_{ij} p'_{ij} } / { sum_j K_{ij} }

        print(f"i: {i} b_i: {bi} sum_j: {sum_j} Aii[i]: {Aii[i]}")

        pi_new = (bi + sum_j) / Aii[i]

        pi_ij = sum_pij / Aii[i]

        if bPlot:
            plt.plot(pi_new[0], pi_new[1], '+k', alpha=0.5)                     # Plot p'_ij point
            plt.plot([pi[0], pi_new[0]], [pi[1], pi_new[1]], 'k--', lw=0.5, alpha=1.0 ) # Draw line from current position to p'_ij

            plt.plot(pi_ij[0], pi_ij[1], '+r', alpha=0.5)  
            plt.plot([pi[0], pi_ij[0]], [pi[1], pi_ij[1]], 'r--', lw=0.5, alpha=1.0 ) 
        
        ps_out[i] = pi_new

        # Calculate residual
        y_i = Aii[i]*ps[i] + sum_j  # This is Ax_i
        r[i] = b[i] - y_i  # Residual r = b - Ax
        #print(f"i: {i} r_i: {np.linalg.norm(r[i]):.6f} pi_: {np.linalg.norm(b[i]/Aii[i]):.6f}")
    
    # Plot bonds between particles
    # for i in range(n):
    #     for j in neighs[i]:
    #         plt.plot([ps[i,0], ps[j,0]], [ps[i,1], ps[j,1]], 'k-', lw=0.5, alpha=1.0 )
    
    if bPlot:
        plt.axis('equal')
        plt.grid(True)
        plt.legend()
        plt.title('Debug visualization of Jacobi iteration')
        #plt.show()
    
    err2 = np.linalg.norm(r)
    return ps_out, err2

class Particle:
    def __init__(self, pos, vel=None, mass=1.0):
        self.pos = np.array(pos, dtype=np.float64)    # position vector [x,y]
        self.vel = np.zeros(2) if vel is None else np.array(vel, dtype=np.float64)   # velocity vector [vx,vy]
        self.m = mass                                 # mass
        self.f = np.zeros(2)                         # force accumulator
        self.neighbors = []                          # list of neighbor particles for collision detection
        
class Bond:
    """A bond between two particles"""
    def __init__(self, i, j, k=1.0, r0=1.0):
        self.i = i   # index of first particle
        self.j = j   # index of second particle
        self.k = k   # spring constant
        self.r0 = r0 # rest length
    
    def get_error(self, pos):
        """Get the error in the bond constraint"""
        d = pos[self.j] - pos[self.i]
        r = np.linalg.norm(d)
        return abs(r - self.r0)

class Molecule:
    def __init__(self, types, pos, bonds=None):
        self.types = np.array(types, dtype=np.float64)  # particle masses
        self.pos = np.array(pos, dtype=np.float64)      # positions
        self.bonds = bonds if bonds else []             # bonds between particles
        self.vel = np.zeros_like(self.pos)             # velocities (initially zero)
        self.charges = np.zeros(len(types))            # charges (initially zero)
        # -- to be done
        self.cog   = None            # np.array[3]
        self.inds  = np.zeros(len(pos), dtype=np.int32)   # list of indices of this p
        
    def update(self, bUpdateCenter):
        """Update center of geometry and radius of the molecule"""
        # Convert pos to numpy array if it's a list
        pos = np.array(self.pos)
        
        # Compute center of geometry
        if bUpdateCenter or self.cog is None:
            self.cog = np.mean(pos, axis=0)
        self.Rg = np.max(np.linalg.norm(pos - self.cog, axis=1))

class World:
    def __init__(self, dt=0.1, D=0.1):
        nmaxneighs = 8
        self.dt = dt
        self.D = D           # collision offset
        self.Kbond = 20.0    # bond stiffness
        self.K_col = 20.0    # collision stiffness
        self.damping = 0.1   # velocity damping
        
        # Non-covalent interaction parameters
        self.eps    = 0.1    # Morse potential well depth
        self.r0_vdw = 0.5    # Morse potential equilibrium distance
        self.alpha  = 2.0    # Morse potential decay rate
        self.k_el   = 1.0    # Coulomb constant (1/4πε０)
        
        # Arrays for positions, velocities, and forces
        self.apos     = None # current positions
        self.vels     = None # velocities
        self.forces   = None # forces
        self.pos_pred = None # predicted positions
        self.pos_solv = None # position buffer used for costraint solver
        
        self.masses   = None # particle masses
        self.charges  = None # particle charges
        
        # Neighbor lists and parameters
        self.neighs = []        # neighbor indices
        self.kngs   = []        # neighbor stiffness
        self.l0s    = []        # rest lengths
        self.neighs_bonds = []  # bonded neighbors
        
        # simulation parameters
        self.Rc     = 1.0       # cutoff radius for non-bonded interactions
        self.Rg     = 2.0       # grouping radius for neighbor search
        self.n_iter = 10        # number of iterations for constraint solver

    def allocate(self, natoms):
        self.apos      = np.zeros((natoms,2))
        self.pos_pred  = np.zeros((natoms,2))
        self.pos_solv  = np.zeros((natoms,2))
        self.vels      = np.zeros((natoms,2))
        self.forces    = np.zeros((natoms,2))
        self.neighs_bonds = [set() for _ in range(natoms)]
        self.neighs       = [[]    for _ in range(natoms)]
        self.kngs         = [[]    for _ in range(natoms)]
        self.l0s          = [[]    for _ in range(natoms)]
        self.b            = np.zeros(natoms)
        self.Aii          = np.zeros(natoms)
        self.masses       = np.ones(natoms)
        self.charges      = np.zeros(natoms)
    
    def add_molecule(self, molecule):
        """Add a molecule to the world"""
        self.molecules.append(molecule)
        self.from_molecules(self.molecules)

    def from_molecules(self, molecules):
        """Initialize world from a list of molecules"""
        # Count total number of atoms
        natoms = sum(len(mol.pos) for mol in molecules)
        self.allocate(natoms)
        
        # Add each molecule
        offset = 0
        for mol in molecules:
            n = len(mol.pos)
            self.apos[offset:offset+n] = mol.pos
            self.vels[offset:offset+n] = mol.vel
            self.masses[offset:offset+n] = mol.types
            if hasattr(mol, 'charges'):
                self.charges[offset:offset+n] = mol.charges
            
            # Add bonds to neighbor list
            for bond in mol.bonds:
                i = bond.i + offset
                j = bond.j + offset
                self.neighs_bonds[i].add((j, bond.r0))
                self.neighs_bonds[j].add((i, bond.r0))
            
            offset += n

    def update_groups(self, bUpdateCenter=False):
        '''
        Goes through atoms of molecule and find the center (optionaly) and group radius Rg
        '''
        for mol in self.molecules: mol.update(bUpdateCenter)
    
    def update_group_neighbors(self):
        """
        this is to accelerate the neighbor search. Each molecule keep its list of atoms which are closer than r<Rcg where Rcg=Rc+Rd and Rg is radius of the group (distance of furthers member atom from center of the group)  
          - this list is update less ofthen for performace
        """
        for im,mol in enumerate(self.molecules):
            cog  = mol.cog
            gngs = mol.mol_neighs
            for jm,molj in enumerate(self.molecules):
                if im == jm: continue
                for jp in molj:
                    if np.linalg.norm(self.apos[jp] - cog) < molj.Rg:
                        gngs[jp].add(im)

    def update_collision_neighbors(self, bCollisions=True):
        """
        Update neighbors list with both bonded and collision neighbors.
        Args:
            bCollisions (bool): If True, include collision neighbors
        """
        print( f"update_collision_neighbors() bCollisions={bCollisions}" )
        n = len(self.apos)
        self.neighs[:] = [None] * n
        self.kngs  [:] = [None] * n
        self.l0s   [:] = [None] * n
        # First add all bonded neighbors
        for i in range(n):
            neighs = []
            kngs   = []
            l0s    = []
            Aii    = 0.0 
            #print(f"Bond neighbors for atom {i}: {self.neighs_bonds[i]}")  # Debug
            for j, r0 in self.neighs_bonds[i]:
                k     =  self.Kbond
                Aii  += k
                kngs  .append(k)
                l0s   .append(r0)
                neighs.append(j)
                #print(f"  Added bond neighbor {i}-{j} with r0={r0}")  # Debug
            # Add collision neighbors if enabled
            if bCollisions:
                for j in range(n):  # we loop over all atoms, not just j>i because it is better for parallelization 
                    if i == j: continue
                    if any(j == x[0] for x in self.neighs_bonds[i]):   # this can be slow, need to optimize (?)
                        continue
                    # Check if within cutoff
                    d = self.apos[j] - self.apos[i]
                    r = np.linalg.norm(d)
                    if r < self.Rc:
                        # Add collision pair (both directions)
                        k    = self.K_col   # perhaps later we can use different stiffnesses for each particle by some mixing ?                        
                        Aii += k
                        kngs  .append(k)
                        l0s   .append(-self.D)    # Negative rest length for collisions
                        neighs.append(j)
            self.neighs[i] = neighs ;
            self.kngs  [i] = kngs   ;
            self.l0s   [i] = l0s    ;
            self.Aii   [i] = Aii    ;
        return self.neighs, self.kngs, self.l0s


    # --------- Iterative Solver  -------------
    def step(self):
        """Perform one time step of the simulation"""
        # Damp velocities
        self.vels[:] *= (1.0 - self.damping)
        self.predict_positions()
        self.solve_constraints()
        self.update_velocities()

    def predict_positions(self):
        """Update predicted positions using leapfrog with non-covalent forces"""
        # Compute non-covalent forces
        self.compute_non_covalent_forces()
        self.vels    [:,:] += self.forces[:,:] * self.dt /        self.masses[:, np.newaxis]
        self.pos_pred[:,:]  = self.apos[:,:]   + self.vels[:,:] * self.dt
        
    def update_velocities(self, pos_pred ):
        self.vels[:,:] = (pos_pred[:,:] - self.pos[:,:])/self.dt
        self.apos[:,:] =  pos_pred[:,:]

    def solve_constraints(self, niter=5, bOnTheFly=False):
        """Solve position constraints using Jacobi iteration"""
        # Update neighbors and constraints
        self.update_collision_neighbors()
        ps_in  = self.pos_pred
        ps_out = self.pos_solv
        if not bOnTheFly:
            b   = self.update_rhs(ps_in, b=self.b  )
            #Aii = self.update_Aii(ps_in, Aii=self.Aii ) # Aii is already updated in update_collision_neighbors()
        for it in range(niter):
            if bOnTheFly:
                ps_out, err2 = jacobi_iteration_sparse_fly(ps_in, ps_out, self.neighs, self.kngs, self.l0s, masses=self.masses, dt=self.dt)
            else:
                ps_out, err2 = jacobi_iteration_sparse(ps, b, self.neighs, self.kngs, Aii, ps_out=ps_out)
            if err2 < 1e-8: break
            ps_in, ps_out = ps_out, ps_in

        self.pos_pred[:,:] = ps
        return ps

    def update_rhs(self, pos_pred ):
        return make_PD_RHS( pos_pred, self.neighs, self.kngs, self.l0s, masses=self.masses, dt=1.0 )

    def update_Aii(self, Aii=None ):
        """Update diagonal elements of PD matrix"""
        n = len(self.apos)
        if Aii is None: Aii = np.zeros(n)
        Aii += 1.0 / (self.dt * self.dt)  # m_i/dt^2 (using unit mass)
        for i in range(n):
            for k in self.kngs[i]:
                Aii[i] += k
        return Aii

    def print_bond_lengths(self, end=" "):
        for i in range(len(self.apos)):
            for j in self.neighs[i]:
                if j > i:  # Only print each distance once
                    d = np.linalg.norm(self.apos[j] - self.apos[i])
                    print(f"{i}-{j}: {d:.6f},", end=end )
        #print()

    def compute_non_covalent_forces(self):
        """Compute forces from non-covalent interactions (Morse + Coulomb)"""
        n = len(self.apos)
        self.forces[:,:] = 0.0
        E = 0.0
        for i in range(n):
            pi = self.apos  [i]
            fi = self.forces[i]
            for j in range(n):
                if i == j: continue
                if j in self.neighs_bonds[i]:  # this can be slow, need to optimize (?)
                    continue
                d = self.apos[j] - pi
                f,E = getMorseQ( d, R0, E0, Qij, alpha=1.7 )
                fi += f
                self.forces[j] -= f
        return E

# ========================================================
# ============ Create initial configurations  ============
# ========================================================

def create_triangle_system(l=1.0):
    """Create three particles in equilateral triangle with bonds"""
    # Create three particles at vertices of equilateral triangle
    pos = [
        [0.0, 0.0],                # First particle at origin
        [l, 0.0],                  # Second particle at distance l along x-axis
        [l/2, l*np.sqrt(3)/2]      # Third particle to form equilateral triangle
    ]
    types = [1.0, 1.0, 1.0]  # Equal masses
    bonds = [
        Bond(0, 1, r0=l),  # Bond between particles 0 and 1
        Bond(1, 2, r0=l),  # Bond between particles 1 and 2
        Bond(2, 0, r0=l)   # Bond between particles 2 and 0
    ]
    return Molecule(types, pos, bonds)

def create_water(pos=[0,0], l0=1.0, ang=np.pi/2, vel=[0,0]):
    """Create a water molecule (H2O) at given position with given velocity"""
    # Create oxygen atom at center
    types = [16.0, 1.0, 1.0]  # O, H, H masses
    charges = [-0.8, 0.4, 0.4]  # O, H, H charges
    
    # Calculate H positions
    r = l0  # O-H bond length
    dx = r * np.cos(ang/2)
    dy = r * np.sin(ang/2)
    
    positions = [
        pos,                    # O at center
        [pos[0]-dx, pos[1]-dy], # H1
        [pos[0]+dx, pos[1]-dy]  # H2
    ]
    
    # Create bonds
    bonds = [
        Bond(0, 1, r0=r),  # O-H1 bond
        Bond(0, 2, r0=r)   # O-H2 bond
    ]
    
    mol = Molecule(types, positions, bonds)
    mol.charges = charges  # Add charges to molecule
    mol.vel = np.array(vel)
    return mol

def create_water_grid(nx=3, ny=3, spacing=2.0, l0=1.0, ang=np.pi/2):
    """Create a grid of water molecules"""
    molecules = []
    for i in range(nx):
        for j in range(ny):
            pos = [i*spacing, j*spacing]
            molecules.append(create_water(pos, l0=l0, ang=ang))
    return molecules

def init_world(molecule_func=create_triangle_system, dt=0.05, D=0.1, **kwargs):
    """Initialize world with molecules created by the provided function
    Args:
        molecule_func: Function to create molecules. Default is create_triangle_system
        dt: Time step for simulation
        D: Collision offset
        **kwargs: Additional arguments passed to molecule_func
    """
    # Create world with specified parameters
    world = World(dt=dt, D=D)
    
    # Create molecules using the provided function
    if molecule_func == create_water_grid:
        # Special case for water grid which returns multiple molecules
        mols = molecule_func(**kwargs)
    else:
        # Single molecule case
        mols = [molecule_func(**kwargs)]
    
    # Initialize world with molecules
    world.from_molecules(mols)
    return world

# ========================================================
# ============ Create initial configurations  ============
# ========================================================

def test_jacobi_debug(world, n_iter=5, bPrint=True, bPlot=True, sz=5, fRnadom=0.3, bUpdateB=True, rseed = 454454 ):
    """Test convergence with three particles and debug visualization"""
    print("\nTesting three-particle system convergence with debug:")
    
    world.update_collision_neighbors(bCollisions=False)  # Comment this out

    # random perturbation to points
    np.random.seed(rseed)
    world.apos[:,:] += np.random.uniform(-fRnadom, fRnadom, world.apos.shape)

    #b = make_PD_RHS(world.apos, world.neighs, world.kngs, world.l0s, masses=world.masses, dt=world.dt)
    b = make_PD_RHS(world.apos, world.neighs, world.kngs, world.l0s )

    print("##### Initial world state: \n")
    for i in range(len(world.apos)):
        print( f"i: {i} Aii: {world.Aii[i]} b: {b[i]} ")  

    if bPlot:
        plt.figure(figsize=(n_iter*sz,sz))

    # Solve iteratively the equation Ap = b, where A is the PD matrix, p is the position, and b is the right-hand side from make_PD_RHS
    ps = world.apos  # Only use x,y components
    for it in range(n_iter):

        if bPlot:
            ax = plt.subplot(1, n_iter, it+1)

        if bUpdateB:
            # NOTE: this goes beyond linear-algebra solution, in poper linear algebra we should not update b during iterative solution of Ap=b    
            #       but with this it converge faster, and it is more physical
            b = make_PD_RHS(world.apos, world.neighs, world.kngs, world.l0s )

        ps_new, r = jacobi_iteration_sparse_debug(ps, b, world.neighs, world.kngs, world.Aii, l0s=world.l0s, bPlot=bPlot)
        ps[:] = ps_new[:]
        if bPrint:
            print(f"\nIteration {it} Distances:", end=" " )
            world.print_bond_lengths( end=" " )

def test_triangle_cl():
    """Compare Python and OpenCL implementations for triangle system"""
    # Initialize world with three particles
    world = init_world()
    n = len(world.apos)

    import constrains_cl
    
    solver = constrains_cl.CLConstrains(max_points=n, max_neighs=2)
    
    # Prepare neighbor data (only bonds, no collisions)
    world.update_collision_neighbors(bCollisions=False)
    
    # Print and perturb initial positions
    print("Initial positions:")
    print(world.apos)
    world.apos += np.random.rand(n,2) * 0.3
    print("\nPerturbed positions:")
    print(world.apos)
    
    # Parameters for both solvers
    dt = 0.1
    inv_dt2 = 1.0/(dt*dt)
    Rd = 0.1
    n_iter = 10
    
    # Run both solvers
    print("\nRunning Python solver...")
    world.solve_constraints(niter=n_iter)
    pos_py = world.apos.copy()
    
    print("\nRunning OpenCL solver...")
    pos_cl = solver.solve_constraints(world.apos, np.ones(n), world.neighs, world.kngs, world.l0s, inv_dt2, Rd, niter=n_iter)
    
    # Compare results
    diff = np.linalg.norm(pos_py - pos_cl)
    print(f"\nDifference between solutions: {diff}")
    
    # Visualize results
    plt.figure(figsize=(12,4))
    plt.subplot(131);  pu.plot_molecules(world, ax=plt.gca(), color='blue',  label='Initial'); plt.title('Initial')
    world.apos = pos_py; plt.subplot(132);  pu.plot_molecules(world, ax=plt.gca(), color='red',   label='Python');  plt.title('Python Solution')
    world.apos = pos_cl; plt.subplot(133);  pu.plot_molecules(world, ax=plt.gca(), color='green', label='OpenCL');  plt.title('OpenCL Solution')
    plt.tight_layout()

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import plot_utils as pu
    
    # Create a system with two water molecules
    #world = init_world(molecule_func=create_water_grid, nx=2, ny=1, spacing=1.5, l0=0.5, ang=104.5*np.pi/180)
    world = init_world()

    
    # Run test with debug visualization
    test_jacobi_debug(world, n_iter=5)
    plt.show()