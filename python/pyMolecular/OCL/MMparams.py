import numpy as np
import math

# Constants (Define these based on your simulation requirements)
Lepair = 1.5    # Example bond length for electron pairs
Kepair = 100.0  # Example bond stiffness for electron pairs
deg2rad = np.pi / 180.0

class ElementType:
    """
    Represents the parameters for each element type based on ElementTypes.dat.
    
    Fields follow the definition from ElementTypes.dat:
    #1        2   3   4   5 6        7         8               9       10               11                    12        13      14        15
    #name     iZ  ne  nv  npi color     Rcov[Ang] RvdW[Ang] EvdW[eV]        Quff[e] Uuff[eV]        Vuff[eV]       QEq{ Eaff[eV] Ehard[eV] Ra[Ang] eta[a.u.] } 
    """
    def __init__(self, name, iZ=0, neval=0, valence=0, piMax=0, color=0, Rcov=0.0, RvdW=0.0, EvdW=0.0, 
                 Quff=0.0, Uuff=0.0, Vuff=0.0, bQEq=False, Eaff=0.0, Ehard=0.0, Ra=0.0, eta=0.0):
        """
        Initializes the ElementType with parameters from ElementTypes.dat.
        Optional parameters are set to default values if not provided.
        """
        self.name = name        # Element symbol
        self.iZ = iZ            # Atomic number
        self.neval = neval      # Number of valence electrons
        self.valence = valence  # Sum of bond orders of all bonds
        self.piMax = piMax      # Max number of pi orbitals
        self.color = color      # Color code (hex)
        self.Rcov = Rcov        # Covalent radius
        self.RvdW = RvdW        # van der Waals radius
        self.EvdW = EvdW        # van der Waals energy well
        self.Quff = Quff        # Effective charge in UFF
        self.Uuff = Uuff        # UFF parameter
        self.Vuff = Vuff        # UFF parameter
        self.bQEq = bQEq        # Flag for QEq parameters
        self.Eaff = Eaff        # Electronegativity
        self.Ehard = Ehard      # Chemical hardness
        self.Ra = Ra            # Atomic size
        self.eta = eta          # Valence orbital exponent

class AtomType:
    """
    Represents the parameters for each atom type based on AtomTypes.dat.
    
    Fields follow the definition from AtomTypes.dat:
    #1     2     3       4    5  6   7   8    9       10      11       12      13       15     16      17  18  19  20
    #name parent element epair nv ne npi sym  Ruff   RvdW    EvdW       Qbase    Hb      Ass    Asp     Kss Ksp Kep Kpp
    """
    def __init__(self, name, parent_name="*", element_name="", epair_name="", valence=0, nepair=0, npi=0, sym=0,
                 Ruff=0.0, RvdW=0.0, EvdW=0.0, Qbase=0.0, Hb=0.0, bMMFF=False, Ass=0.0, Asp=0.0, Kss=0.0, Ksp=0.0, Kep=0.0, Kpp=0.0,
                 element=-1, parrent=-1, ePairType=-1, iZ=0, color=0, subTypes=None):
        """
        Initializes the AtomType with parameters from AtomTypes.dat.
        Optional parameters are set to default values if not provided.
        """
        self.name = name                 # Atom type name
        self.parent_name = parent_name   # Parent type name
        self.element_name = element_name # Element type name
        self.epair_name = epair_name     # Electron pair type name
        self.valence = valence           # Sum of bond orders of all bonds
        self.nepair = nepair             # Number of electron pairs
        self.npi    = npi                # Number of pi orbitals
        self.sym    = sym                # Symmetry type (tetrahedral, triangle, etc.)
        self.Ruff   = Ruff               # UFF natural bond radius
        self.RvdW   = RvdW               # LJ distance parameter
        self.EvdW   = EvdW               # LJ energy parameter
        self.Qbase  = Qbase              # Atomic charge
        self.Hb     = Hb                 # Hydrogen bond correction
        self.bMMFF  = bMMFF              # Flag for MMFF parameters
        self.Ass    = Ass                # Equilibrium angle for sigma-sigma
        self.Asp    = Asp                # Equilibrium angle for sigma-pi
        self.Kss    = Kss                # Force constant for sigma-sigma
        self.Ksp    = Ksp                # Force constant for sigma-pi
        self.Kep    = Kep                # Force constant for pi-lone pair
        self.Kpp    = Kpp                # Force constant for lone pair-lone pair
        
        # References to other types (to be filled after all types are loaded)
        self.element = element         # Index of corresponding element type
        self.parrent = parrent         # Index of parent type
        self.ePairType = ePairType     # Type of lone pair owned
        self.iZ = iZ                   # Atomic number from element
        self.color = color             # Color from element
        self.subTypes = subTypes if subTypes is not None else [0, 0, 0]  # sp1, sp2, sp3

def read_element_types(filepath):
    """
    Read element types from ElementTypes.dat file.
    
    Parameters:
    - filepath (str): Path to the ElementTypes.dat file
    
    Returns:
    - dict: Dictionary mapping element names to ElementType objects
    """
    print("read_element_types() filepath=", filepath)
    element_types = {}
    with open(filepath, 'r') as f:
        lines = f.readlines()
        
        for line in lines:
            line = line.strip()
            # Skip comments and empty lines
            if line.startswith('#') or not line:
                continue
                
            parts = line.split()
            if not parts:
                continue
                
            name = parts[0]
            et = ElementType(name=name)
            
            # Parse required parameters (minimum 12 values as per C++ code)
            try:
                if len(parts) >= 2:
                    et.iZ = int(parts[1])
                if len(parts) >= 3:
                    et.neval = int(parts[2])
                if len(parts) >= 4:
                    et.valence = int(parts[3])
                if len(parts) >= 5:
                    et.piMax = int(parts[4])
                if len(parts) >= 6:
                    et.color = int(parts[5], 16) if parts[5].startswith('0x') else int(parts[5])
                if len(parts) >= 7:
                    et.Rcov = float(parts[6])
                if len(parts) >= 8:
                    et.RvdW = float(parts[7])
                if len(parts) >= 9:
                    et.EvdW = float(parts[8])
                if len(parts) >= 10:
                    et.Quff = float(parts[9])
                if len(parts) >= 11:
                    et.Uuff = float(parts[10])
                if len(parts) >= 12:
                    et.Vuff = float(parts[11])
                
                # Parse optional QEq parameters
                if len(parts) >= 16:
                    et.bQEq = True
                    et.Eaff = float(parts[12])
                    et.Ehard = float(parts[13])
                    et.Ra = float(parts[14])
                    et.eta = float(parts[15])
            except (ValueError, IndexError) as e:
                print(f"Error parsing line: {line}")
                print(f"Error: {e}")
                continue
                
            element_types[name] = et
            
    return element_types

def read_atom_types(filepath, element_types=None):
    """
    Read atom types from AtomTypes.dat file.
    
    Parameters:
    - filepath (str): Path to the AtomTypes.dat file
    - element_types (dict, optional): Dictionary of element types for reference
    
    Returns:
    - dict: Dictionary mapping atom type names to AtomType objects
    """
    atom_types = {}
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
        
        # First pass: Create atom types with basic parameters
        for line in lines:
            line = line.strip()
            # Skip comments and empty lines
            if line.startswith('#') or not line:
                continue
                
            parts = line.split()
            if not parts or len(parts) < 5:  # Need at least name, parent, element, epair, valence
                continue
                
            name = parts[0]
            parent_name = parts[1]
            element_name = parts[2]
            epair_name = parts[3]
            
            at = AtomType(
                name=name,
                parent_name=parent_name,
                element_name=element_name,
                epair_name=epair_name
            )
            
            # Parse the rest of the parameters
            try:
                if len(parts) >= 5:
                    at.valence = int(parts[4])
                if len(parts) >= 6:
                    at.nepair = int(parts[5])
                if len(parts) >= 7:
                    at.npi = int(parts[6])
                if len(parts) >= 8:
                    at.sym = int(parts[7])
                if len(parts) >= 9:
                    at.Ruff = float(parts[8])
                if len(parts) >= 10:
                    at.RvdW = float(parts[9])
                if len(parts) >= 11:
                    at.EvdW = float(parts[10])
                if len(parts) >= 12:
                    at.Qbase = float(parts[11])
                if len(parts) >= 13:
                    at.Hb = float(parts[12])
                
                # MMFF parameters (optional)
                if len(parts) >= 19:
                    at.bMMFF = True
                    at.Ass = float(parts[13])
                    at.Asp = float(parts[14])
                    at.Kss = float(parts[15])
                    at.Ksp = float(parts[16])
                    at.Kep = float(parts[17])
                    at.Kpp = float(parts[18])
            except (ValueError, IndexError) as e:
                print(f"Error parsing line: {line}")
                print(f"Error: {e}")
                continue
            
            atom_types[name] = at
    
    # Second pass: Resolve references if element_types is provided
    if element_types and atom_types:
        for name, at in atom_types.items():
            # Resolve element reference
            if at.element_name in element_types:
                et = element_types[at.element_name]
                at.iZ = et.iZ
                at.color = et.color
            
            # Resolve parent reference
            if at.parent_name in atom_types and at.parent_name != "*":
                at.parrent = list(atom_types.keys()).index(at.parent_name)
            
            # Resolve epair reference
            if at.epair_name in atom_types and at.epair_name != "*":
                at.ePairType = list(atom_types.keys()).index(at.epair_name)
    
    return atom_types

def read_AtomAndElementTypes(path, felement_types='ElementTypes.dat', fatom_types='AtomTypes.dat'):
    print("read_AtomAndElementTypes() path", path)
    path1 = path + felement_types; print("path1", path1)
    path2 = path + fatom_types;    print("path2", path2)
    element_types = read_element_types(path1)
    atom_types    = read_atom_types   (path2, element_types)
    return element_types, atom_types

def generate_REQs_from_atom_types(mol, atom_types):
    """
    Generate REQs array for the molecule based on atom types.
    
    Parameters:
    - mol (AtomicSystem): The atomic system
    - atom_types (dict): Dictionary of AtomType objects
    
    Returns:
    - np.array: REQs array with shape (natoms, 4)
    """
    natom = len(mol.apos)
    REQs = np.zeros((natom, 4), dtype=np.float32)
    
    for ia in range(natom):
        atom_name = mol.enames[ia]
        if mol.qs is not None:
            q=mol.qs[ia]
        else:
            q=0
        if atom_name in atom_types:
            at = atom_types[atom_name]
            REQs[ia, 0] = at.RvdW
            REQs[ia, 1] = at.EvdW
            REQs[ia, 2] = q
            REQs[ia, 3] = at.Hb   # H-bond correction
        else:
            print(f"Warning: Atom type {atom_name} not found in atom_types")
            # Set default values
            REQs[ia, 0] = 1.0  # Default radius
            REQs[ia, 1] = 0.01  # Default energy
            REQs[ia, 2] = q
            REQs[ia, 3] = 0.0   # H-bond correction
    
    return REQs

