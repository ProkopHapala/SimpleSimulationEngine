---@meta
---This file provides Lua Language server anotations of functions declared in EditSpacecraft.cpp

-- Global vectors commonly used
---@type number[] Position vector {0,0,0}
origin = {}
---@type number[] X-axis vector {1,0,0}
xvec = {}
---@type number[] Z-axis vector {0,0,1}
zvec = {}

--- Creates a new material in C++
---@param props table Material properties
---@param props.name string Name of the material
---@param props.density number Material density
---@param props.Spull number Pull strength
---@param props.Spush number Push strength
---@param props.Kpull number Pull coefficient
---@param props.Kpush number Push coefficient
---@param props.reflectivity number Material reflectivity
---@param props.Tmelt number Melting temperature
---@return number materialId ID of the created material in C++
function Material(props) end

--- Creates a new stick material in C++
---@param name string Name of the stick material
---@param material string Base material name
---@param length number Length of the stick
---@param width number Width of the stick
---@param angle number Angle of the stick
---@return number materialId ID of the created stick material in C++
function StickMaterial(name, material, length, width, angle) end

--- Creates a new node at specified position in C++
---@param pos number[] Position vector {x,y,z}
---@return number nodeId ID of the created node in C++
function Node(pos) end

--- Creates a bound node attached to a girder or rope
---@param component number Component ID (girder/rope)
---@param position number Position along component (0.0-1.0)
---@param offset number[] Offset vector {x,y,z}
---@return number nodeId ID of the created bound node in C++
function BoundNode(component, position, offset) end

--- Creates a slider component
---@param comp1 number First component ID
---@param comp2 number Second component ID
---@param kind1 number First component kind
---@param kind2 number Second component kind
---@param nseg1 number Segments on first component
---@param nseg2 number Segments on second component
---@param offset1 number Offset on first component
---@param offset2 number Offset on second component
---@return number sliderId ID of the created slider in C++
function Slider(comp1, comp2, kind1, kind2, nseg1, nseg2, offset1, offset2) end

--- Creates a rope component
---@param node1 number First node ID
---@param node2 number Second node ID
---@param nseg number Number of segments
---@param material string Material name
---@param radius number Rope radius
---@return number ropeId ID of the created rope in C++
function Rope(node1, node2, nseg, material, radius) end

--- Creates a rope component (alternative version)
---@param nodes number[] Array of node IDs
---@param nseg number Number of segments
---@param material string Material name
---@param radius number Rope radius
---@return number ropeId ID of the created rope in C++
function Rope2(nodes, nseg, material, radius) end

--- Creates a girder component
---@param from number Starting node ID
---@param to number Ending node ID
---@param up number[] Up vector {x,y,z}
---@param nseg integer Number of segments
---@param mseg integer Number of minor segments
---@param dims number[] Width and height {w,h}
---@param material string Material name
---@param sticks number[] Array of stick material IDs
---@return number girderId ID of the created girder in C++
function Girder(from, to, up, nseg, mseg, dims, material, sticks) end

--- Creates a ring component
---@param node number Center node ID
---@param dir number[] Direction vector {x,y,z}
---@param up number[] Up vector {x,y,z}
---@param radius number Ring radius
---@param nseg number Number of segments
---@param dims number[] Width and height {w,h}
---@param material string Material name
---@param sticks number[] Array of stick material IDs
---@return number ringId ID of the created ring in C++
function Ring(node, dir, up, radius, nseg, dims, material, sticks) end

--- Creates a ring component (alternative version)
---@param girders number[] Array of girder IDs (3-4 girders)
---@param coords number[] Array of coordinates along girders
---@param center number[] Center position {x,y,z}
---@param nseg number Number of segments
---@param dims number[] Width and height {w,h}
---@param material string Material name
---@param sticks number[] Array of stick material IDs
---@return number ringId ID of the created ring in C++
function Ring2(girders, coords, center, nseg, dims, material, sticks) end

--- Creates a weld between components
---@param comp1 number First component ID
---@param comp2 number Second component ID
---@param maxRadius number Maximum radius
---@param material number Material ID
---@return number weldId ID of the created weld in C++
function Weld(comp1, comp2, maxRadius, material) end

--- Creates a radiator component
---@param girder1 number First girder ID
---@param span1_start number Start position on first girder (0.0-1.0)
---@param span1_end number End position on first girder (0.0-1.0)
---@param girder2 number Second girder ID
---@param span2_start number Start position on second girder (0.0-1.0)
---@param span2_end number End position on second girder (0.0-1.0)
---@param temperature number Operating temperature
---@return number radiatorId ID of the created radiator in C++
function Radiator(girder1, span1_start, span1_end, girder2, span2_start, span2_end, temperature) end

--- Creates a shield component
---@param girder1 number First girder ID
---@param span1_start number Start position on first girder (0.0-1.0)
---@param span1_end number End position on first girder (0.0-1.0)
---@param girder2 number Second girder ID
---@param span2_start number Start position on second girder (0.0-1.0)
---@param span2_end number End position on second girder (0.0-1.0)
---@return number shieldId ID of the created shield in C++
function Shield(girder1, span1_start, span1_end, girder2, span2_start, span2_end) end

--- Creates a tank component
---@param position number[] Position vector {x,y,z}
---@param direction number[] Direction vector {x,y,z}
---@param dimensions number[] Size vector {x,y,z}
---@param material string Material name
---@return number tankId ID of the created tank in C++
function Tank(position, direction, dimensions, material) end

--- Creates a thruster component
---@param position number[] Position vector {x,y,z}
---@param direction number[] Direction vector {x,y,z}
---@param dimensions number[] Size vector {x,y,z}
---@param kind string Thruster kind/type
---@return number thrusterId ID of the created thruster in C++
function Thruster(position, direction, dimensions, kind) end

--- Creates a gun component
---@param support number Support component ID
---@param span_start number Start position on support (0.0-1.0)
---@param span_end number End position on support (0.0-1.0)
---@param kind string Gun kind/type
---@return number gunId ID of the created gun in C++
function Gun(support, span_start, span_end, kind) end

--- Creates a rock component
---@param position number[] Position vector {x,y,z}
---@param direction number[] Direction vector {x,y,z}
---@param up number[] Up vector {x,y,z}
---@param dimensions number[] Size vector {x,y,z}
---@return number rockId ID of the created rock in C++
function Rock(position, direction, up, dimensions) end

--- Creates a balloon component
---@param position number[] Position vector {x,y,z}
---@param direction number[] Direction vector {x,y,z}
---@param up number[] Up vector {x,y,z}
---@param dimensions number[] Size vector {x,y,z}
---@return number balloonId ID of the created balloon in C++
function Balloon(position, direction, up, dimensions) end


ComponetKind = {Node=0,ShipComponent=1,StructuralComponent=2,NodeLinker=3,Girder=4,Ring=5,Rope=6,Pipe=7,Plate=8,Radiator=9,Shield=10,Collector=11,Thruster=12,Rotor=13,Slider=14,Accelerator=15,Gun=16,Modul=17,Tank=18,Balloon=19,Rock=20};

Material{ name="Kevlar", density=1.44e+3, Spull=3.6e+9, Spush=0.0,    Kpull=154.0e+9, Kpush=0.0,      reflectivity=0.6,  Tmelt=350 }
Material{ name="Steel" , density=7.89e+3, Spull=1.2e+9, Spush=1.2e+9, Kpull=200.0e+9, Kpush=200.0e+9, reflectivity=0.85, Tmelt=800 }
--Material{ name="Titanium" , density=7.89e+3, Spull=3.6e+9, Spush=0.0, Kpull=154.0e+9, Kpush=0.0, reflectivity=0.7, Tmelt=450 }

origin = {0.0,0.0,0.0}
xvec   = {1.0,0.0,0.0}
yvec   = {0.0,1.0,0.0}
zvec   = {0.0,0.0,1.0}

function tanks( n,aOff, R, r, L, z0 )
    for i=1,n do
        local a = (i/n + aOff)*2*math.pi 
        Tank( {math.cos(a)*R,math.sin(a)*R,z0}, zvec, {r,r,L}, "H2")
    end
end

function girderFan( n, aOff, w, h, z0,z1,   nseg, thick, shielded, st )
    local tipNode = Node({0.0,0.0,z1});
    local a = aOff*2*math.pi 
    local gd0  = Girder( tipNode, Node({math.cos(a)*w,math.sin(a)*h,z0}), xvec, nseg, 2, {thick,thick}, "steel", st )
    local ogd  = gd0;
    for i=1,n do
        local gd=gd0
        if i<n then
            local a = (i/n + aOff)*2*math.pi
            gd = Girder( tipNode, Node({math.cos(a)*w,math.sin(a)*h,z0}), xvec, nseg, 2, {thick,thick},  "steel", st )
        end
        if shielded then Shield( ogd,0.0,1.0, gd,0.0,1.0 ) end
        ogd = gd;
    end
end
