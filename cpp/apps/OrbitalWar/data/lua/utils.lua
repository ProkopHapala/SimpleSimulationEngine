Material{ name="Kevlar", density=1.44e+3, Spull=3.6e+9, Spush=0.0, Kpull=154.0e+9, Kpush=0.0, reflectivity=0.6,  Tmelt=350 }
Material{ name="Steel" , density=7.89e+3, Spull=1.2e+9, Spush=0.0, Kpull=200.0e+9, Kpush=0.0, reflectivity=0.85, Tmelt=800 }
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

function girderFan( n, aOff, w, h, z0,z1,   nseg, thick, shielded )
    local tipNode = Node({0.0,0.0,z1});
    local a = aOff*2*math.pi 
    local gd0  = Girder( tipNode, Node({math.cos(a)*w,math.sin(a)*h,z0}), xvec, nseg, 2, {thick,thick}, "steel")
    local ogd  = gd0;
    for i=1,n do
        local gd=gd0
        if i<n then
            local a = (i/n + aOff)*2*math.pi
            gd = Girder( tipNode, Node({math.cos(a)*w,math.sin(a)*h,z0}), xvec, nseg, 2, {thick,thick},  "steel" )
        end
        if shielded then Shield( ogd,0.0,1.0, gd,0.0,1.0 ) end
        ogd = gd;
    end
end
