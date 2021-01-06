#!/usr/bin/python
import numpy as np

const_airDensGround = 1.22;
const_GravAccel     = 9.81;

def trj( ts, omega, R=1, x0=0, y0=0 ):
    phi = ts*omega
    xs = R*np.sin(phi)          + x0
    ys = R*(1-np.cos(phi))      + y0
    return xs,ys, phi

def trunRate_simple( area_wing, CD0, aspect, power, mass,     area_hull=0, CLmax=1.5 ):
    area_tot = area_wing + area_hull;
    c3       = area_wing/(np.pi*aspect);
    c4       = CD0 * area_tot;
    CL       = np.sqrt(3*c4/c3);   # Best Lift
    CL       = np.minimum( CL, CLmax )
    P_Fd = power/( 0.5*const_airDensGround * ( CL*CL*c3 + c4 ) )
    speed  = np.power( P_Fd,  1./3. );
    Flift  = speed*speed * CL * area_wing * 0.5 * const_airDensGround
    accel  = Flift/mass
    #printf( "CL %g c3 %g c4 %g speed %g P_Fd %g Fl %g[kN] accel %g[g] \n", CL, c3, c4, speed, P_Fd, Flift*1e-3, accel/const_GravAccel );
    Rturn  = (speed*speed)/accel;
    omega = speed/Rturn
    return omega, Rturn  # Omega

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    
    o1, R1 = trunRate_simple( area_wing=16, CD0=0.03, aspect=6.0, power=1000e+3, mass=3400 ) # Bf109 G6
    o2, R2 = trunRate_simple( area_wing=22, CD0=0.03, aspect=6.0, power=500e+3,  mass=1960 ) # I-153 (M76)
    
    print R1," [m] ",o1," [rad/s] ", o1*R1*3.6," [km/h] "
    print R2," [m] ",o2," [rad/s] ", o2*R2*3.6," [km/h] "
    
    ts = np.arange(0,10.,0.1)
    x1,y1,a1 = trj( ts, 0.1,  R=500, x0=0  , y0=0 )
    x2,y2,a2 = trj( ts, 0.15, R=200, x0=300, y0=0 )
    
    tga = np.tan(a1)
    tgb = (y2-y1)/(x2-x1)
    
    plt.figure()
    plt.plot(x1,y1, label='hunter' )
    plt.plot(x2,y2, label='pray' )
    plt.axis('equal')
    plt.legend()
    plt.grid()
    
    plt.figure()
    plt.plot(ts,tga, label='hunter' )
    plt.plot(ts,tgb, label='pray' )
    plt.ylim(-3.,3.)
    plt.legend()
    plt.grid()
    
    plt.show()
