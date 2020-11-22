
import numpy as np

'''

Speed / Climb / Dive:

SimpleModel:
 - Neglect Lift-Induced drag
 - Lift automatically compensate perpendicular component of grafity
 - Thrust along flight direction
Fdrag  =  Fthrust + Fgrav*sin(a)
Cd*v^2 =  Fthrust + Fgrav*sin(a)
Cd*v^2 =  Fthrust + m*g*sin(a)
v = sqrt( ( Fthrust - m*g*sin(a) )/Cd )



Lift = LL * sa
Drag = DD * sa^2




Fx = cos(a)
Fy = sin(a)

'''

# Bf109 as reference : https://en.wikipedia.org/wiki/Messerschmitt_Bf_109
def climbModel1( angle, area=16, power=1085e+3, mass=3000, Cdrag = 0.01, v0=100 ):
    Fthrust = power/v0      ; print "Thrust[N]: ", Fthrust  
    Fg      = 9.81 * mass
    sa      = np.sin(angle)
    v = np.sqrt( (Fthrust - sa * Fg)/(0.5*Cdrag*1.2*area) )
    return v

''''
# Bf109 as reference : https://en.wikipedia.org/wiki/Messerschmitt_Bf_109
def climbModel2( angle, area=16, power=1085e+3, mass=3000, Cdrag = 0.01, v0=100 ):
    Fthrust = power/v0      ; print "Thrust[N]: ", Fthrust  
    Fg      = 9.81 * mass
    sa      = np.sin(angle)
    v = np.sqrt( (Fthrust - sa * Fg)/(0.5*Cdrag*1.2*area) )
    return v
'''

if __name__ == "__main__":  
    import matplotlib.pyplot as plt
    angles = np.linspace(-np.pi*0.5,np.pi*0.5,100)
    v      = climbModel1( angles, area=16, power=1085e+3, mass=3000, Cdrag = 0.01, v0=100 )
    #plt.plot(angles, v ); plt.xlabel("angle[rad]"); plt.ylabel("v[m/s]")
    ca = np.cos(angles)
    sa = np.sin(angles)
    plt.plot( ca*v, sa*v ); plt.xlabel("vx[m/s]"); plt.ylabel("vy[m/s]"); plt.axis('equal')
    
    plt.grid()
    plt.minorticks_on()
    plt.grid(which='minor', linestyle=':', linewidth='0.5', color='gray')
    plt.show()