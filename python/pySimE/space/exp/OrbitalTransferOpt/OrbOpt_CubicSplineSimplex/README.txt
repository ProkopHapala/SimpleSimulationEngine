

   Motivation
  ------------
  
	Purpose of this program is to provide very simple and illustrative algorithm to optimize orbital manuever with low thrust hi-specific impulse engine, such as ion-thuster or fission-framgnet rocket or solar sail. It can be also used for optimization of intercepting trajectory of one spaceship by an other. 
	
   Theory 
  --------
  
	Mathematicaly is this problem formulated as 4th-order boundary problem of second order non-linear diferential equation:
		r''(t) = T(t) + G(r(t))       Eq.(1)
		where r(t) is position vector of the spaceship, T(t) is time dependent vector of engine thrust, G(r(t)) is gravitational accleration due to massive central body (Sun, Earth ... ) at position r(t)
	
	At time t=0 is the ship in point r0 with velocity v0. After time t1 we want be at point r1 (position of target at time t1 ) with velocity v1. Note that target is also moving (for example a planet orbiting the Sun). If we are not interested just in fast fly-bay around our target, but we want to match it orbit, we have to bound the trajectory also by tho speed at and v1.
	
	So we have 4 boundary conditions
		r(0)   = r0
		r(t1)  = r1
		r'(0)  = v(0) = v0
		r'(t1) = v(t1) = v1
		
	Now we just choose some suitable description of curve r(t) which fulfill our boundary conditions. Aside exact fullfilment of boundary conditions we alo want:
		1) r(t), r'(t) be smooth
		2) simple analytical derivatives r'(t) r''(t)
		3) low number of parameters to optimize
		4) low entanglement of parameters ( this means that change of one parameter doesn't change much optimal value of other parameter  ) 
    	5) fast to numericaly evaluate ( including derivatives  r'(t), r''(t) )
	In this program we use clamped cubic spline for simplicity which fulfill these requrements quite well. ( "Clampled spline" mean that exact derivatives at ends are given )
	
	position vector r(t) can be expressed in cartesin or in polar coordinates. We choose polar because gravitational acceleration G(r(t)) can be easily computed in polar coordinates and because trajectories in polar coordinates are slowly varying curves ( less waves and fetures ) which lower the number of control points requred for it's description.
	
	Drawback of polar coordinate system (which is non-uniform and non-inertial ) is that expressin of kinematic accelerations ( r''(t) ) is a bit more complicated (however still easy to evaluate numerically ). Lets define unit vector in radial and angulr direction uR and uPhi.
		uR   =  (  cos( Phi ), sin( Phi ) )
		uPhi =  ( -sin( Phi ), cos( Phi ) )
		r(t)    =  R(t)*uR
	then the kinemetic acceleration can be expressed as
		r''(t) =   (  R''(t) -  R(t) * Phi'(t)  ) * uR   +   (  R(t) * Phi''(t) +  2 * R'(t) * Phi'(t)  )* uPhi
	gravitational acceleration is just in radial direction
		G(R(t)) = 1/R(t)^2
	Then the diferential equation Eq.(1) can be rewriten as 
		TPhi(t)  = R(t) * Phi''(t) +  2 * R'(t) * Phi'(t)
		TR  (t)  = R''(t)          -       R(t) * Phi'(t)  - 1/R(t)^2
	Where TPhi(t) and TR(t) is engine thrust which the ship must exert to move along this trajectory. We want to minimize total thrust
		| T(t) | = sqrt( TPhi(t)^2 + TR(t)^2 )	
	
	To do this we vary parameters (control points) of cubic splines R(t) and Phi(t). We have choosen the Simplex minimization method (Nelder–Mead method) which is very robust and reasonably fast. Even simpler choice could some sort of pattern search or stochastic minimization or even gentic algorithms or particle swarm optimization - Basically all methos which doesn't need gradietns.
	
	Because variational derivatives of | T(t) | acording to spline parameters are not so easy to evaluate analytically
	
	In each step the fittness function Integral(| T(t) |) is evaluated numerically at some sampling points in between control points. The sampled values are than integrated by trapezoid rule ( more sophisticated integration rules such as Simpson-rule or gauss-quadrature doesn't make sense there because |T(t)| is not smooth enought so it is better to increase number of sampling poins and keep order of integration method low. )  
	
	Because we use very simple minimization method (without gradietns ) we can alternatively formulate more sophisticated fittness function, which expres consumption model of our engine ( for example some polynom of |T(t)| ) or mecha our other objectives. 
	For example - when intercepting other ship we are not sure about it's manuvering in future ( = unpredicatable changes of position and velocity ). Because of this we would prefer spend less propelant at the beginning of the trajectory and keep more propelat for more agresive manuevers at the end. In this case the objective function can be modified by some time dependent weight  FF(t) = W(t)*|T(t)|

   Limitations and To Do
  -----------------------
  
  Cubic spline is not smooth enought accleration are not smooth at controlpoints ( it derivatives a'(t) = r'''(t) is not continuous here ).
  Even though it is physically possible to change engine thrust even sepwise, we can expect that optimal trajectory would be continuous to higher degree. Because of lower smoothens of the trajectory the objective function | T(t) | has sharp spikes which for sure doesnt coresponf to ideal trajectory.
  
  Higher number of controlpoinst cause the curves to be kinky, which makes parameters so entangled that even robust method like (Nelder–Mead method  or random wakl minimization) has hard time to converge.

	This can be possibly overcome by alternative representation of contropoints in input vector. For example hierarchical one where parametres such that:
	x( ti + dt/2 ) = ( x( ti ) + x( ti + dt )  )  +  dx( ti + dt/2 )  
	Where x(ti) are higher level of hirerarchy, and dx at halfstep are lower level of hierarchy.
	Because dx makes the curve kinky ( which affect the fittness strongly ), it can be rescaled in the parametr space correspondingly in order to make the steep and curvature of fittnes surface more uniform.  
	
	|___x___|   - top level   x
	|_x_|_x_|   - 2nd level  dx  ( scale x 4  )
	|x|x|x|x|   - 3rd level ddx  ( scale x 16 )
	
	
http://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method
http://en.wikipedia.org/wiki/Pattern_search_(optimization)
  http://en.wikipedia.org/wiki/Polar_coordinate_system
  http://en.wikipedia.org/wiki/Ion_thruster
  http://en.wikipedia.org/wiki/Fission-fragment_rocket
  http://en.wikipedia.org/wiki/Trajectory_optimization
  http://en.wikipedia.org/wiki/Orbital_maneuver
